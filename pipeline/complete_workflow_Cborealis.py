# -*- coding: utf-8 -*-
"""
Created on Mon Sep  8 14:51:30 2025

@author: lafields2
"""

import pandas as pd
import re
import json
import gzip
import os
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict

alignment_results_path = r"D:\Manuscripts\2025_inSilicoPrediction\Final_Data_Analysis\Final_Cbor_Predictions\NP_query_retry_v02"
jsonl_path = r"C:\Users\lawashburn\Documents\CborealisGenome\cborealis_report\ncbi_dataset\data\GCA_041682235.1\sequence_report.jsonl"
gff3_path       = r"D:\Manuscripts\2025_inSilicoPrediction\Cancer borealis prediction\JC_geneModels_braker3+manual.gff"  # or .gff.gz
genome_fasta = r"D:\Manuscripts\2025_inSilicoPrediction\Cancer borealis prediction\JC-genome_v3-230128\JC-genome_v3-230128.fasta"


e_value_threshold = 6
identities = 20


out_dir = f"{alignment_results_path}\\results"
try:
    os.mkdir(out_dir)
    print(f"Directory '{out_dir}' created successfully.")
except FileExistsError:
    pass

blosum_out_low = f'{alignment_results_path}\\BLOSUM62_001.csv'
blosum_out_high = f'{alignment_results_path}\\BLOSUM62_6.csv'
pam_out_low = f'{alignment_results_path}\\PAM30_001.csv'
pam_out_high = f'{alignment_results_path}\\PAM30_6.csv'

default_cols = [
    "qseqid","sseqid","pident","length","mismatch","gapopen",
    "qstart","qend","sstart","send","evalue","bitscore",
    "positives","-","frame"
]

STD = ["qseqid","sseqid","pident","length","mismatch","gapopen",
       "qstart","qend","sstart","send","evalue","bitscore"]
STD_POS = STD + ["positives"]
STD_QS  = STD + ["qframe","sframe"]
STD_POS_QS = STD + ["positives","qframe","sframe"]

def read_blast_table(path):
    # try to read with sniffed sep, no header
    df = pd.read_csv(path, sep=None, engine="python", header=None, comment="#")
    n = df.shape[1]
    # Assign column names based on column count
    if n == 12:
        df.columns = STD
    elif n == 13:
        df.columns = STD_POS
    elif n == 14:
        # most common is std + qframe + sframe
        df.columns = STD_QS
    elif n == 15:
        df.columns = STD_POS_QS
    else:
        df.columns = [f"col_{i}" for i in range(n)]
        print(f"Warning: unexpected number of columns ({n}) in {path}")
    return df

blosum_high = read_blast_table(blosum_out_high)
pam_high    = read_blast_table(pam_out_high)
blosum_low = read_blast_table(blosum_out_low)
pam_low   = read_blast_table(pam_out_low)

merge_alignments = pd.concat([blosum_high, pam_high, blosum_low, pam_low], ignore_index=True)

# Clean/convert pident and evalue to numeric
# Handle % signs, comma decimals, stray text
def to_float_series(s):
    s = s.astype(str)
    s = s.str.strip()
    s = s.str.replace("%", "", regex=False)
    s = s.str.replace(",", ".", regex=False)  # 72,5 -> 72.5
    s = s.str.extract(r'([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)')[0]  # keep numeric token
    return pd.to_numeric(s, errors="coerce")

merge_alignments["pident"] = to_float_series(merge_alignments["pident"])
merge_alignments["evalue"] = to_float_series(merge_alignments["evalue"])

# DEBUG: inspect distribution
#print(merge_alignments[["pident","evalue"]].describe())

# Use a loose E-value while testing pident logic
# (Your current e_value_threshold = 0 will drop almost everything.)
#e_value_threshold = 1E-6  # <-- TEMP for debugging; set back later
identities = identities   # from your variable above

merge_alignments_filtered = merge_alignments[
    (merge_alignments["evalue"] <= e_value_threshold) &
    (merge_alignments["pident"] >= identities)
]


scaffold_extract_df = merge_alignments_filtered.copy()
# scaffold_extract_df = scaffold_extract_df.drop(columns=["pident","length","mismatch","gapopen",
#                                                         "qstart","qend","evalue","bitscore",
#                                                         "positives","-","frame"])
# scaffold_extract_df.rename(columns=scaffold_extract_df.iloc[0]).drop(scaffold_extract_df.index[0]).reset_index(drop=True) #remove df header
###########



default_cols = [
    "qseqid","sseqid","pident","length","mismatch","gapopen",
    "qstart","qend","sstart","send","evalue","bitscore"
]

acc2name = {}

# regex for accessions (RefSeq CM… or WGS JAHFWG… etc.)
acc_pat = re.compile(r"(CM\d{6,}\.\d+|[A-Z]{4,}\d{6,}\.\d+)")

# keys that often contain the *scaffold/sequence* name that matches your GFF seqid
name_keys = [
    "sequence_name", "sequenceName", "name", "seq_name", "seqName",
    "assigned_molecule", "assignedMolecule", "chromosome"
]

# keys that can contain accessions
acc_keys = [
    "accession", "genbank_accession", "genbankAccession",
    "refseq_accession", "refseqAccession"
]

def first_present(d, keys):
    for k in keys:
        if k in d and d[k]:
            return d[k]
    return None

with open(jsonl_path, "r", encoding="utf-8") as f:
    for line in f:
        if not line.strip():
            continue
        obj = json.loads(line)
        # collect candidate accessions on this line
        accs = set()
        # explicit fields
        for k in acc_keys:
            v = obj.get(k)
            if isinstance(v, str) and acc_pat.fullmatch(v):
                accs.add(v)
        # also scan all string values for accession-like tokens
        for v in obj.values():
            if isinstance(v, str):
                for m in acc_pat.findall(v):
                    accs.add(m)

        # pick a sequence/scaffold name
        seq_name = first_present(obj, name_keys)

        # As a fallback, try to find a HiC-style name anywhere in the object
        if not seq_name:
            for v in obj.values():
                if isinstance(v, str):
                    m = re.search(r"(HiC_scaffold_\d+)", v)
                    if m:
                        seq_name = m.group(1)
                        break

        # If still no name, optionally use "assigned molecule" formats like "chr1" etc.
        # (We already considered "assigned_molecule" above.)

        # Record mapping(s)
        if seq_name and accs:
            for a in accs:
                acc2name[a] = seq_name

if not acc2name:
    raise SystemExit("No accession→name mappings found in JSONL. Inspect the file and adjust key lists above.")

# --- 2) Read BLAST table ---
# if blast_has_header:
#     df = pd.read_csv(blast_path, sep=",")
# else:
#     # headerless: assign defaults (adjust if your file has more/less columns)
#     df = pd.read_csv(blast_path, sep=",", header=None)
#     if df.shape[1] < len(default_cols):
#         raise SystemExit(f"BLAST table has {df.shape[1]} cols; expected ≥ {len(default_cols)}. "
#                          "Edit 'default_cols' to match your file.")
#     df = df.iloc[:, :len(default_cols)]
#     df.columns = default_cols

# Try to find subject id / start / end in a robust way
lower = [c.lower() for c in scaffold_extract_df.columns]
def colname(cands):
    for c in cands:
        if c in lower:
            return scaffold_extract_df.columns[lower.index(c)]
    raise KeyError(f"Missing columns: tried {cands}")

col_sseqid = colname(["sseqid","subjectid","subject id","subject"])
col_sstart = colname(["sstart","s.start","s start"])
col_send   = colname(["send","s.end","s end"])
# Query id for the output (optional)
col_qseqid = colname(["qseqid","queryid","query id","query"]) if any(
    c in lower for c in ["qseqid","queryid","query id","query"]) else None

# --- 3) Apply mapping & normalize coordinates ---
scaffold_extract_df["scaffold"] = scaffold_extract_df[col_sseqid].map(acc2name).fillna(scaffold_extract_df[col_sseqid])  # if already a name, keep it
sstart = scaffold_extract_df[col_sstart].astype(int)
send   = scaffold_extract_df[col_send].astype(int)
scaffold_extract_df["start"] = sstart.where(sstart <= send, send)
scaffold_extract_df["end"]   = send.where(sstart <= send, sstart)

if col_qseqid:
    out = scaffold_extract_df[[col_qseqid, "scaffold", "start", "end"]].copy()
    out.columns = ["qseqid","scaffold","start","end"]
else:
    out = scaffold_extract_df[["scaffold","start","end"]].copy()
    out.insert(0, "qseqid", "(unknown)")

# Report any unmapped accessions
unmapped = scaffold_extract_df[scaffold_extract_df["scaffold"] == scaffold_extract_df[col_sseqid]]
if len(unmapped):
    print(f"{len(unmapped)} rows had SubjectIDs not found in the JSONL mapping.")
    print(unmapped[[col_sseqid]].drop_duplicates().head())
out_csv = f"{out_dir}\\blast_hits_with_scaffolds.csv"
out.to_csv(out_csv, index=False)
print(f"Wrote scaffold-based hits to: {out_csv}")

scaffolds = out.copy()

scaffolds_filtered = scaffolds[scaffolds['scaffold'].str.contains('HiC_scaffold')]

out_csv = f"{out_dir}\\blast_hits_with_scaffolds_filtered.csv"
scaffolds_filtered.to_csv(out_csv, index=False)

blast_hits = scaffolds_filtered.copy()
blast_hits['true_start'] = blast_hits[['start', 'end']].min(axis=1)
blast_hits['true_end']   = blast_hits[['start', 'end']].max(axis=1)

# --- helpers ---
def open_text(path):
    if path.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")

def parse_attrs(attr_str):
    d = {}
    for item in attr_str.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            d[k] = v
    return d

def norm_gene_id(s):
    # normalize 'g1234.t1' -> 'g1234'; 'g1234-RA' -> 'g1234'
    s = str(s).strip()
    # strip common transcript suffixes but keep base gene id
    for suf in (".t1", ".t2", ".t3"):
        if s.endswith(suf):
            return s[: -len(suf)]
    if s.endswith("-RA") or s.endswith("-RB") or s.endswith("-RC"):
        return s[:-3]
    return s

# === Load GFF (gene + transcript/mRNA) ===
records = []
wanted = {"gene", "mrna", "transcript"}  # case-insensitive check

with open_text(gff3_path) as f:
    for line in f:
        if not line or line.startswith("#"):
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue
        seqid, source, ftype, start, end, score, strand, phase, attrs = parts
        ftype_l = ftype.lower()
        if ftype_l not in wanted:
            continue
        start, end = int(start), int(end)
        A = parse_attrs(attrs)

        if ftype_l == "gene":
            gene_id = A.get("ID") or A.get("gene_id") or A.get("Name")
            if not gene_id:
                continue
            base_gene = norm_gene_id(gene_id)
            records.append([seqid, start, end, strand, base_gene, "gene"])
        else:  # transcript or mRNA
            # transcript ID is A.get("ID"); parent points to gene
            parent_gene = A.get("Parent") or A.get("gene_id") or A.get("gene")
            if not parent_gene:
                continue
            base_gene = norm_gene_id(parent_gene)
            # keep the transcript id for context (optional)
            tx_id = A.get("ID")
            records.append([seqid, start, end, strand, base_gene, "transcript", tx_id])

# unify to DataFrame
cols = ["scaffold", "start", "end", "strand", "gene_id", "feature_type", "transcript_id"]
gff_df = pd.DataFrame(records, columns=cols)
if "transcript_id" not in gff_df.columns:
    gff_df["transcript_id"] = ""

# === Match BLAST hits to overlapping gene/transcript features ===
matches = []
for _, hit in blast_hits.iterrows():
    scaffold = str(hit["scaffold"])
    hit_start, hit_end = int(hit["true_start"]), int(hit["true_end"])

    sub = gff_df[
        (gff_df["scaffold"] == scaffold) &
        (gff_df["start"]   <= hit_end) &
        (gff_df["end"]     >= hit_start)
    ]

    if sub.empty:
        matches.append([
            hit.get("qseqid", hit.get("query_id", "")),
            scaffold, hit_start, hit_end,
            "NO_MATCH", "", "", ""
        ])
        continue

    # Prefer 'gene' rows over 'transcript' if both overlap
    gene_rows = sub[sub["feature_type"] == "gene"]
    rows_to_use = gene_rows if not gene_rows.empty else sub

    for _, row in rows_to_use.iterrows():
        matches.append([
            hit.get("qseqid", hit.get("query_id", "")),
            scaffold, hit_start, hit_end,
            row["gene_id"], row["feature_type"],
            int(row["start"]), int(row["end"]),
        ])

# === Save ===
matches_df = pd.DataFrame(matches, columns=[
    "query_id", "scaffold", "hit_start", "hit_end",
    "matched_gene_id", "feature_type", "gene_start", "gene_end"
])
out_csv = f"{out_dir}\\gene_matches.csv"
matches_df.to_csv(out_csv, index=False)
print(f"Gene matches written to {out_csv}")

scaffolds_genes_filtered = matches_df[matches_df["feature_type"] == "gene"]
scaffolds_of_interest_dups = scaffolds_genes_filtered['scaffold'].values.tolist()
scaffolds_of_interest = set(scaffolds_of_interest_dups)

scaffold_output_dir = f"{out_dir}\\scaffolds_extracted"
# Ensure output directory exists
os.makedirs(scaffold_output_dir, exist_ok=True)

# Create a dictionary to hold open file handles
output_files = {
    scaffold: open(os.path.join(scaffold_output_dir, f"{scaffold}.gff"), 'w')
    for scaffold in scaffolds_of_interest
}

# Process input GFF
with open(gff3_path, 'r') as fin:
    for line in fin:
        if line.startswith('#') or line.strip() == '':
            for fout in output_files.values():
                fout.write(line)
        else:
            scaffold = line.split('\t')[0]
            if scaffold in scaffolds_of_interest:
                output_files[scaffold].write(line)

# Close all file handles
for fout in output_files.values():
    fout.close()

print(f"Wrote {len(scaffolds_of_interest)} GFF files to {scaffold_output_dir}")


output_fasta = f"{scaffold_output_dir}\\cbor_predicted_proteins.fa"

# ---- Load genome ----
genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))

# ---- Load target genes ----
hits = scaffolds_genes_filtered.copy()
hits = hits[(hits["feature_type"] == "gene") & (hits["matched_gene_id"] != "NO_MATCH")]

def norm_gene(s):
    s = str(s).strip()
    return re.split(r'(-RA|\.t\d+|\.\d+)$', s)[0]

target_genes = {norm_gene(x) for x in hits["matched_gene_id"].tolist()}

# ---- GFF helpers ----
def parse_attrs(attr):
    d = {}
    for item in attr.split(";"):
        if not item: continue
        if "=" in item:
            k,v = item.split("=",1)
            d[k] = v
    return d

def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if path.endswith(".gz") \
           else open(path, "r", encoding="utf-8", errors="replace")

# 1) Build transcript -> gene map (handles 'transcript' and 'mRNA')
t2g = {}
gff_lines = []
with open_text(gff3_path) as fh:
    for ln in fh:
        gff_lines.append(ln)
        if ln.startswith("#"): continue
        parts = ln.rstrip("\n").split("\t")
        if len(parts) < 9: continue
        ftype = parts[2].lower()
        if ftype in ("transcript","mrna"):
            attrs = parse_attrs(parts[8])
            tid = attrs.get("ID") or attrs.get("transcript_id")
            gid = attrs.get("Parent") or attrs.get("gene_id") or attrs.get("gene")
            if tid and gid:
                t2g[tid] = norm_gene(gid)

# 2) Collect CDS pieces keyed by gene, including phase (frame)
cds_by_gene = defaultdict(list)
for ln in gff_lines:
    if ln.startswith("#"): continue
    parts = ln.rstrip("\n").split("\t")
    if len(parts) < 9: continue
    seqid, source, ftype, start, end, score, strand, phase, attr = parts
    if ftype != "CDS": continue
    attrs = parse_attrs(attr)
    parent = (attrs.get("Parent") or attrs.get("transcript_id") or "").split(",")[0]
    gid = t2g.get(parent) or norm_gene(attrs.get("gene_id",""))
    if not gid or gid not in target_genes: 
        continue
    # phase is '0','1','2' or '.'
    ph = 0
    if phase in ("0","1","2"):
        ph = int(phase)
    cds_by_gene[gid].append({
        "chrom": seqid,
        "start": int(start),
        "end": int(end),
        "strand": strand,
        "phase": ph
    })

print(f"Found CDS for {len(cds_by_gene)} / {len(target_genes)} target genes.")

# 3) Build CDS with phase-correct trimming and translate
n_written = 0
with open(output_fasta, "w") as out_fa:
    for gid, pieces in cds_by_gene.items():
        if not pieces:
            continue
        strand = pieces[0]["strand"]
        if strand == "+":
            pieces.sort(key=lambda r: r["start"])
            seq_parts = []
            for p in pieces:
                s = p["start"] - 1
                e = p["end"]
                sub = genome[p["chrom"]].seq[s:e]
                # For + strand, trim 'phase' bases from the LEFT
                if p["phase"] in (1,2):
                    sub = sub[p["phase"]:]
                seq_parts.append(sub)
            cds = Seq("").join(seq_parts)
        else:
            # For - strand, process in reverse genomic order, trim 'phase' from the RIGHT of each piece
            pieces.sort(key=lambda r: r["end"], reverse=True)
            seq_parts = []
            for p in pieces:
                s = p["start"] - 1
                e = p["end"]
                sub = genome[p["chrom"]].seq[s:e]
                if p["phase"] in (1,2):
                    sub = sub[:len(sub) - p["phase"]]
                seq_parts.append(sub)
            cds = Seq("").join(seq_parts).reverse_complement()

        # safety: length multiple of 3 improves translation
        rem = len(cds) % 3
        if rem != 0:
            # keep in-frame; drop trailing extra bases
            cds = cds[:len(cds)-rem]

        if len(cds) == 0:
            continue

        # translate; don't require canonical start/stop
        prot = cds.translate(to_stop=True)
        if len(prot) == 0:
            # fallback: translate full to inspect internal stops
            prot = cds.translate()
            # you can comment the next line if you prefer to keep even with '*'
            prot = str(prot).split("*")[0]

        if len(prot) == 0:
            continue

        out_fa.write(f">{gid}\n{prot}\n")
        n_written += 1

print(f"Wrote {n_written} proteins to: {output_fasta}")
