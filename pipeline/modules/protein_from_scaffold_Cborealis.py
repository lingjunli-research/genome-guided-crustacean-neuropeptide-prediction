from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import pandas as pd
import re
import os
import gzip

# --- INPUTS (adjust as before) ---
blast_gene_matches_path = r"D:\Manuscripts\2025_inSilicoPrediction\Cancer borealis prediction\cbor_blast_gene_matches_formatted_filtered.csv"
genome_fasta = r"D:\Manuscripts\2025_inSilicoPrediction\Cancer borealis prediction\JC-genome_v3-230128.fasta\JC-genome_v3-230128.fasta"
gff_path = r"D:\Manuscripts\2025_inSilicoPrediction\Cancer borealis prediction\JC_geneModels_braker3+manual.gff"   # use the single full GFF
output_fasta = r"D:\Manuscripts\2025_inSilicoPrediction\Cancer borealis prediction\cbor_predicted_proteins.fa"

# ---- Load genome ----
genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))

# ---- Load target genes ----
hits = pd.read_csv(blast_gene_matches_path)
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
with open_text(gff_path) as fh:
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

print(f"✅ Wrote {n_written} proteins to: {output_fasta}")
