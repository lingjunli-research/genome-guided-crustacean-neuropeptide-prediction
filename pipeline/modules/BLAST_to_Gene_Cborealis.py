# -*- coding: utf-8 -*-
"""
Created on Wed Aug 20 12:28:45 2025

@author: lafields2
"""

import pandas as pd
import gzip

# === CONFIG ===
blast_hits_path = r"D:\Manuscripts\2025_inSilicoPrediction\Cancer borealis prediction\blast_hits_with_scaffolds_filtered.csv"
gff3_path       = r"D:\Manuscripts\2025_inSilicoPrediction\Cancer borealis prediction\JC_geneModels_braker3+manual.gff"  # or .gff.gz
output_path     = r"D:\Manuscripts\2025_inSilicoPrediction\Cancer borealis prediction\cbor_blast_gene_matches_formatted.csv"

# === Load and preprocess BLAST hits ===
blast_hits = pd.read_csv(blast_hits_path)
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

matches_df.to_csv(output_path, index=False)
print(f"✔ Done! Gene matches written to {output_path}")
