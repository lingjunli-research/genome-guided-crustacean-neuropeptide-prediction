# -*- coding: utf-8 -*-
"""
Created on Wed Aug 13 09:39:41 2025

@author: lafields2
"""

import re
import json
import pandas as pd
from pathlib import Path

blast_path = r"D:\Manuscripts\2025_inSilicoPrediction\Cancer borealis prediction\HitsFiltered.csv"
jsonl_path = r"C:\Users\lawashburn\Documents\CborealisGenome\cborealis_report\ncbi_dataset\data\GCA_041682235.1\sequence_report.jsonl"
blast_has_header = False
out_csv = r"D:\Manuscripts\2025_inSilicoPrediction\Cancer borealis prediction\blast_hits_with_scaffolds.csv"

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

# Quick sanity check
print(f"Built mapping for {len(acc2name)} accessions. Example:",
      next(iter(acc2name.items())) if acc2name else "(none)")

if not acc2name:
    raise SystemExit("No accession→name mappings found in JSONL. Inspect the file and adjust key lists above.")

# --- 2) Read BLAST table ---
if blast_has_header:
    df = pd.read_csv(blast_path, sep=",")
else:
    # headerless: assign defaults (adjust if your file has more/less columns)
    df = pd.read_csv(blast_path, sep=",", header=None)
    if df.shape[1] < len(default_cols):
        raise SystemExit(f"BLAST table has {df.shape[1]} cols; expected ≥ {len(default_cols)}. "
                         "Edit 'default_cols' to match your file.")
    df = df.iloc[:, :len(default_cols)]
    df.columns = default_cols

# Try to find subject id / start / end in a robust way
lower = [c.lower() for c in df.columns]
def colname(cands):
    for c in cands:
        if c in lower:
            return df.columns[lower.index(c)]
    raise KeyError(f"Missing columns: tried {cands}")

col_sseqid = colname(["sseqid","subjectid","subject id","subject"])
col_sstart = colname(["sstart","s.start","s start"])
col_send   = colname(["send","s.end","s end"])
# Query id for the output (optional)
col_qseqid = colname(["qseqid","queryid","query id","query"]) if any(
    c in lower for c in ["qseqid","queryid","query id","query"]) else None

# --- 3) Apply mapping & normalize coordinates ---
df["scaffold"] = df[col_sseqid].map(acc2name).fillna(df[col_sseqid])  # if already a name, keep it
sstart = df[col_sstart].astype(int)
send   = df[col_send].astype(int)
df["start"] = sstart.where(sstart <= send, send)
df["end"]   = send.where(sstart <= send, sstart)

# Build simple output
if col_qseqid:
    out = df[[col_qseqid, "scaffold", "start", "end"]].copy()
    out.columns = ["qseqid","scaffold","start","end"]
else:
    out = df[["scaffold","start","end"]].copy()
    out.insert(0, "qseqid", "(unknown)")

# Report any unmapped accessions
unmapped = df[df["scaffold"] == df[col_sseqid]]
if len(unmapped):
    print(f"⚠ {len(unmapped)} rows had SubjectIDs not found in the JSONL mapping.")
    print(unmapped[[col_sseqid]].drop_duplicates().head())

out.to_csv(out_csv, index=False)
print(f"✅ Wrote scaffold-based hits to: {out_csv}")