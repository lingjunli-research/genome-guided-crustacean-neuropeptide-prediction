# -*- coding: utf-8 -*-
"""
Created on Wed Sep 17 11:48:49 2025

@author: lafields2
"""

from pathlib import Path
import sys
import re
import pandas as pd
from pandas.errors import EmptyDataError
# ---- config you can tweak ----
ALLOWED_TAXIDS = {6763}   # still supported if your CSVs have staxid
#ACCESSION_PREFIXES = ("CM", "JAHF")
ACCESSION_PREFIXES = ("CM", "JBCDW")
DEFAULT_TITLE_TOKENS = {"Callinectes sapidus", "blue crab"}  # used only if a title/organism column exists

# ---------- helpers ----------
def smart_read_csv(path: Path) -> pd.DataFrame:
    # Try headered first
    try:
        df = pd.read_csv(path)
        # Heuristic: if the first "header" looks like data, fall back to header=None
        if df.columns[0].startswith(("sp|", "tr|")):
            raise ValueError("Likely headerless; retrying with header=None")
        return df
    except Exception:
        pass
    # Headerless (BLAST outfmt 10 style)
    return pd.read_csv(path, header=None)

def detect_cols(df: pd.DataFrame):
    """Return a dict of canonical fields -> column label (str or int)."""
    cols_lower = {c.lower(): c for c in df.columns if isinstance(c, str)}

    def pick(*names):
        for n in names:
            if n in cols_lower:
                return cols_lower[n]
        return None

    colmap = {
        "taxid":     pick("staxid", "taxid", "subject tax id", "subject_tax_id"),
        "organism":  pick("organism", "scientific name", "sci name", "sci_name", "sorganism"),
        "title":     pick("stitle", "salltitles", "subject title", "definition", "title"),
        "accession": pick("saccver", "sacc", "subject accession", "accession", "sseqid", "subject id"),
    }
    if any(colmap.values()):
        return colmap, False  # headered

    # Headerless BLAST outfmt 10: qseqid, saccver, pident, ...
    ncol = df.shape[1]
    return {"taxid": None, "organism": None, "title": None, "accession": 1 if ncol > 1 else None}, True

def norm(x):
    return "" if pd.isna(x) else str(x).strip()

# ---------- core ----------
def run_filter(family,
               in_dir: Path,
               out_dir: Path,
               
               accession_prefixes=ACCESSION_PREFIXES,
               ):

    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "filtered").mkdir(exist_ok=True)
    (out_dir / "reports").mkdir(exist_ok=True)

    summary, off_rows = [], []

    for csv_path in sorted(in_dir.glob("*.csv")):
        try:
            df = smart_read_csv(csv_path)
            colmap, headerless = detect_cols(df)
    
            keep = pd.Series(False, index=df.index)
    
            tax_col = colmap.get("taxid")
            org_col = colmap.get("organism")
            tit_col = colmap.get("title")
            acc_col = colmap.get("accession")
    
            # 1) If a taxid column exists, trust it
            if tax_col is not None:
                try:
                    keep |= (
                        df[tax_col].astype(str)
                        .str.extract(r"(\d+)", expand=False)
                        .astype("Int64")
                        .isin(list(ALLOWED_TAXIDS))
                        .fillna(False)
                    )
                except Exception:
                    pass
    
            # 2) If organism/title exist, allow “Callinectes sapidus” / “blue crab”
            for c in [org_col, tit_col]:
                if c is not None:
                    keep |= df[c].astype(str).str.lower().str.contains(
                        "|".join(x.lower() for x in DEFAULT_TITLE_TOKENS), regex=True, na=False
                    )
    
            # 3) Your rule: accession STARTS WITH CM or JAHF
            if acc_col is not None:
                keep |= df[acc_col].astype(str).str.upper().str.startswith(
                    tuple(p.upper() for p in accession_prefixes)
                )
    
            # Split & save
            df_csap = df[keep].copy()
            df_off  = df[~keep].copy()
    
            df_csap.to_csv(out_dir / "filtered" / csv_path.name, index=False)
    
            # Off-target report (best-effort fields)
            for _, r in df_off.iterrows():
                off_rows.append({
                    "source_file": csv_path.name,
                    "accession": r.get(acc_col) if acc_col in df.columns else None,
                    "organism":  r.get(org_col) if org_col in df.columns else None,
                    "title":     r.get(tit_col) if tit_col in df.columns else None,
                    "taxid":     r.get(tax_col) if tax_col in df.columns else None,
                })
    
            summary.append({
                "file": csv_path.name,
                "total_rows": int(len(df)),
                "csap_rows": int(len(df_csap)),
                "off_target_rows": int(len(df_off)),
                "headerless_detected": bool(headerless),
            })
        except EmptyDataError:
            pass

    pd.DataFrame(summary).to_csv(out_dir / "reports" / "summary.csv", index=False)
    pd.DataFrame(off_rows).to_csv(out_dir / "reports" / "off_target_hits.csv", index=False)
    print("Done.")
    print(f"- Filtered: {out_dir/'filtered'}")
    print(f"- Reports:  {out_dir/'reports'}")
    print(f"- # False IDs for {family}:  {len(off_rows)}")

# ---------- Spyder entry point ----------

import pathlib
def list_directories_pathlib(path):
    """Lists all directories within the given path using pathlib module."""
    p = pathlib.Path(path)
    directories = [item.name for item in p.iterdir() if item.is_dir()]
    return directories

#working_dir = r"D:\Manuscripts\2025_inSilicoPrediction\round2_FASTAs\Blue"
working_dir = r"D:\Manuscripts\2025_inSilicoPrediction\round2_FASTAs\Jonah"
working_dir_subfolders = list_directories_pathlib(working_dir)
print(working_dir_subfolders)


if __name__ == "__main__":
    for x in working_dir_subfolders:
    
        working_family_dir = f"{working_dir}\{x}"
        input_dir  = Path(working_family_dir)
        output_dir = Path(working_family_dir)
    
        run_filter(
            x,
            in_dir=input_dir,
            out_dir=output_dir,
            #accession_prefixes=("CM", "JAHF"),
            accession_prefixes=("CM", "JBCDW")
            )
