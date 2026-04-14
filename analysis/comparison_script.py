# -*- coding: utf-8 -*-
"""
Created on Mon Sep  8 10:27:21 2025

@author: lafields2
"""

import pandas as pd
import re

# ===============================
# User inputs – paste your paths here
# ===============================
putative_file = r"D:\Manuscripts\2025_inSilicoPrediction\comparison_w_CNPDB\putative.txt"
verified_file = r"D:\Manuscripts\2025_inSilicoPrediction\comparison_w_CNPDB\verified.txt"
output_excel = r"D:\Manuscripts\2025_inSilicoPrediction\comparison_w_CNPDB\peptide_comparison.xlsx"

ignore_case = True     # True: case-insensitive comparison
il_equivalent = True   # True: treat I and L as equivalent (also maps J->L)
strip_nonletters = True  # True: remove non A–Z chars before comparing
# ===============================


def normalize_seq(seq: str, ignore_case=True, il_equivalent=True, strip_nonletters=True) -> str:
    """Normalize a peptide for comparison while preserving an original copy for reporting."""
    s = seq.strip()
    if strip_nonletters:
        s = re.sub(r"[^A-Za-z]", "", s)
    if ignore_case:
        s = s.upper()
    if il_equivalent:
        # Treat I/L as indistinguishable. Also coerce J (I/L ambiguous) to L.
        s = s.replace("I", "L").replace("J", "L")
        # If you prefer map both to 'I' instead, swap to: s = s.replace("L","I").replace("J","I")
    return s


def read_peptides(txt_path):
    with open(txt_path, "r", encoding="utf-8") as f:
        return [line.strip() for line in f if line.strip()]

# ---- Read lists ----
putative_raw = read_peptides(putative_file)
verified_raw = read_peptides(verified_file)

# ---- Normalize for comparison (keep originals for output) ----
put_norm = [normalize_seq(s, ignore_case, il_equivalent, strip_nonletters) for s in putative_raw]
ver_norm = [normalize_seq(s, ignore_case, il_equivalent, strip_nonletters) for s in verified_raw]

# Collapse to unique while preserving first original appearance
put_unique = list(dict.fromkeys(zip(put_norm, putative_raw)))  # [(norm, orig), ...]
ver_unique = list(dict.fromkeys(zip(ver_norm, verified_raw)))

# Lookups (normalized -> original)
put_lookup = {p_norm: p_orig for p_norm, p_orig in put_unique}
ver_lookup = {v_norm: v_orig for v_norm, v_orig in ver_unique}

put_norm_only = [p[0] for p in put_unique]
ver_norm_only  = [v[0] for v in ver_unique]

# ---- Exact matches on normalized strings ----
exact_norm = set(put_norm_only).intersection(ver_norm_only)
df_exact = pd.DataFrame(
    [{"Peptide": put_lookup[n], "Matched_in_verified": ver_lookup[n]} for n in sorted(exact_norm)],
    columns=["Peptide", "Matched_in_verified"]
)

# ---- Substring matches (exclude exact) ----
put_in_ver_rows = []   # Putative ⊂ Verified
for p_norm, p_orig in put_unique:
    for v_norm, v_orig in ver_unique:
        if p_norm == v_norm:
            continue
        if p_norm and p_norm in v_norm:
            put_in_ver_rows.append({"Putative": p_orig, "Verified_match": v_orig})

ver_in_put_rows = []   # Verified ⊂ Putative
for v_norm, v_orig in ver_unique:
    for p_norm, p_orig in put_unique:
        if v_norm == p_norm:
            continue
        if v_norm and v_norm in p_norm:
            ver_in_put_rows.append({"Verified": v_orig, "Putative_match": p_orig})

df_put_in_ver = pd.DataFrame(put_in_ver_rows, columns=["Putative", "Verified_match"])
df_ver_in_put = pd.DataFrame(ver_in_put_rows, columns=["Verified", "Putative_match"])

# ---- Determine "only" sets (no exact and no substring in either direction) ----
# Build sets of normalized sequences that have ANY relationship
related_put_norms = set(exact_norm)
related_ver_norms = set(exact_norm)

related_put_norms |= set(
    normalize_seq(row["Putative"], ignore_case, il_equivalent, strip_nonletters) for row in put_in_ver_rows
)
related_ver_norms |= set(
    normalize_seq(row["Verified_match"], ignore_case, il_equivalent, strip_nonletters) for row in put_in_ver_rows
)
related_ver_norms |= set(
    normalize_seq(row["Verified"], ignore_case, il_equivalent, strip_nonletters) for row in ver_in_put_rows
)
related_put_norms |= set(
    normalize_seq(row["Putative_match"], ignore_case, il_equivalent, strip_nonletters) for row in ver_in_put_rows
)

putative_only = [p_orig for p_norm, p_orig in put_unique if p_norm not in related_put_norms]
verified_only = [v_orig for v_norm, v_orig in ver_unique if v_norm not in related_ver_norms]

df_put_only = pd.DataFrame(putative_only, columns=["Putative_only"])
df_ver_only = pd.DataFrame(verified_only, columns=["Verified_only"])

# ---- Summary ----
df_summary = pd.DataFrame([{
    "Putative_count (unique)": len(put_unique),
    "Verified_count (unique)": len(ver_unique),
    "Exact_matches (I/L-eq)": len(df_exact),
    "Putative_in_Verified (substring)": len(df_put_in_ver),
    "Verified_in_Putative (substring)": len(df_ver_in_put),
    "Putative_only": len(df_put_only),
    "Verified_only": len(df_ver_only),
    "Case_insensitive": bool(ignore_case),
    "I/L_equivalent": bool(il_equivalent),
    "Stripped_nonletters": bool(strip_nonletters),
}])

# ---- Write Excel ----
with pd.ExcelWriter(output_excel) as writer:
    df_summary.to_excel(writer, sheet_name="Summary", index=False)
    df_exact.to_excel(writer, sheet_name="Exact_matches", index=False)
    df_put_in_ver.to_excel(writer, sheet_name="Putative_in_Verified", index=False)
    df_ver_in_put.to_excel(writer, sheet_name="Verified_in_Putative", index=False)
    df_put_only.to_excel(writer, sheet_name="Putative_only", index=False)
    df_ver_only.to_excel(writer, sheet_name="Verified_only", index=False)

print(f"Comparison complete. Results written to:\n{output_excel}")
