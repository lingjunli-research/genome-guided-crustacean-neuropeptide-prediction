# -*- coding: utf-8 -*-
"""
Created on Sat Oct 11 13:20:46 2025

@author: lafields2
"""

import sys
from Bio import SeqIO
from rapidfuzz import fuzz, process
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

predicted_sequences = r"D:\Manuscripts\2025_denovo_sequencing\EG_searches\round2\output\Brain_TR1\Brain_TR1\final_results_EG_score_Brain_TR1_round2.csv"
output_dir = r"D:\Manuscripts\2025_denovo_sequencing\EG_searches\round2\output\Brain_TR1\Brain_TR1"

motif_db_path = r"D:\Manuscripts\2025_inSilicoPrediction\motif_DB.fasta"
base = Path("D:/Manuscripts")
np_seqs_path = base / "2025_inSilicoPrediction" / "np_database.fasta"
new_np_seqs_path = r"D:\Manuscripts\2025_inSilicoPrediction\new_NPs_since_2010.fasta"
# np_seqs_path = r"D:\Manuscripts\2025_inSilicoPrediction\np_database.fasta"


predicted = pd.read_csv(predicted_sequences)
predicted_list = predicted['Unmodified sequence'].values.tolist()

fasta_peptides = []
for record in SeqIO.parse(np_seqs_path, "fasta"):
    seq = str(record.seq).strip()
    if seq:
        fasta_peptides.append(seq)


results = []
for query in predicted_list:
    match, score, _ = process.extractOne(query, fasta_peptides, scorer=fuzz.ratio)
    results.append((query, match, score))

df = pd.DataFrame(results, columns=["DeNovo_Peptide", "Matched_Database_Seq", "Similarity_Score"])
df.to_csv(f"{output_dir}\\cNP_db_matches.csv", index=False)

plt.figure(figsize=(8, 5))
plt.hist(df["Similarity_Score"], bins=np.arange(0, 101, 2), color="#2A9D8F", edgecolor="black")
plt.title("Distribution of Similarity Scores")
plt.xlabel("Similarity Score")
plt.ylabel("Frequency")
plt.xlim(0, 100)
plt.tight_layout()

plt.savefig(f"{output_dir}\\cNP_match_similarity_score_histogram.svg", format="svg")
plt.show()

fasta_peptides = []
for record in SeqIO.parse(new_np_seqs_path, "fasta"):
    seq = str(record.seq).strip()
    if seq:
        fasta_peptides.append(seq)

results = []
for query in predicted_list:
    match, score, _ = process.extractOne(query, fasta_peptides, scorer=fuzz.ratio)
    results.append((query, match, score))

df = pd.DataFrame(results, columns=["DeNovo_Peptide", "Matched_Database_Seq", "Similarity_Score"])
df.to_csv(f"{output_dir}\\NPs_since_2010_db_matches.csv", index=False)

plt.figure(figsize=(8, 5))
plt.hist(df["Similarity_Score"], bins=np.arange(0, 101, 2), color="#2A9D8F", edgecolor="black")
plt.title("Distribution of Similarity Scores")
plt.xlabel("Similarity Score")
plt.ylabel("Frequency")
plt.xlim(0, 100)
plt.tight_layout()

plt.savefig(f"{output_dir}\\since2010_match_similarity_score_histogram.svg", format="svg")
plt.show()

fasta_peptides = []
for record in SeqIO.parse(motif_db_path, "fasta"):
    seq = str(record.seq).strip()
    if seq:
        fasta_peptides.append(seq)


results = []
for query in predicted_list:
    match, score, _ = process.extractOne(query, fasta_peptides, scorer=fuzz.ratio)
    results.append((query, match, score))

df = pd.DataFrame(results, columns=["DeNovo_Peptide", "Matched_Database_Seq", "Similarity_Score"])
df.to_csv(f"{output_dir}\\motif_db_matches.csv", index=False)

plt.figure(figsize=(8, 5))
plt.hist(df["Similarity_Score"], bins=np.arange(0, 101, 2), color="#2A9D8F", edgecolor="black")
plt.title("Distribution of Similarity Scores")
plt.xlabel("Similarity Score")
plt.ylabel("Frequency")
plt.xlim(0, 100)
plt.tight_layout()

plt.savefig(f"{output_dir}\\motifDB_similarity_score_histogram.svg", format="svg")
plt.show()