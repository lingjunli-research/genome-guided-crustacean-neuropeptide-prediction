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

predicted_sequences = r"D:\Manuscripts\2025_inSilicoPrediction\Manuscript\Motifs\Csap_full_motifs.csv"
predicted = pd.read_csv(predicted_sequences)
predicted_list = predicted['motif'].values.tolist()

output_dir = r"D:\Manuscripts\2025_inSilicoPrediction\Manuscript\Motifs\fuzzy_reps"

motif_db_path = r"D:\Manuscripts\2025_inSilicoPrediction\Manuscript\Motifs\Manual_motif_DB.csv"
motifs = pd.read_csv(motif_db_path)
motif_list = motifs['Sequence'].values.tolist()

results = []
for query in predicted_list:
    match, score, _ = process.extractOne(query, motif_list, scorer=fuzz.ratio)
    results.append((query, match, score))

df = pd.DataFrame(results, columns=["Predicted_Peptide", "Matched_Motif_Seq", "Similarity_Score"])
df.to_csv(f"{output_dir}\\csap_cNP_db_matches.csv", index=False)

plt.figure(figsize=(8, 5))
plt.hist(df["Similarity_Score"], bins=np.arange(0, 101, 2), color="#2A9D8F", edgecolor="black")
plt.title("Distribution of Similarity Scores")
plt.xlabel("Similarity Score")
plt.ylabel("Frequency")
plt.xlim(0, 100)
plt.tight_layout()

plt.savefig(f"{output_dir}\\csap_cNP_match_similarity_score_histogram.svg", format="svg")
plt.show()

# fasta_peptides = []
# for record in SeqIO.parse(new_np_seqs_path, "fasta"):
#     seq = str(record.seq).strip()
#     if seq:
#         fasta_peptides.append(seq)

# results = []
# for query in predicted_list:
#     match, score, _ = process.extractOne(query, fasta_peptides, scorer=fuzz.ratio)
#     results.append((query, match, score))

# df = pd.DataFrame(results, columns=["DeNovo_Peptide", "Matched_Database_Seq", "Similarity_Score"])
# df.to_csv(f"{output_dir}\\NPs_since_2010_db_matches.csv", index=False)

# plt.figure(figsize=(8, 5))
# plt.hist(df["Similarity_Score"], bins=np.arange(0, 101, 2), color="#2A9D8F", edgecolor="black")
# plt.title("Distribution of Similarity Scores")
# plt.xlabel("Similarity Score")
# plt.ylabel("Frequency")
# plt.xlim(0, 100)
# plt.tight_layout()

# plt.savefig(f"{output_dir}\\since2010_match_similarity_score_histogram.svg", format="svg")
# plt.show()

# fasta_peptides = []
# for record in SeqIO.parse(motif_db_path, "fasta"):
#     seq = str(record.seq).strip()
#     if seq:
#         fasta_peptides.append(seq)


# results = []
# for query in predicted_list:
#     match, score, _ = process.extractOne(query, fasta_peptides, scorer=fuzz.ratio)
#     results.append((query, match, score))

# df = pd.DataFrame(results, columns=["DeNovo_Peptide", "Matched_Database_Seq", "Similarity_Score"])
# df.to_csv(f"{output_dir}\\motif_db_matches.csv", index=False)

# plt.figure(figsize=(8, 5))
# plt.hist(df["Similarity_Score"], bins=np.arange(0, 101, 2), color="#2A9D8F", edgecolor="black")
# plt.title("Distribution of Similarity Scores")
# plt.xlabel("Similarity Score")
# plt.ylabel("Frequency")
# plt.xlim(0, 100)
# plt.tight_layout()

# plt.savefig(f"{output_dir}\\motifDB_similarity_score_histogram.svg", format="svg")
# plt.show()