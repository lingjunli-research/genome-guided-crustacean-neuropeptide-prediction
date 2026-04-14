# -*- coding: utf-8 -*-
"""
Created on Sat Oct 11 17:44:45 2025

@author: lafields2
"""

import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt

round2fasta_dir_CB = r"D:\Manuscripts\2025_inSilicoPrediction\Final_Data_Analysis\Final_Cbor_Predictions\Round2Fasta"
round2fasta_dir_CS = r"D:\Manuscripts\2025_inSilicoPrediction\Final_Data_Analysis\Final_Csap_Predictions\round2_FASTAs\Blue"
confident_score = 0.65

others_CB = [r"D:\Manuscripts\2025_inSilicoPrediction\Final_Data_Analysis\Final_Cbor_Predictions\HT_retry_v02\Signal.txt",
          r"D:\Manuscripts\2025_inSilicoPrediction\Final_Data_Analysis\Final_Cbor_Predictions\NP_query_retry_v02\SignalP_results.txt"]
others_CS = [r"D:\Manuscripts\2025_inSilicoPrediction\Final_Data_Analysis\Final_Csap_Predictions\Np_query_retry\Signal.txt",
          r"D:\Manuscripts\2025_inSilicoPrediction\Final_Data_Analysis\Final_Csap_Predictions\HT_retry\signal.txt"]
output_dir = r"D:\Manuscripts\2025_inSilicoPrediction\Manuscript"

def get_dir_names_with_strings_list(str_list,results_directory):
    full_list = [name for name in os.listdir(results_directory)
                 if os.path.isdir(os.path.join(results_directory, name))]
    return [nm for ps in str_list for nm in full_list if ps in nm]

str_list = ['']
results_dir_list_CB = get_dir_names_with_strings_list(str_list,round2fasta_dir_CB)

signal_df_CB = pd.DataFrame()

for x in results_dir_list_CB:
    signal_path_CB = f'{round2fasta_dir_CB}\\{x}\\SignalP.txt'
    try:
        signal_path_CB = f'{round2fasta_dir_CB}\\{x}\\SignalP.txt'
        signal_peps_CB = pd.read_csv(signal_path_CB, sep='\t',skiprows=1)
    except FileNotFoundError:
        signal_path_CB = f'{round2fasta_dir_CB}\\{x}\\Signal.txt'
        signal_peps_CB = pd.read_csv(signal_path_CB, sep='\t',skiprows=1)
    signal_peps_filtered_CB = signal_peps_CB[signal_peps_CB['Prediction'] == 'SP']
    signal_peps_filtered_CB = signal_peps_filtered_CB[signal_peps_filtered_CB['SP(Sec/SPI)'] >= confident_score]
    signal_df_CB = pd.concat([signal_df_CB, signal_peps_filtered_CB])
    
for y in others_CB:
    signal_peps_others_CB = pd.read_csv(signal_path_CB, sep='\t',skiprows=1)
    signal_peps_others_filtered_CB = signal_peps_others_CB[signal_peps_others_CB['Prediction'] == 'SP']
    signal_peps_others_filtered_CB = signal_peps_others_filtered_CB[signal_peps_others_filtered_CB['SP(Sec/SPI)'] >= confident_score]
    signal_df_CB = pd.concat([signal_df_CB, signal_peps_others_filtered_CB])

signal_df_nodups_CB = signal_df_CB.drop_duplicates()
signal_df_nodups_CB['SP_length'] = signal_df_nodups_CB['CS Position'].str.extract(r'pos:\s*(\d+)-')
s_CB = signal_df_nodups_CB['SP_length'].dropna().astype(int)

results_dir_list_CS = get_dir_names_with_strings_list(str_list,round2fasta_dir_CS)

signal_df_CS = pd.DataFrame()

for x in results_dir_list_CS:
    signal_path_CS = f'{round2fasta_dir_CS}\\{x}\\SignalP.txt'
    try:
        signal_path_CS = f'{round2fasta_dir_CS}\\{x}\\SignalP.txt'
        signal_peps_CS = pd.read_csv(signal_path_CS, sep='\t',skiprows=1)
    except FileNotFoundError:
        signal_path_CS = f'{round2fasta_dir_CS}\\{x}\\Signal.txt'
        signal_peps_CS = pd.read_csv(signal_path_CS, sep='\t',skiprows=1)
    signal_peps_filtered_CS = signal_peps_CS[signal_peps_CS['Prediction'] == 'SP']
    signal_peps_filtered_CS = signal_peps_filtered_CS[signal_peps_filtered_CS['SP(Sec/SPI)'] >= confident_score]
    signal_df_CS = pd.concat([signal_df_CS, signal_peps_filtered_CS])
    
for y in others_CS:
    signal_peps_others_CS = pd.read_csv(signal_path_CS, sep='\t',skiprows=1)
    signal_peps_others_filtered_CS = signal_peps_others_CS[signal_peps_others_CS['Prediction'] == 'SP']
    signal_peps_others_filtered_CS = signal_peps_others_filtered_CS[signal_peps_others_filtered_CS['SP(Sec/SPI)'] >= confident_score]
    signal_df_CS = pd.concat([signal_df_CS, signal_peps_others_filtered_CS])
signal_df_nodups_CS = signal_df_CS.drop_duplicates()
signal_df_nodups_CS['SP_length'] = signal_df_nodups_CS['CS Position'].str.extract(r'pos:\s*(\d+)-')
s_CS = signal_df_nodups_CS['SP_length'].dropna().astype(int)


def use_prism_style():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 12,
        "axes.linewidth": 1.2,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 5,
        "ytick.major.size": 5,
        "xtick.major.width": 1.2,
        "ytick.major.width": 1.2,
        "figure.dpi": 120,
    })

global_min = min(s_CB.min(), s_CS.min())
global_max = max(s_CB.max(), s_CS.max())
bins = np.arange(global_min - 0.5, global_max + 1.5, 1)

use_prism_style()
fig, ax = plt.subplots(figsize=(6.8, 4.2))

counts_list, edges, patches = ax.hist(
    [s_CB, s_CS],
    bins=bins,
    stacked=True,
    edgecolor="black",
    linewidth=2.0,
    alpha=0.5,  # 🔹 slightly transparent bars
    label=["C. borealis", "C. sapidus"]
)

ax.set_xlabel("Signal peptide length")
ax.set_ylabel("Frequency")
#ax.set_title("Start position distribution (stacked)", pad=8)
ax.grid(axis="y", linestyle="-", alpha=0.25)
ax.set_xlim(bins[0], bins[-1])
tick_step = max(1, (global_max - global_min) // 10 or 1)
ax.set_xticks(np.arange(global_min, global_max + 1, tick_step))
ax.legend(frameon=False)
#annotate_totals(ax, counts_list, edges)

fig.tight_layout()
plt.show()
fig.savefig(f"{output_dir}\\signal_peptide_distribution.svg",
            format="svg", dpi=300)