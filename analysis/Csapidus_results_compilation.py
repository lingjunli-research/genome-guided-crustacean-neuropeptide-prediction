# -*- coding: utf-8 -*-
"""
Created on Wed Oct  8 15:42:25 2025

@author: lafields2
"""

import pandas as pd
import os

results_folder = r"D:\Manuscripts\2025_inSilicoPrediction\Final_Data_Analysis\EG_validation\csap_predictions_search\output_data"
neuropred_results_path = r"D:\Manuscripts\2025_inSilicoPrediction\Final_Data_Analysis\Final_Csap_Predictions\full_working_prediction_database.csv"
neuropred_results_pt2_path = r"D:\Manuscripts\2025_inSilicoPrediction\Final_Data_Analysis\Final_Csap_Predictions\putative Blue_combined_round2FASTAs.xlsx"

merged_gene_paths = [
    r"D:\Manuscripts\2025_inSilicoPrediction\Final_Data_Analysis\Final_Csap_Predictions\HT_retry\results\gene_matches.csv",
    r"D:\Manuscripts\2025_inSilicoPrediction\Final_Data_Analysis\Final_Csap_Predictions\Np_query_retry\results\gene_matches.csv"
                     ]

outside_gene_dir = r"D:\Manuscripts\2025_inSilicoPrediction\Final_Data_Analysis\Final_Csap_Predictions\round2_FASTAs\Blue"


def get_dir_names_with_strings_list(str_list,results_directory):
    full_list = [name for name in os.listdir(results_directory)
                 if os.path.isdir(os.path.join(results_directory, name))]
    return [nm for ps in str_list for nm in full_list if ps in nm]

str_list = ['']
results_dir_list = get_dir_names_with_strings_list(str_list,results_folder)
merged_df = pd.DataFrame()

for folder in results_dir_list:
    folder_path = f'{results_folder}\\{folder}'
    subdir_list = get_dir_names_with_strings_list(str_list,folder_path)
    subdir = subdir_list[0]
    results_path = f'{folder_path}\\{subdir}\\final_results_EG_score.csv'
    results = pd.read_csv(results_path)
    results_df = pd.DataFrame({'Peptide': results['Peptide'],
                               f'{subdir}_intensity': results['precursor_intensity']})
    if len(merged_df) == 0:
        merged_df = results_df
    else:
        merged_df = pd.merge(left=merged_df,right=results_df,on='Peptide',how='outer')
        
neuropred_results = pd.read_csv(neuropred_results_path)
neuropred_results_pt2 = pd.read_excel(neuropred_results_pt2_path, sheet_name=None)
neuropred_results_pt2_all = pd.concat([df.assign(_sheet=sheet) for sheet, df in neuropred_results_pt2.items()],ignore_index=True,sort=False)
combined = pd.concat([neuropred_results, neuropred_results_pt2_all],ignore_index=True,sort=False)
neuropred_results = pd.concat([neuropred_results,neuropred_results_pt2_all])

results_neuropred_merge = pd.merge(left=merged_df, right=neuropred_results, left_on='Peptide', right_on='Peptide sequence', how='left')

str_list = ['']
gene_dir_list = get_dir_names_with_strings_list(str_list,outside_gene_dir)

gene_matches_full = pd.DataFrame()

for x in gene_dir_list:
    gene_rep_path = f'{outside_gene_dir}\\{x}\\results\\gene_matches.csv'
    gene_rep = pd.read_csv(gene_rep_path)
    gene_matches_full = pd.concat([gene_matches_full, gene_rep],ignore_index=True)
    
for y in merged_gene_paths:
    gene_rep = pd.read_csv(y)
    gene_matches_full = pd.concat([gene_matches_full, gene_rep],ignore_index=True)
    
gene_matches_full_filtered = gene_matches_full[(gene_matches_full['matched_gene_id'] != 'NO_MATCH')]
gene_matches_full_filtered = gene_matches_full_filtered.drop(columns=['hit_start','hit_end','feature_type','gene_start','gene_end'])

results_neuropred_gene_merge = pd.merge(left=results_neuropred_merge, right=gene_matches_full_filtered, left_on='Transcript', right_on='matched_gene_id',how='left')
results_neuropred_gene_merge.to_csv(f"{results_folder}\\compiled_results_redo.csv", index=False)