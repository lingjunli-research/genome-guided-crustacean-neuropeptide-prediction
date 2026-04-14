# -*- coding: utf-8 -*-
"""
Created on Sun Oct 12 09:29:33 2025

@author: lafields2
"""

import pandas as pd
import os
# sample_name_list =['Brain_Fed_TR2',
#                    'Brain_Fed_TR3',
#                    'Brain_Fed_TR4',
#                    'Brain_Unfed_TR1',
#                    'Brain_Unfed_TR2',
#                    'Brain_Unfed_TR3',
#                    'COG_Fed_TR1',
#                    'COG_Fed_TR2',
#                    'COG_Fed_TR3',
#                    'COG_Unfed_TR1',
#                    'COG_Unfed_TR2',
#                    'COG_Unfed_TR3',
#                    'PO_Fed_TR1',
#                    'PO_Fed_TR2',
#                    'PO_Fed_TR3',
#                    'PO_Unfed_TR1',
#                    'PO_Unfed_TR2',
#                    'PO_Unfed_TR3',
#                    'SG_Fed_TR1',
#                    'SG_Fed_TR2',
#                    'SG_Fed_TR3',
#                    'SG_Unfed_TR1',
#                    'SG_Unfed_TR2',
#                    'SG_Unfed_TR3',
#                    'TG_Fed_TR1',
#                    'TG_Fed_TR2',
#                    'TG_Fed_TR3',
#                    'TG_Unfed_TR1',
#                    'TG_Unfed_TR2',
#                    'TG_Unfed_TR3']

parent_dir = r"D:\Manuscripts\2024_Multiplexed_Feeding\LumosRun_20250830\EG_search_out\QuantEval\updated_EG"

def get_dir_names_with_strings_list(str_list,results_directory):
    full_list = [name for name in os.listdir(results_directory)
                 if os.path.isdir(os.path.join(results_directory, name))]
    return [nm for ps in str_list for nm in full_list if ps in nm]


str_list = ['']
sample_name_list = get_dir_names_with_strings_list(str_list,parent_dir)

def fdr_to_EG(parent_dir, sample):
    results_dir = f"{parent_dir}\\{sample}\\{sample}"
    score_threshold = 750
    
    scores_path = f"{results_dir}\\IDed_peptide_scores.csv"
    scores = pd.read_csv(scores_path)
    scores_filtered = scores[scores['Final score, run:1'] >= score_threshold]
    
    psm_rep_path = f"{results_dir}\\final_psm_report_out.csv"
    psm_rep = pd.read_csv(psm_rep_path)
    
    merged_rep = pd.merge(left=scores_filtered,right=psm_rep, 
                          left_on=['Peptide','Scan'],
                          right_on=['Sequence with mod','Scan'])

    merged_rep.rename(columns={'Sequence': 'Unmodified sequence'}, inplace=True)
    
    merge_rep_format = merged_rep.drop(
        columns=[
            'Sequence coverage',
            'Correlation value',
            'Sequence with mod',
            'count',
            'Motif',
            'Motif status',
            'Motif Length',
            'Motif:Seq Ratio',
            '# Matching motifs',
            'Sequence/Motif Coverage',
            'Seq Length',
            'Step assigned',
        ]
    )
    
    merge_rep_target= merge_rep_format[merge_rep_format['Status'] == 'Target']
    
    spectra_path = f'{results_dir}_formatted.txt'
    spectra = pd.read_csv(spectra_path, sep=',')
    spectra_drop_dups = spectra.drop_duplicates(subset='ms2_scan')
    spectra_filtered = spectra_drop_dups.drop(
        columns=[
            'fragment_mz',
            'fragment_intensity',
            'fragment_z',
            'fragment_resolution',
            'precursor_mz',
            'precursor_z'
        ]
    )
    
    merge_rep_w_spectra = pd.merge(left=merge_rep_target,right=spectra_filtered,left_on='Scan',right_on='ms2_scan')
    merge_rep_w_spectra_filtered = merge_rep_w_spectra.drop(
        columns=[
            'ms2_scan',
        ]
    )
    output_path = f'{results_dir}\\final_results_EG_score.csv'
    merge_rep_w_spectra_filtered.to_csv(output_path)
    
for x in sample_name_list:
    fdr_to_EG(parent_dir, x)