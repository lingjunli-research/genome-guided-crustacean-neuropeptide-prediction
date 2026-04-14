# -*- coding: utf-8 -*-
"""
Created on Wed Jun 11 09:46:21 2025

@author: lafields2
"""

import os

# File paths
scaffold_list_path = r"D:\Manuscripts\2025_inSilicoPrediction\Cancer borealis prediction\matched_genes.txt"
input_gff = r"D:\Manuscripts\2025_inSilicoPrediction\Cancer borealis prediction\JC_geneModels_braker3+manual.gff"
output_dir = r"D:\Manuscripts\2025_inSilicoPrediction\Cancer borealis prediction\scaffolds"

# Load scaffold names from text file
with open(scaffold_list_path, 'r') as f:
    scaffolds_of_interest = set(line.strip() for line in f if line.strip())

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Create a dictionary to hold open file handles
output_files = {
    scaffold: open(os.path.join(output_dir, f"{scaffold}.gff"), 'w')
    for scaffold in scaffolds_of_interest
}

# Process input GFF
with open(input_gff, 'r') as fin:
    for line in fin:
        if line.startswith('#') or line.strip() == '':
            for fout in output_files.values():
                fout.write(line)
        else:
            scaffold = line.split('\t')[0]
            if scaffold in scaffolds_of_interest:
                output_files[scaffold].write(line)

# Close all file handles
for fout in output_files.values():
    fout.close()

print(f"✔ Wrote {len(scaffolds_of_interest)} GFF files to {output_dir}")
