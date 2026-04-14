# -*- coding: utf-8 -*-
"""
Created on Fri May 16 08:55:57 2025

@author: lafields2
"""

import pandas as pd

# === CONFIG ===
blast_hits_path = r"D:\Manuscripts\2025_inSilicoPrediction\Cancer borealis prediction\blast_hits_with_scaffolds_filtered.csv" # Your formatted BLAST result file
gff3_path = r"D:\Manuscripts\2025_inSilicoPrediction\Cancer borealis prediction\JC_geneModels_braker3+manual.gff"  # Decompressed GFF
output_path = r"D:\Manuscripts\2025_inSilicoPrediction\Cancer borealis prediction\cbor_blast_gene_matches.csv"

# === Load and preprocess BLAST hits ===
blast_hits = pd.read_csv(blast_hits_path)

# Fix flipped coordinates
blast_hits['true_start'] = blast_hits[['start', 'end']].min(axis=1)
blast_hits['true_end'] = blast_hits[['start', 'end']].max(axis=1)

# === Load GFF (only gene/mRNA entries) ===
gff_rows = []
with open(gff3_path, 'r') as f:
    for line in f:
        if line.startswith('#') or line.strip() == '':
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9:
            continue
        feature_type = parts[2]
        if feature_type not in ['gene', 'mRNA']:
            continue
        scaffold = parts[0]
        start = int(parts[3])
        end = int(parts[4])
        strand = parts[6]
        attributes = parts[8]
        gene_id = attributes.split(';')[0].replace('ID=', '')
        gff_rows.append([scaffold, start, end, strand, gene_id, feature_type])

gff_df = pd.DataFrame(gff_rows, columns=['scaffold', 'start', 'end', 'strand', 'gene_id', 'feature_type'])
#%%
# === Match BLAST hits to overlapping genes ===
matches = []

for _, hit in blast_hits.iterrows():
    scaffold = hit['scaffold']
    hit_start = hit['true_start']
    hit_end = hit['true_end']
    
    overlapping = gff_df[
        (gff_df['scaffold'] == scaffold) &
        (gff_df['start'] <= hit_end) &
        (gff_df['end'] >= hit_start)
    ]
    
    if overlapping.empty:
        matches.append([
            hit['qseqid'], scaffold, hit_start, hit_end,
            'NO_MATCH', '', '', ''
        ])
    else:
        for _, gene in overlapping.iterrows():
            matches.append([
                hit['qseqid'], scaffold, hit_start, hit_end,
                gene['gene_id'], gene['feature_type'],
                gene['start'], gene['end']
            ])

# === Save results ===
matches_df = pd.DataFrame(matches, columns=[
    'query_id', 'scaffold', 'hit_start', 'hit_end',
    'matched_gene_id', 'feature_type', 'gene_start', 'gene_end'
])
#%%
matches_df.to_csv(output_path, index=False)
print(f"✔ Done! Gene matches written to {output_path}")