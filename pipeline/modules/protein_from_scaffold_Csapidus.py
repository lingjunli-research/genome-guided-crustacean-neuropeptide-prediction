# -*- coding: utf-8 -*-
"""
Created on Wed Jun 11 09:52:02 2025

@author: lafields2
"""

from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import pandas as pd

blast_gene_matches_path = r"D:\Manuscripts\2025_inSilicoPrediction\Cancer borealis prediction\cbor_blast_gene_matches_filtered.csv"
genome_fasta = r"D:\Manuscripts\2025_inSilicoPrediction\Cancer borealis prediction\JC-genome_v3-230128.fasta\JC-genome_v3-230128.fasta"
gff_directory = r"D:\Manuscripts\2025_inSilicoPrediction\Cancer borealis prediction\scaffolds"
output_fasta_directory = r"D:\Manuscripts\2025_inSilicoPrediction\Cancer borealis prediction\scaffolds"

genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))

blast_matches = pd.read_csv(blast_gene_matches_path)
blast_matches = blast_matches[blast_matches['feature_type'] == 'gene']
blast_matches = blast_matches[blast_matches['matched_gene_id'] != 'NO_MATCH']

scaffolds = blast_matches['scaffold'].values.tolist()
scaffolds_nodups = set(scaffolds)

for x in scaffolds_nodups:
    blast_filtered = blast_matches[blast_matches['scaffold'] == x]
    genes = blast_filtered['matched_gene_id'].values.tolist()
    gene_ids_to_extract = set(genes)

    gff_file = f'{gff_directory}\\{x}.GFF'


    # === 3. Parse GFF and collect CDS features by gene ===
    cds_by_gene = defaultdict(list)
    
    with open(gff_file) as gff:
        for line in gff:
            if line.startswith("#"):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
    
            chrom, _, feature_type, start, end, _, strand, _, attributes = parts
    
            if feature_type != "CDS":
                continue
    
            gene_id = None
            for attr in attributes.split(';'):
                if attr.startswith("Parent="):
                    # remove transcript suffix like -RA if present
                    gene_id = attr.replace("Parent=", "").split("-")[0]
                    break
    
            if gene_id and gene_id in gene_ids_to_extract:
                cds_by_gene[gene_id].append({
                    'chrom': chrom,
                    'start': int(start),
                    'end': int(end),
                    'strand': strand
                })
    
    print(f"✔ Found CDS entries for {len(cds_by_gene)} matching genes.")
    output_fasta = f'{output_fasta_directory}\\{x}.fa'
    # === 4. Extract and translate protein sequences ===
    with open(output_fasta, "w") as out_fasta:
        for gene_id, regions in cds_by_gene.items():
            regions.sort(key=lambda r: r['start'])
    
            seq_parts = []
            for r in regions:
                subseq = genome[r['chrom']].seq[r['start'] - 1 : r['end']]
                seq_parts.append(subseq)
    
            full_seq = Seq("").join(seq_parts)
    
            if regions[0]['strand'] == "-":
                full_seq = full_seq.reverse_complement()
    
            try:
                protein_seq = full_seq.translate(to_stop=True)
            except Exception as e:
                print(f"⚠ Translation failed for {gene_id}: {e}")
                continue
    
            out_fasta.write(f">{gene_id}\n{protein_seq}\n")

print(f"✅ Protein FASTA written to: {output_fasta}")
