# Genome-Guided Neuropeptide Prediction in *Cancer borealis* and *Callinectes sapidus*

This repository contains the Python code used for genome-guided neuropeptide prediction, result compilation, and figure generation for our ACS Chemical Biology manuscript on in silico neuropeptide discovery and LC-MS/MS validation in two decapod crustaceans: *Cancer borealis* (Jonah crab) and *Callinectes sapidus* (blue crab).

The code in this repository supports three major parts of the study:

1. Genome-guided candidate discovery from public genome assemblies and annotation resources.
2. Post-processing and validation of predicted peptides against LC-MS/MS datasets.
3. Preparation of manuscript figures and summary tables.

## Repository Contents

```text
GitHubRepo/
|-- pipeline/
|   |-- complete_workflow_Cborealis.py
|   |-- complete_workflow_Csapidus.py
|   `-- modules/
|       |-- scaffold_from_alignment.py
|       |-- scaffold_gff_extract.py
|       |-- BLAST_to_Gene_Cborealis.py
|       |-- BLAST_to_Gene_Csapidus.py
|       |-- protein_from_scaffold_Cborealis.py
|       `-- protein_from_scaffold_Csapidus.py
|-- analysis/
|   |-- Cborealis_results_compilation.py
|   |-- Csapidus_results_compilation.py
|   |-- comparison_script.py
|   |-- full_FDR_to_EG.py
|   |-- fuzzy_align.py
|   |-- motif_fuzzy_align.py
|   `-- verification.py
|-- figures/
|   |-- alignment_histogram.py
|   |-- signal_peptide_distribution.py
|   `-- venn_diagram_generate.py
`-- utils/
    `-- find_crz_candidates.py
```

## What Each Directory Does

### `pipeline/`
Primary scripts for genome-guided neuropeptide discovery.

- `complete_workflow_Cborealis.py`: integrates alignment outputs, genome scaffold mapping, gene extraction, and protein-level candidate generation for *C. borealis*.
- `complete_workflow_Csapidus.py`: corresponding workflow for *C. sapidus*.
- `modules/scaffold_from_alignment.py`: maps BLAST hit accessions to scaffold identifiers using NCBI sequence report metadata.
- `modules/scaffold_gff_extract.py`: extracts scaffold-specific GFF entries for downstream sequence parsing.
- `modules/BLAST_to_Gene_Cborealis.py` and `modules/BLAST_to_Gene_Csapidus.py`: match BLAST hit coordinates to overlapping gene models in GFF annotations.
- `modules/protein_from_scaffold_Cborealis.py` and `modules/protein_from_scaffold_Csapidus.py`: reconstruct CDS regions and translate predicted proteins from genome and annotation files.

### `analysis/`
Scripts used to compile, compare, and validate peptide prediction results.

- `Cborealis_results_compilation.py`: compiles EndoGenius search results and merges them with *C. borealis* NeuroPred and gene-match outputs.
- `Csapidus_results_compilation.py`: same general workflow for *C. sapidus*.
- `full_FDR_to_EG.py`: filters peptide-spectrum matches by score and reformats outputs for EndoGenius-based downstream analysis.
- `fuzzy_align.py`: performs sequence similarity comparisons between predicted peptides and known neuropeptide databases.
- `motif_fuzzy_align.py`: motif-level fuzzy matching against a curated motif database.
- `comparison_script.py`: compares predicted and verified peptide lists, including exact and substring matches with optional I/L normalization.
- `verification.py`: quality-control and verification utilities for sequence and annotation outputs.

### `figures/`
Scripts used to generate figures for manuscript preparation.

- `alignment_histogram.py`: plots similarity score distributions.
- `signal_peptide_distribution.py`: summarizes signal peptide length distributions.
- `venn_diagram_generate.py`: generates Venn diagrams comparing peptide sets.

### `utils/`
Small supporting utilities.

- `find_crz_candidates.py`: searches translated reading frames for corazonin-like sequence motifs.

## Workflow Summary

The overall workflow used in the manuscript is:

1. Assemble a neuropeptide precursor query database from known crustacean and related precursor sequences.
2. Search genome assemblies and/or annotated gene models for candidate neuropeptide loci.
3. Map significant hits to scaffolds and genes using NCBI metadata and genome annotation files.
4. Reconstruct candidate protein sequences from matched loci.
5. Identify signal peptide-containing candidates and annotate likely mature peptides.
6. Compare predicted peptides to known neuropeptide databases and motif collections.
7. Validate predicted peptides using LC-MS/MS data processed with EndoGenius.
8. Compile final candidate tables and generate manuscript figures.

Some steps in the full study relied on manual review, external software, and iterative filtering in addition to the Python scripts deposited here.

## External Software and Previously Published Code

This repository contains the custom scripts written for the manuscript, but the complete analysis also depended on external software and previously published tools.

### Required external tools

- NCBI BLAST+ for sequence similarity searches.
- SignalP for signal peptide prediction.
- NeuroPred for neuropeptide cleavage prediction.
- Standard Python scientific packages such as `pandas`, `biopython`, `matplotlib`, `numpy`, `rapidfuzz`, and `regex`.

### Previously published tools used in the analysis

- EndoGenius: [lingjunli-research/EndoGenius](https://github.com/lingjunli-research/EndoGenius)
- MotifQuest: [lingjunli-research/MotifQuest](https://github.com/lingjunli-research/MotifQuest)

EndoGenius and MotifQuest are not introduced here as new software in this repository; they are separate previously published tools used during the MS-based validation and motif-oriented analysis stages. Please cite those repositories and their associated publications as appropriate when reusing that portion of the workflow.

## Important Reproducibility Notes

Most scripts in this repository were developed for direct use on the project data directories used during the study. As a result:

- Many scripts currently contain hard-coded absolute file paths.
- Several scripts assume specific intermediate file names produced in earlier pipeline steps.
- Some downstream steps depend on files generated outside this repository, including genome resources, EndoGenius outputs, and manuscript working files.
- A subset of the analysis involved manual curation of candidate lists and interpretation of outputs.

Accordingly, this repository is best understood as the deposited analysis code used for the manuscript rather than a fully packaged one-command pipeline.

## Input Data Used by the Workflow

The workflow draws from several classes of input data:

- Public genome assemblies and gene annotations for *C. borealis* and *C. sapidus*.
- Query neuropeptide precursor FASTA databases.
- Signal peptide and cleavage prediction outputs.
- LC-MS/MS search outputs used for validation and feeding-study analyses.
- Curated peptide and motif reference databases used for similarity comparisons.

Large raw and intermediate data files are not included directly in this code repository.

## How To Reuse This Repository

Users interested in reproducing or adapting the workflow should generally:

1. Review the scripts in `pipeline/` to identify the organism-specific entry point.
2. Replace hard-coded paths with local paths to genome FASTA, GFF, BLAST outputs, and result folders.
3. Run the relevant module scripts in the order implied by the workflow.
4. Use the `analysis/` scripts to merge prediction outputs with EndoGenius and comparison results.
5. Use the `figures/` scripts to regenerate summary plots after updating input paths.

Because the workflow was executed in a manuscript project environment, some adaptation will be required before running it on a new system.

## Suggested Python Environment

A typical environment for these scripts includes:

- Python 3.10+
- pandas
- biopython
- matplotlib
- numpy
- rapidfuzz
- regex
- openpyxl

Additional system-level dependencies such as BLAST+ and access to SignalP/NeuroPred are required for the full end-to-end workflow.

## Citation

If you use this repository, please cite:

- the associated ACS Chemical Biology manuscript describing this workflow and dataset
- EndoGenius
- MotifQuest
- any public genome and annotation resources used in your analysis

A manuscript DOI and full citation can be added here once available.

## Contact

For questions about this repository or the manuscript analysis workflow, please contact the corresponding authors of the study.
