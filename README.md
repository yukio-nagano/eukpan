# EukPan

EukPan is a command-line pipeline for eukaryotic pangenome analysis.

It uses genome FASTA files and annotation files in GFF, GFF3, or GTF format to generate protein datasets, infer orthogroups, construct a concatenated single-copy core protein alignment, and analyze shared accessory genome presence/absence patterns.

## Main features

EukPan performs the following steps:

1. Converts GFF/GFF3/GTF files using AGAT
2. Keeps the longest isoform for each gene
3. Extracts protein sequences using gffread
4. Simplifies protein FASTA headers
5. Runs OrthoFinder in a separate conda environment
6. Extracts single-copy core orthologues
7. Aligns single-copy orthologues using MAFFT
8. Trims alignments using trimAl
9. Concatenates single-copy core protein alignments
10. Extracts shared accessory orthogroups
11. Visualizes accessory genome presence/absence patterns
12. Provides optional downstream scripts for group-specific orthogroup extraction, representative sequence extraction, and PCA of accessory genome presence/absence profiles

## Repository structure

    eukpan.sh
    scripts/
      analyze_shared_accessory.py
      extract_group_presence_absence.py
      extract_representative_sequences.py
      extract_shared_accessory.py
      pca_accessory_presence_absence.py

## Installation

Download EukPan from GitHub:

    git clone https://github.com/yukio-nagano/eukpan.git
    cd eukpan

## Requirements

Create a dedicated OrthoFinder environment:

    conda create -n orthofinder \
      -c conda-forge -c bioconda \
      orthofinder

Create the main EukPan environment:

    conda create -n eukpan \
      -c conda-forge -c bioconda \
      python \
      agat \
      gffread \
      seqkit \
      mafft \
      parallel \
      trimal \
      iqtree \
      modeltest-ng \
      pandas \
      numpy \
      scipy \
      matplotlib

Optional: install Inkscape for EMF export.

    sudo apt update
    sudo apt install inkscape

Optional: install R and ggplot2 for PCA visualization.

    sudo apt update
    sudo apt install r-base

Then, in R:

    install.packages("ggplot2")

## Input

Prepare two input directories.

One directory should contain annotation files:

    gff/
      sample1.gff3
      sample2.gff3
      sample3.gff3

The annotation files can be obtained from public genome databases or generated from genome FASTA files using gene prediction tools.

For fungal genomes, GFF files can also be generated using ANNEVO. For example:

    python ~/ANNEVO/annotation.py \
      --genome assembled_genomes/sample1.fna \
      --model_path ~/ANNEVO/ANNEVO_model/ANNEVO_Fungi.pt \
      --output gff/sample1.gff \
      --threads 32

The other directory should contain genome FASTA files:

    assembled_genomes/
      sample1.fna
      sample2.fna
      sample3.fna

The file basenames must match.

For example:

    gff/sample1.gff3
    assembled_genomes/sample1.fna

or, when ANNEVO is used:

    gff/sample1.gff
    assembled_genomes/sample1.fna

Supported annotation extensions:

    .gff
    .gff3
    .gtf

Supported genome FASTA extensions:

    .fna
    .fa
    .fasta

## Usage

Activate the main environment:

    conda activate eukpan

Run EukPan:

    bash eukpan.sh \
      --gff-dir gff \
      --genome-dir assembled_genomes \
      --outdir pangenome_results \
      --threads 16 \
      --orthofinder-env orthofinder \
      --scripts-dir scripts

## Main outputs

    pangenome_results/proteomes/
    pangenome_results/proteomes/OrthoFinder/
    pangenome_results/phylo/concat_all.fa
    pangenome_results/phylo/concat_partitions.txt
    pangenome_results/accessory_results/shared_accessory_presence_absence.tsv
    pangenome_results/accessory_analysis/

## Accessory genome visualizations

The following files are generated in:

    pangenome_results/accessory_analysis/

Main figures:

    01_presence_distribution.png
    02_presence_absence_heatmap_reordered.png
    03_jaccard_distance_matrix_reordered.png
    04_hierarchical_clustering.png
    05_dendrogram_heatmap_combined.png

The reordered Jaccard distance matrix is drawn so that the matrix panel itself is square, with the color scale placed outside the matrix.

## Group-specific orthogroup analysis

Prepare GROUP_SAMPLES.txt with one sample name per line.

Example:

    sample1
    sample2
    sample3

Then run:

    python scripts/extract_group_presence_absence.py \
      -i pangenome_results/proteomes/OrthoFinder/Results_*/Orthogroups/Orthogroups.GeneCount.tsv \
      -g GROUP_SAMPLES.txt \
      -o group_biased_presence_absence.tsv \
      --group-min-pct 90 \
      --group-max-pct 100 \
      --nongroup-min-pct 0 \
      --nongroup-max-pct 10

This example extracts orthogroups present in 90-100% of the focal group and 0-10% of the non-group samples.

## Representative sequence extraction

After group-specific orthogroup extraction, representative sequences can be extracted using:

    python scripts/extract_representative_sequences.py \
      -i group_biased_presence_absence.tsv \
      -s pangenome_results/proteomes/OrthoFinder/Results_*/Orthogroup_Sequences \
      -o group_biased_representatives.fa

## PCA of accessory orthogroup presence/absence profiles

EukPan also provides an optional script for principal component analysis of shared accessory orthogroup presence/absence profiles.

This script is independent of the main pipeline and can be run after EukPan has generated:

    pangenome_results/accessory_results/shared_accessory_presence_absence.tsv

Run PCA without group labels:

    python scripts/pca_accessory_presence_absence.py \
      -i pangenome_results/accessory_results/shared_accessory_presence_absence.tsv \
      -o pangenome_results/accessory_pca

Run PCA with group labels:

    python scripts/pca_accessory_presence_absence.py \
      -i pangenome_results/accessory_results/shared_accessory_presence_absence.tsv \
      -o pangenome_results/accessory_pca \
      -g GROUP_SAMPLES.txt \
      --group-name Group_A \
      --nongroup-name Non_group_A

The PCA script generates R-ready output files:

    pangenome_results/accessory_pca/accessory_pca_scores.tsv
    pangenome_results/accessory_pca/accessory_pca_loadings.tsv
    pangenome_results/accessory_pca/accessory_pca_variance.tsv
    pangenome_results/accessory_pca/accessory_pca_metadata.tsv
    pangenome_results/accessory_pca/accessory_presence_absence_matrix_for_R.tsv
    pangenome_results/accessory_pca/accessory_pca_used_orthogroups.tsv
    pangenome_results/accessory_pca/plot_accessory_pca.R

To visualize the PCA in R:

    cd pangenome_results/accessory_pca
    Rscript plot_accessory_pca.R

This generates:

    accessory_pca_PC1_PC2.pdf
    accessory_pca_PC1_PC2.png

## Notes

No example genome or annotation data are included in this repository.

Users should prepare their own genome FASTA and annotation files.

If annotation files are not available, they may be generated using appropriate gene prediction tools, including ANNEVO for fungal genomes.

The PCA script is intended as an optional downstream analysis. It does not change the main EukPan pipeline output.

## License

This project is released under the MIT License.
