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

## Repository structure

    eukpan.sh
    scripts/
      analyze_shared_accessory.py
      extract_group_presence_absence.py
      extract_representative_sequences.py
      extract_shared_accessory.py

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

## Notes

No example genome or annotation data are included in this repository.

Users should prepare their own genome FASTA and annotation files.

If annotation files are not available, they may be generated using appropriate gene prediction tools, including ANNEVO for fungal genomes.

## License

This project is released under the MIT License.
