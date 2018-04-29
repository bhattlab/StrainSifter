# StrainSifter

A straightforward bioinformatic pipeline for detecting the presence of a bacterial strain in one or more metagenome(s).

Strainsifter is based on [Snakemake](https://snakemake.readthedocs.io/en/stable/). This pipeline allows you to output phylogenetic trees showing strain relatedness of input strains, as well as pairwise counts of single-nucleotide variants (SNVs) between input samples.

## Usage

### Phylogeny
snakemake tree_plot/strain.tree.pdf --configfile config.yaml

### SNV counts
snakemake snv_distances/strain.snps.tsv --configfile config.yaml

## Installation

Download and install [miniconda3](https://conda.io/miniconda.html):

For Linux:

    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh

### Dependencies

The following tools must be installed and in your system PATH to run StrainSifter:
* SAMtools
* BamTools
* BedTools
* MUSCLE
* FastTree
* R
