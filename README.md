# StrainSifter

A straightforward bioinformatic pipeline for detecting the presence of a bacterial strain in one or more metagenome(s).

Can output phylogenetic trees showing strain relatedness of input strains, as well as pairwise counts of single-nucleotide variants (SNVs) between input samples.

## Usage

### Phylogeny
snakemake tree_plot/strain.tree.pdf

### SNV counts
snakemake snv_distances/strain.snps.tsv

## Dependencies

The following tools must be installed and in your system PATH to run StrainSifter:
1. Python 3
2. SAMtools
3. BamTools
4. BedTools
5. MUSCLE
6. FastTree
