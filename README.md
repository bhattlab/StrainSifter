# StrainSifter

A straightforward bioinformatic pipeline for detecting the presence of a bacterial strain in one or more metagenome(s).

StrainSifter is based on [Snakemake](https://snakemake.readthedocs.io/en/stable/). This pipeline allows you to output phylogenetic trees showing strain relatedness of input strains, as well as pairwise counts of single-nucleotide variants (SNVs) between input samples.

## Installation

#### (One time only)
1. Download and install [miniconda3](https://conda.io/miniconda.html):

    For Linux:
    
        wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
        bash Miniconda3-latest-Linux-x86_64.sh

2. Clone the StrainSifter workflow to the directory where you wish to run the pipeline:

        git clone https://github.com/bhattlab/strainsifter

3. Create the new conda environment:

        cd strainsifter
        conda env create -f envs/environment.yaml

#### (Every time you use StrainSifter)
Activate the environment:

    source activate ssift
    
### Dependencies

If you wish to run StrainSifter without using the conda environgment, the following tools must be installed and in your system PATH:
* [Burrows-Wheeler Aligner (BWA)](http://bio-bwa.sourceforge.net)
* [Samtools](http://www.htslib.org)
* [Bamtools](https://github.com/pezmaster31/bamtools)
* [Bedtools](http://bedtools.readthedocs.io/en/latest/)
* [MUSCLE](https://www.drive5.com/muscle/)
* [FastTree](http://www.microbesonline.org/fasttree/)
* [Python3](https://www.python.org/downloads/)
* [R](https://www.r-project.org)

## Usage

### Input files

* Reference genome assembly in fasta format (can be a draft genome or a finished reference genome)
* Two or more short read datasets in fastq format (metagenomic reads or isolate reads)

### Config file

You must update the config.yaml file as follows:

*reference:* Path to reference genome (fasta format)
*samples:* List of input samples (fastq or fastq.gz format)
*reads_dir:* Directory containing the input samples
*paired_end:* Y for yes; N for no
*min_cvg:* Minimum coverage for calling a base

Example:

    reference: /home/tamburin/fiona/crassphage/strainsifter/ref/B_uniformis.fna
    samples: /home/tamburin/fiona/crassphage/all_samples.list
    reads_dir: /home/tamburin/fiona/crassphage/readlinks
    paired_end: Y
    min_cvg: 10

### Running StrainSifter

You should then be able to run StrainSifter as follows:

#### Phylogeny

To generate a phylogenetic tree showing all of the input samples that contain your strain of interest at sufficient coverage to profile:

    snakemake strain.tree.pdf --configfile config.yaml

#### SNV counts

To generate a list of pairwise SNV counts between all input samples:

    snakemake strain.snps.tsv --configfile config.yaml
