# StrainSifter

A straightforward bioinformatic pipeline for detecting the presence of a bacterial strain in one or more metagenome(s).

StrainSifter is based on [Snakemake](https://snakemake.readthedocs.io/en/stable/). This pipeline allows you to output phylogenetic trees showing strain relatedness of input strains, as well as pairwise counts of single-nucleotide variants (SNVs) between input samples.

## Installation

To run StrainSifter, you must have miniconda3 and Snakemake installed.

#### Install instructions (One time only)
1. Download and install [miniconda3](https://conda.io/miniconda.html):

    For Linux:

        wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
        bash Miniconda3-latest-Linux-x86_64.sh

2. Clone the StrainSifter workflow to the directory where you wish to run the pipeline:

        git clone https://github.com/bhattlab/strainsifter

3. Create the new conda environment:

        cd strainsifter
        conda env create -f envs/environment.yaml

4. Check that Snakemake is installed and up to date:

  Snakemake should come packaged with miniconda3. To verify, type:

    snakemake --version

  If Snakemake is not installed, you will see the error:

    snakemake: command not found

  If Snakemake is not installed, install with:

    conda install snakemake

  StrainSifter requires Snakemake version 5.1.4 -- it may be necessary to update Snakemake, which you can do by typing:

    conda update snakemake

#### Activate the conda environment (Every time you use StrainSifter)

    source activate ssift

### Dependencies

We recommend running StrainSifter in the provided conda environment. If you wish to run StrainSifter without using the conda environment, the following tools must be installed and in your system PATH:
* [Burrows-Wheeler Aligner (BWA)](http://bio-bwa.sourceforge.net)
* [Samtools](http://www.htslib.org)
* [Bamtools](https://github.com/pezmaster31/bamtools)
* [Bedtools](http://bedtools.readthedocs.io/en/latest/)
* [MUSCLE](https://www.drive5.com/muscle/)
* [FastTree](http://www.microbesonline.org/fasttree/)
* [Python3](https://www.python.org/downloads/)
* [R](https://www.r-project.org)

## Usage

Due to the computing demands of the StrainSifter pipeline, we recommend running on a computing cluster if possible.
Instructions to enable Snakemake to schedule cluster jobs with SLURM can be found at https://github.com/bhattlab/slurm

### Input files

* Reference genome assembly in fasta format (can be a draft genome or a finished reference genome)
Acceptable file extensions: ".fasta", ".fa", ".fna"

* Two or more short read datasets in fastq format (metagenomic reads or isolate reads), optionally gzipped
Acceptable file extensions: ".fq", ".fastq", ".fq.gz", ".fastq.gz"

Short read data can be paired- or single-end.

You will need to indicate input files in the config file for each sample you wish to run StrainSifter on. This is described below in the *Config* section:

### Config

You must update the config.yaml file as follows:

**reference:** Path to reference genome (fasta format)

**reads:** Samples and the file path(s) to the input reads.

Optionally, you can update the following parameters:

**prefix:** (optional) desired filename for output files. If blank, the name of the reference genome will be used.

**mapq:** minimum mapping quality score to evaluate a read aligment

**n_mismatches:** consider reads with this many mismatches or fewer

**min_cvg:** minimum read depth to determine the nucleotide at any given postion

**min_genome_percent:** the minimum fraction of bases that must be covered at min_cvg or greater to process an sample

**base_freq:** minimum frequency of a nucleotide to call a base at any position

Example config.yaml:

    ##### input files #####

    # reference genome (required)
    reference: /path/to/ref.fna

    # short read data (at least two samples required)
    reads:
    sample1:
    [/path/to/sample1_R1.fq,
    /path/tp/sample1_R2.fq]
    sample2:
    [/path/to/sample2_R1.fq,
    /path/to/sample2_R2.fq]
    sample3: /path/to/sample3.fq

    # prefix for output files (optional - can leave blank)
    prefix:


    ##### workflow parameters #####

    # alignment parameters:
    mapq: 60
    n_mismatches: 5

    # variant calling parameters:
    min_cvg: 5
    min_genome_percent: 0.5
    base_freq: 0.8


### Running StrainSifter

To run StrainSifter, the config file must be present in the directory in which you wish to run the workflow.
You should then be able to run StrainSifter as follows:

#### Phylogeny

To generate a phylogenetic tree showing all of the input samples that contain your strain of interest at sufficient coverage to profile:

    snakemake {prefix}.tree.pdf

#### SNV counts

To generate a list of pairwise SNV counts between all input samples:

    snakemake {prefix}.dist.tsv

### FAQ

Q: Can StrainSifter be used for non-bacterial genomes (e.g. yeast)?

A: At present, we recommend StrainSifter for bacteria only. As yeast and other fungi can be diploid, adjustments would likely need to be made to this workflow.
