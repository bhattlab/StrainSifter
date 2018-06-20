import os
import re
from snakemake.utils import min_version

##### set minimum snakemake version #####
min_version("5.1.4")

##### load config file and sample list #####
configfile: "config.yaml"

##### read list of input samples #####
samples = []
with open(config['samples']) as samples_list:
	for line in samples_list:
		sample = line.rstrip('\n')
		samples += [sample]

##### list of samples that meet coverage criteria -- initially empty #####
passed_samples = []
try:
	f = open("passed_samples.list")
	for line in f:
		sample = line.rstrip('\n')
		passed_samples += [sample]
except (IOError, OSError) as e:
	print("Samples meeting {cvg}X coverage not yet determined".format(cvg=config['min_cvg']))
else:
	f.close()

##### prefix for phylogenetic tree and SNV distance files #####
if config['prefix'] is None:
	prefix = re.split("/|\.", config['reference'])[-2]
else:
	prefix = config['prefix']

##### handle paired-end vs. single-end files #####
if 'reads2' in config and not config['reads2'] is None:
	reads = [config['reads1'], config['reads2']]
else:
	reads = [config['reads1']]

##### rules #####

rule bwa_index:
	input: config['reference']
	output:
		"{ref}.amb".format(ref=config['reference']),
		"{ref}.ann".format(ref=config['reference']),
		"{ref}.bwt".format(ref=config['reference']),
		"{ref}.pac".format(ref=config['reference']),
		"{ref}.sa".format(ref=config['reference'])
	shell:
		"bwa index {input}"

rule bwa_align:
	input:
		ref = config['reference'],
		ref_index = rules.bwa_index.output,
		r = lambda wildcards: expand(reads, sample=wildcards.sample)
	output:
		"filtered_bam/{sample}.filtered.bam"
	resources:
		mem=32,
		time=6
	threads: 8
	params:
		qual=config['mapq'],
		nm=config['n_mismatches']
	shell:
		"bwa mem -t {threads} {input.ref} {input.r} | "\
		"samtools view -b -q {params.qual} | "\
		"bamtools filter -tag 'NM:<={params.nm}' | "\
		"samtools sort --threads {threads} -o {output}"

rule genomecov:
	input:
		rules.bwa_align.output
	output:
		"genomecov/{sample}.tsv"
	resources:
		mem=16,
		time=1,
	threads: 1
	shell:
		"bedtools genomecov -ibam {input} > {output}"

rule calc_coverage:
	input:
		rules.genomecov.output
	output:
		"coverage/{sample}.cvg"
	resources:
		mem=16,
		time=1,
	threads: 1
	params:
		cvg=config['min_cvg']
	script:
		"scripts/getCoverage.py"

rule filter_samples:
	input: expand("coverage/{sample}.cvg", sample = samples)
	output:
		dynamic("passed_samples/{sample}.bam")
		# "passed_samples/{sample}.bam"
	resources:
		mem=1,
		time=1
	threads: 1
	params:
		min_cvg=config['min_cvg'],
		min_perc=config['min_genome_percent']
	run:
		samps = input
		# out_file = open("passed_samples.list", "a")
		for samp in samps:
			with open(samp) as s:
				cvg, perc = s.readline().rstrip('\n').split('\t')
			if (float(cvg) >= params.min_cvg and float(perc) > params.min_perc):
				shell("ln -s $PWD/filtered_bam/{s}.filtered.bam passed_samples/{s}.bam".format(s=os.path.basename(samp).rstrip(".cvg")))
				# shell("ln -s $PWD/filtered_bam/{s}.filtered.bam passed_samples/{s}.bam".format(s=os.path.basename(samp).rstrip(".cvg")))
				# out_file.write(os.path.basename(samp).rstrip(".cvg") + '\n')
				# shell("echo {p} >> output".format(p=samp.rstrip(".cvg")))

rule faidx:
	input: config['reference']
	output: "{ref}.fai".format(ref=config['reference'])
	resources:
		mem=8,
		time=1
	shell:
		"samtools faidx {input}"

rule pileup:
	input:
		# bam="filtered_bam/{sample}.filtered.bam",
		bam="passed_samples/{sample}.bam",
		ref=config['reference'],
		index=rules.faidx.output
	output: "pileup/{sample}.pileup"
	resources:
		mem=32,
		time=1
	threads: 16
	shell:
		"samtools mpileup -f {input.ref} -B -aa -o {output} {input.bam}"

rule call_snps:
	input: rules.pileup.output
	output: "snp_calls/{sample}.tsv"
	resources:
		mem=32,
		time=2
	threads: 16
	params:
		min_cvg=5,
		min_freq=0.8,
		min_qual=20
	script:
		"scripts/callSNPs.py"

rule snp_consensus:
	input: rules.call_snps.output
	output: "consensus/{sample}.txt"
	resources:
		mem=32,
		time=2
	threads: 1
	params:
		min_cvg=5,
		min_freq=0.8,
		min_phred=20
	shell:
		"(echo {wildcards.sample}; cut -f4 {input}) > {output}"

# rule all:
#     input: dynamic("consensus/{sample}.txt")
# 	output: "done"
# 	shell:
# 		"touch output"

# rule all_bam:
#     input: expand("coverage/{sample}.cvg", sample = samples)
# 	output: "done2"
# 	shell:
# 		"touch output"

rule combine:
	input:
		dynamic("consensus/{sample}.txt")
		# dynamic(expand("consensus/{sample}.txt", sample = passed_samples))
		# samps=dynamic(expand(rules.call_snps.output, sample = passed_samples)),
		# filt=dynamic(rules.filter_samples.output)
	output: "{name}.cns.tsv".format(name = prefix)
	resources:
		mem=2,
		time=1
	threads: 1
	shell:
		"paste consensus/* > {output}"

rule core_snps:
	input: rules.combine.output
	output: "{name}.core_snps.tsv".format(name = prefix)
	resources:
		mem=16,
		time=1
	threads: 1
	script:
		"scripts/findCoreSNPs.py"

rule core_snps_to_fasta:
	input: rules.core_snps.output
	output: "{name}.fasta".format(name = prefix)
	resources:
		mem=16,
		time=1
	threads: 1
	script:
		"scripts/coreSNPs2fasta.py"

rule multi_align:
	input: rules.core_snps_to_fasta.output
	output: "{name}.afa".format(name = prefix)
	resources:
		mem=200,
		time=12
	threads: 1
	shell:
		"muscle -in {input} -out {output}"

rule build_tree:
	input: rules.multi_align.output
	output: "{name}.tree".format(name = prefix)
	resources:
		mem=16,
		time=1
	threads: 1
	shell:
		"fasttree -nt {input} > {output}"

rule plot_tree:
	input: rules.build_tree.output
	output: "{name}.tree.pdf".format(name = prefix)
	resources:
		mem=8,
		time=1
	threads: 1
	script:
		"scripts/renderTree.R"

rule pairwise_snvs:
	input: expand("consensus/{sample}.txt", sample = passed_samples)
	output: "{name}.dist.tsv".format(name = prefix)
	resources:
		mem=8,
		time=1
	threads: 1
	script:
		"scripts/pairwiseDist.py"
