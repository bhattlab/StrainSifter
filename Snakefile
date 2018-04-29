import os

configfile: "config.yaml"

# read list of samples
samples = []
with open(config['samples']) as samples_list:
	for line in samples_list:
		sample = line.rstrip('\n')
		samples += [sample]

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
		reads = "input_samples/{sample}.fq",
		# r2 = "input_samples/{sample}_PE2.fq"
	output:
		"filtered_bam/{sample}.filtered.bam"
	resources:
		mem=32,
		time=6
	threads: 8
	shell:
		"bwa mem -t {threads} {input.ref} {input.reads} | samtools view -b -q 60 | bamtools filter -tag 'NM:<2' | samtools sort --threads {threads} -o {output}"

rule calc_coverage:
	input:
		"filtered_bam/{sample}.filtered.bam"
		# expand("filtered_bam/{sample}.{organism}.filtered.bam", sample = samples, organism = organisms)
	output:
		"coverage/{sample}.cvg"
	resources:
		mem=16,
		time=1,
	threads: 1
	params:
		cvg=5
	shell:
		"bedtools genomecov -ibam {input} | python3 scripts/getCoverage.py {params.cvg} {input} > {output}"

rule filter_samples:
	input:
		rules.calc_coverage.output
		# expand("calc_coverage/{sample}.{organism}.list", sample = samples, organism = organisms)
	output:
		"passed_samples/{sample}.bam"
		# dynamic("passed_samples/{sample}.{organism}.passed")
	resources:
		mem=1,
		time=1
	threads: 1
	params:
		min_cvg=5,
		min_perc=0.5
	shell:
		"if ($(cat {input} | awk '{{if ($2 >= {params.min_cvg} && $3 >= {params.min_perc}) print \"true\"; else print \"false\"}}') == true);"\
		"then ln -s $PWD/filtered_bam/{wildcards.sample}.filtered.bam {output}"\
		" && touch -h {output};"\
		"else touch {output}; fi"

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
		bam=rules.filter_samples.output,
		ref=config['reference'],
		index=rules.faidx.output
	output: "pileup/{sample}.pileup"
	resources:
		mem=32,
		time=1
	threads: 1
	shell:
		"samtools mpileup -f {input.ref} -B -aa -o {output} {input.bam}"

rule call_snps:
	input: rules.pileup.output
	output: "snp_calls/{sample}.tsv"
	resources:
		mem=32,
		time=2
	threads: 1
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

rule core_snps:
	input:
		lambda wildcards: expand("consensus/{sample}.txt", sample=samples)
	output: "core_snps/all.tsv"
	resources:
		mem=16,
		time=1
	threads: 1
	shell:
		"paste $(find consensus/*.txt -size +1) | python3 snv_wf/scripts/findCoreSNPs.py > {output}"

rule core_snps_to_fasta:
	input: rules.core_snps.output
	output: "fasta/all.fasta"
	resources:
		mem=16,
		time=1
	threads: 1
	script:
		"scripts/coreSNPs2fasta.py"

rule multi_align:
	input: rules.core_snps_to_fasta.output
	output: "afa/all.afa"
	resources:
		mem=200,
		time=12
	threads: 1
	shell:
		"muscle -in {input} -out {output}"

rule build_tree:
	input: rules.multi_align.output
	output: "tree/all.tree"
	resources:
		mem=16,
		time=1
	threads: 1
	shell:
		"fasttree -nt {input} > {output}"

rule plot_tree:
	input:
		rules.build_tree.output,
		"/home/tamburin/fiona/bacteremia/bacteremia_metadata_all.csv",
	output: "tree_plot/all.tree.pdf"
	resources:
		mem=8,
		time=1
	threads: 1
	script:
		"scripts/renderTree.R"

rule pairwise_snvs:
	input: lambda wildcards: expand("consensus/{sample}.filtered.tmp", sample=samples)
	output: "snv_distances/all.snps.tsv"
	resources:
		mem=8,
		time=1
	threads: 1
	script:
		"scripts/pairwiseDist.py"
