import os

# extra credit:
## bwa index rule

# To do:
## rule to rename trees based on labels
## fix calc_coverage filename input

# read list of stool samples
organism_to_samples = {}
samples_to_organism = {}
with open('samples_with_enough_cvg.list') as samples_list:
	for line in samples_list:
		sample, organism = line.rstrip('\n').split('.')

		if not sample.startswith('S'):
			samples_to_organism[sample] = organism

		if organism in organism_to_samples:
			organism_to_samples[organism] += [sample]
		else:
			organism_to_samples[organism] = [sample]

rule bwa_index:
	input: "/home/tamburin/fiona/bacteremia/14.assemble_trimmed/filtered_assemblies/{assembly}.fasta"
	output:
		"/home/tamburin/fiona/bacteremia/14.assemble_trimmed/filtered_assemblies/{assembly}.fasta.amb",
		"/home/tamburin/fiona/bacteremia/14.assemble_trimmed/filtered_assemblies/{assembly}.fasta.ann",
		"/home/tamburin/fiona/bacteremia/14.assemble_trimmed/filtered_assemblies/{assembly}.fasta.bwt",
		"/home/tamburin/fiona/bacteremia/14.assemble_trimmed/filtered_assemblies/{assembly}.fasta.pac",
		"/home/tamburin/fiona/bacteremia/14.assemble_trimmed/filtered_assemblies/{assembly}.fasta.sa"
	shell:
		"bwa index {input}"

rule bwa_align_stool:
	input:
		ref = "/home/tamburin/fiona/bacteremia/14.assemble_trimmed/filtered_assemblies/{assembly}.fasta",
		ref_index = rules.bwa_index.output,
		r1 = "/home/tamburin/fiona/bacteremia/stool_reads_all/{sample}_PE1.fq",
		r2 = "/home/tamburin/fiona/bacteremia/stool_reads_all/{sample}_PE2.fq"
	output:
		"filtered_bam/{sample}.{assembly}.filtered.bam"
	resources:
		mem=32,
		time=6
	threads: 8
	shell:
		"bwa mem -t {threads} {input.ref} {input.r1} {input.r2} | samtools view -b -q 60 | bamtools filter -tag 'NM:<2' | samtools sort --threads {threads} -o {output}"

rule bwa_align_isolate:
	input:
		ref = "/home/tamburin/fiona/bacteremia/14.assemble_trimmed/filtered_assemblies/{assembly}.fasta",
		ref_index = rules.bwa_index.output,
		r1 = "/home/tamburin/fiona/bacteremia/isolate_reads_all/{sample}_PE1.fq",
		r2 = "/home/tamburin/fiona/bacteremia/isolate_reads_all/{sample}_PE2.fq"
	output:
		"filtered_bam/{sample}.{assembly}.filtered.bam"
	resources:
		mem=32,
		time=6
	threads: 8
	shell:
		"bwa mem -t {threads} {input.ref} {input.r1} {input.r2} | samtools view -b -q 60 | bamtools filter -tag 'NM:<2' | samtools sort --threads {threads} -o {output}"

# rule calc_coverage:
# 	input:
# 		"filtered_bam/{sample}.{organism}.filtered.bam"
# 		# expand("filtered_bam/{sample}.{organism}.filtered.bam", sample = samples, organism = organisms)
# 	output:
# 		"calc_coverage/{sample}.{organism}.list"
# 	resources:
# 		mem=16,
# 		time=1,
# 	threads: 1
# 	shell:
# 		"bedtools genomecov -ibam {input} | python3 getCoverage.py 5 {input} > {output}"

# rule filter_samples:
# 	input:
# 		rules.calc_coverage.output
# 		# expand("calc_coverage/{sample}.{organism}.list", sample = samples, organism = organisms)
# 	output:
# 		"passed_samples/{sample}.{organism}.bam"
# 		# dynamic("passed_samples/{sample}.{organism}.passed")
# 	resources:
# 		mem=1,
# 		time=1
# 	threads: 1
# 	params:
# 		min_cvg=5,
# 		min_perc=0.5
# 	shell:
# 		"if ($(cat {input} | awk '{{if ($2 >= {params.min_cvg} && $3 >= {params.min_perc}) print \"true\"; else print \"false\"}}') == true);"\
# 		"then ln -s ~/fiona/bacteremia/10.call_variants_v2/filtered_bam/{wildcards.sample}.{wildcards.organism}.filtered.bam {output}"\
# 		" && touch -h {output};"\
# 		"else touch {output}; fi"
		# "if ($(cat {input} | awk '{{if ($2 >= 5 && $3 >= 0.5) print \"true\"; else print \"false\"}}') == true); then touch {output}; fi"

rule faidx:
	input: "/home/tamburin/fiona/bacteremia/isolate_assemblies_all/original/{assembly}.fasta"
	output: "/home/tamburin/fiona/bacteremia/isolate_assemblies_all/original/{assembly}.fasta.fai"
	resources:
		mem=8,
		time=1
	shell:
		"samtools faidx {input}"

rule pileup:
	input:
		bam="filtered_bam/{sample}.{assembly}.filtered.bam",
		index=rules.faidx.output
	output: "pileup/{sample}.{assembly}.pileup"
	resources:
		mem=32,
		time=1
	threads: 1
	shell:
		"samtools mpileup -f /home/tamburin/fiona/bacteremia/14.assemble_trimmed/filtered_assemblies/{wildcards.assembly}.fasta -B -aa -o {output} {input.bam}"

rule call_snps:
	input: rules.pileup.output
	output: "snp_calls/{sample}.{assembly}.tsv"
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
	output: "tmp/{sample}.{assembly}.filtered.tmp"
	resources:
		mem=32,
		time=2
	threads: 1
	params:
		min_cvg=5,
		min_freq=0.8,
		min_phred=20
	shell:
		"(echo {wildcards.sample}.{wildcards.assembly}; cut -f4 {input}) > {output}"

rule core_snps:
	input:
		lambda wildcards: expand("tmp/{sample}.{assembly}.filtered.tmp", sample=organism_to_samples[samples_to_organism[wildcards.assembly]], assembly=wildcards.assembly)
	output: "core_snps/{assembly}.tsv"
	resources:
		mem=16,
		time=1
	threads: 1
	shell:
		"paste $(find tmp/*{wildcards.assembly}.filtered.tmp -size +1) | python3 snv_wf/scripts/findCoreSNPs.py > {output}"

rule core_snps_to_fasta:
	input: rules.core_snps.output
	output: "fasta/{assembly}.fasta"
	resources:
		mem=16,
		time=1
	threads: 1
	script:
		"scripts/coreSNPs2fasta.py"

rule multi_align:
	input: rules.core_snps_to_fasta.output
	output: "afa/{assembly}.afa"
	resources:
		mem=200,
		time=12
	threads: 1
	shell:
		"muscle -in {input} -out {output}"

rule build_tree:
	input: rules.multi_align.output
	output: "tree/{assembly}.tree"
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
	output: "tree_plot/{assembly}.tree.pdf"
	resources:
		mem=8,
		time=1
	threads: 1
	script:
		"scripts/renderTree.R"

rule pairwise_snvs:
	input: lambda wildcards: expand("tmp/{sample}.{assembly}.filtered.tmp", sample=organism_to_samples[samples_to_organism[wildcards.assembly]], assembly=wildcards.assembly)
	output: "snv_distances/{assembly}.snps.tsv"
	resources:
		mem=8,
		time=1
	threads: 1
	script:
		"scripts/pairwiseDist.py"
