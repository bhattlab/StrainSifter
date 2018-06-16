import os
import re
# configfile: "strainsifter/config.yaml"

localrules: do_filter

# read list of samples
samples = []
with open(config['samples']) as samples_list:
	for line in samples_list:
		sample = line.rstrip('\n')
		samples += [sample]

passed_samples = []
try:
	f = open("passed_samples.list")
	for line in f:
		sample = line.rstrip('\n')
		passed_samples += [sample]
except (IOError, OSError) as e:
	print("File {f} not yet created".format(f="passed_samples.list"))
else:
	f.close()

refname = re.split("/|\.", config['reference'])[-2]

if config['paired_end'] == "Y":
	reads = expand("{reads}/{{sample}}_{pe}.fq.gz", pe=["1", "2"], reads=config['reads_dir'])
else:
	reads = "{reads}.fq".format(reads=config['reads_dir'])

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
		r = reads
	output:
		"filtered_bam/{sample}.filtered.bam"
	resources:
		mem=32,
		time=6
	threads: 8
	params:
		qual=40
	shell:
		"bwa mem -t {threads} {input.ref} {input.r} | "\
		"samtools view -b -q {params.qual} | bamtools filter -tag 'NM:<6' | "\
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
		cvg=5
	script:
		"scripts/getCoverage.py"

rule filter_samples:
	input: expand("coverage/{sample}.cvg", sample = samples)
	output: dynamic("passed_samples/{sample}.bam")
		# "passed_samples/{sample}.bam"
	resources:
		mem=1,
		time=1
	threads: 1
	params:
		min_cvg=config['min_cvg'],
		min_perc=0.5
	run:
		samps = os.listdir("coverage")
		for samp in samps:
			with open("coverage/" + samp) as s:
				cvg, perc = s.readline().rstrip('\n').split('\t')
			if (float(cvg) >= params.min_cvg and float(perc) > params.min_perc):
				# passed = "passed_samples/" + samp.rstrip(".cvg") + ".bam"
				shell("ln -s $PWD/filtered_bam/{p}.filtered.bam passed_samples/{p}.bam; "\
				"echo {p} >> passed_samples.list".format(p=samp.rstrip(".cvg")))
	# shell:
	# 	"if ($(cat {input} | awk '{{if ($1 >= {params.min_cvg} && $2 >= {params.min_perc}) print \"true\"; else print \"false\"}}') == true);"\
	# 	"then ln -s $PWD/filtered_bam/{wildcards.sample}.filtered.bam {output}"\
	# 	" && touch -h {output}"\
		# "else touch {output}; fi"
rule do_filter:
	input: rules.filter_samples.output
	output: "filtered"
	shell:
		"touch {output}"

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
		bam="passed_samples/{sample}.bam",
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

rule combine:
	input:
		samps=expand("consensus/{sample}.txt", sample = passed_samples),
		filt=rules.filter_samples.output
	output: "{ref}.cns.tsv".format(ref=refname)
	resources:
		mem=2,
		time=1
	threads: 1
	shell:
		"paste {input.samps} > {output}"

rule core_snps:
	input: rules.combine.output
	output: "{ref}.tsv".format(ref=refname)
	resources:
		mem=16,
		time=1
	threads: 1
	script:
		"scripts/findCoreSNPs.py"

rule core_snps_to_fasta:
	input: rules.core_snps.output
	output: "{ref}.fasta".format(ref=refname)
	resources:
		mem=16,
		time=1
	threads: 1
	script:
		"scripts/coreSNPs2fasta.py"

rule multi_align:
	input: rules.core_snps_to_fasta.output
	output: "{ref}.afa".format(ref=refname)
	resources:
		mem=200,
		time=12
	threads: 1
	shell:
		"muscle -in {input} -out {output}"

rule build_tree:
	input: rules.multi_align.output
	output: "{ref}.tree".format(ref=refname)
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
	output: "{ref}.tree.pdf".format(ref=refname)
	resources:
		mem=8,
		time=1
	threads: 1
	script:
		"scripts/renderTree.R"

rule pairwise_snvs:
	input: expand("consensus/{sample}.txt", sample = passed_samples)
	output: "{ref}.dist.tsv".format(ref=refname)
	resources:
		mem=8,
		time=1
	threads: 1
	script:
		"scripts/pairwiseDist.py"
