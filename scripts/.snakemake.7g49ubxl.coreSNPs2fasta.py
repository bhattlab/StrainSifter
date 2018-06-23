
######## Snakemake header ########
import sys; sys.path.append("/home/tamburin/tools/miniconda3/lib/python3.6/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x05\x00\x00\x00inputq\x03csnakemake.io\nInputFiles\nq\x04)\x81q\x05X\x14\x00\x00\x00E_coli.core_snps.tsvq\x06a}q\x07X\x06\x00\x00\x00_namesq\x08}q\tsbX\x06\x00\x00\x00outputq\ncsnakemake.io\nOutputFiles\nq\x0b)\x81q\x0cX\x0c\x00\x00\x00E_coli.fastaq\ra}q\x0eh\x08}q\x0fsbX\x06\x00\x00\x00paramsq\x10csnakemake.io\nParams\nq\x11)\x81q\x12}q\x13h\x08}q\x14sbX\t\x00\x00\x00wildcardsq\x15csnakemake.io\nWildcards\nq\x16)\x81q\x17}q\x18h\x08}q\x19sbX\x07\x00\x00\x00threadsq\x1aK\x01X\t\x00\x00\x00resourcesq\x1bcsnakemake.io\nResources\nq\x1c)\x81q\x1d(K\x01K\x01K\x10K\x01e}q\x1e(h\x08}q\x1f(X\x06\x00\x00\x00_coresq K\x00N\x86q!X\x06\x00\x00\x00_nodesq"K\x01N\x86q#X\x03\x00\x00\x00memq$K\x02N\x86q%X\x04\x00\x00\x00timeq&K\x03N\x86q\'uh K\x01h"K\x01h$K\x10h&K\x01ubX\x03\x00\x00\x00logq(csnakemake.io\nLog\nq))\x81q*}q+h\x08}q,sbX\x06\x00\x00\x00configq-}q.(X\t\x00\x00\x00referenceq/X\x0e\x00\x00\x00isolate2.fastaq0X\x05\x00\x00\x00readsq1}q2(X\x05\x00\x00\x00meta1q3]q4(X\x16\x00\x00\x00reads/meta1_1.fastq.gzq5X\x16\x00\x00\x00reads/meta1_1.fastq.gzq6eX\x05\x00\x00\x00meta2q7]q8(X\x16\x00\x00\x00reads/meta2_1.fastq.gzq9X\x16\x00\x00\x00reads/meta2_1.fastq.gzq:eX\x05\x00\x00\x00meta3q;]q<(X\x16\x00\x00\x00reads/meta3_1.fastq.gzq=X\x16\x00\x00\x00reads/meta3_1.fastq.gzq>eX\x08\x00\x00\x00isolate1q?]q@(X\x19\x00\x00\x00reads/isolate1_1.fastq.gzqAX\x19\x00\x00\x00reads/isolate1_2.fastq.gzqBeX\x08\x00\x00\x00isolate2qC]qD(X\x19\x00\x00\x00reads/isolate2_1.fastq.gzqEX\x19\x00\x00\x00reads/isolate2_2.fastq.gzqFeuX\x06\x00\x00\x00prefixqGX\x06\x00\x00\x00E_coliqHX\x04\x00\x00\x00mapqqIK<X\x0c\x00\x00\x00n_mismatchesqJK\x05X\x07\x00\x00\x00min_cvgqKK\x05X\x12\x00\x00\x00min_genome_percentqLG?\xe0\x00\x00\x00\x00\x00\x00X\t\x00\x00\x00base_freqqMG?\xe9\x99\x99\x99\x99\x99\x9auX\x04\x00\x00\x00ruleqNX\x12\x00\x00\x00core_snps_to_fastaqOub.'); from snakemake.logging import logger; logger.printshellcmds = False
######## Original script #########
import sys

HEADER_POS = 0
FASTA_POS = 1
#coreSNPs_file = sys.argv[1]

snp_file = snakemake.input[0]
out_file = snakemake.output[0]

# hold sequences as we go
fastas = []

with open(snp_file, "r") as core_snps:

    # process header
    header = core_snps.readline().rstrip("\n").split("\t")

    # add header to dict with empty string
    for h in range(len(header)):
        fastas += [[header[h], ""]]

    for line in core_snps:
        line = line.rstrip("\n").split("\t")

        for pos in range(len(line)):
            # print(fastas[pos])
            # print(fastas[pos][FASTA_POS])
            fastas[pos][FASTA_POS] += line[pos]

with open(out_file, "w") as sys.stdout:
    for f in range(len(fastas)):
        print(">" + fastas[f][0])
        print(fastas[f][1])
