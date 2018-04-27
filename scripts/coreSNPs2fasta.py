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
