#!/usr/bin/env python3
# find core snvs shared between samples

import sys

snp_file = snakemake.input[0]
out_file = snakemake.output[0]

# read each line and print if base is called for every sample and
# bases are not the same in every sample
with open(snp_file, "r") as snps:
    with open(out_file, "w") as sys.stdout:

        # print header
        print(snps.readline().rstrip('\n'))

        # process each line and output positions with no unknown bases ("N")
        # and where all samples do not have same base
        for line in snps:
            positions = line.rstrip('\n').split('\t')
            if (len(set(positions)) > 1) and ("N" not in positions):
                print("\t".join(positions))
