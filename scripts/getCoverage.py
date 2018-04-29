# Fiona Tamburini
# Jan 26 2018

# calculate average coverage and % bases covered above min threshold from bedtools genomecov output
# usage python3 getCoverage.py 5 sample < depth.tsv

import sys

sample_file = snakemake.input[0]
out_file = snakemake.output[0]

minCvg = int(snakemake.params[0])

totalBases = 0
coveredBases = 0
weightedAvg = 0
with open(sample_file, 'r') as sample:
    for line in sample:
        if line.startswith("genome"):
            chr, depth, numBases, size, fraction = line.rstrip('\n').split('\t')

            depth = int(depth)
            numBases = int(numBases)

            totalBases += numBases
            weightedAvg += depth * numBases

            if depth >= minCvg:
                coveredBases += numBases
avgCvg = 0
percCovered = 0

if totalBases > 0:
    avgCvg = float(weightedAvg) / float(totalBases)
    percCovered = float(coveredBases) / float(totalBases)

#print("Total bases" + "\t" + str(totalBases))
#print("Average coverage:" + "\t" + str(round(avgCgv, 2)))
#print("% bases with at least " + str(minCvg) + "X coverage:" + "\t" + str(round(percCovered, 2)))
with open(out_file, 'w') as out:
    print("\t".join([str(round(avgCvg, 2)), str(round(percCovered, 2))]), file = out)
