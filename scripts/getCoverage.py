# Fiona Tamburini
# Jan 26 2018

# calculate average coverage and % bases covered above min threshold from bedtools genomecov output
# usage python3 getCoverage.py 5 sample < depth.tsv

import sys

minCvg = int(sys.argv[1])
sample = sys.argv[2]

totalBases = 0
coveredBases = 0
weightedAvg = 0

for line in sys.stdin:
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

print("\t".join([sample, str(round(avgCvg, 2)), str(round(percCovered, 2))]))
