# find core snvs shared between samples

import sys

#genotype_file = sys.argv[1]

# counter = 0

#with open(genotype_file) as genotypes:
#    for line in genotypes:
for line in sys.stdin:
    positions = line.rstrip('\n').split('\t')

    # if position is undefined for any sample or if the positions are the same in every file
    if (len(set(positions)) > 1) and ("N" not in positions):
        print("\t".join(positions))
#        if "N" not in positions:
#            counter +=1

#print(counter)
