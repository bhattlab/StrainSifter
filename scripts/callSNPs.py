#!/usr/bin/env python3
# callSNPS.py
# call SNPs from a samtools pileup file

import sys
import gzip
import re

DEBUG = 0

ASCII_OFFSET = ord("!")
INDEL_PATTERN = re.compile(r"[+-](\d+)")
INDEL_STRING = "[+-]INTEGER[ACGTN]{INTEGER}"
ROUND = 3  # places to right of decimal point

min_coverage, min_proportion, min_qual = snakemake.params
pileup_file = snakemake.input[0]
out_file = snakemake.output[0]

min_coverage = int(min_coverage)
min_proportion = float(min_proportion)
min_qual = int(min_qual)

## function to process each line of pileup

def parse_pileup(line):

    # read fields from line of pileup
    chromosome, position, reference, coverage, pileup, quality = line.rstrip("\r\n").split("\t")

    # proportion is 0 if no consensus base, top base is N by default
    parsed = {
        "proportion": 0.0,
        "chromosome": chromosome,
        "position": int(position),
        "reference": reference,
        "coverage": int(coverage),
        "pileup": pileup,
        "quality": quality,
        "top_base": "N"}

    # if the base coverage is below the limit or above our acceptable max, call N
    if parsed["coverage"] < min_coverage:
        return parsed

    # uppercase pileup string for processing
    pileup = pileup.upper()

    # Remove start and stop characters from pileup string
    pileup = re.sub(r"\^.|\$", "", pileup)

    # Remove indels from pileup string
    start = 0

    while True:
        match = INDEL_PATTERN.search(pileup, start)

        if match:
            integer = match.group(1)
            pileup = re.sub(INDEL_STRING.replace("INTEGER", integer), "", pileup)
            start = match.start()
        else:
            break

    # get total base count and top base count
    total = 0
    top_base = "N"
    top_base_count = 0

    # uppercase reference base for comparison
    reference = reference.upper()

    base_counts = {"A": 0, "C": 0, "G": 0, "T": 0}

    quality_length = len(quality)

    for i in range(quality_length):

        # convert ASCII character to phred base quality
        base_quality = ord(quality[i]) - ASCII_OFFSET

        # only count high-quality bases
        if base_quality >= min_qual:

            currBase = pileup[i]

            if currBase in base_counts:
                base_counts[currBase] += 1

            else:
                base_counts[reference] += 1

            total += 1

    parsed["total"] = total

    # find top base
    for base in base_counts:
        if base_counts[base] > top_base_count:
            top_base = base
            top_base_count = base_counts[base]

    # if more that 0 bases processed
    if total > 0:
        prop = top_base_count / float(total)

        if prop >= min_proportion:
            parsed["proportion"] = prop
            parsed["top_base"] = top_base

    return parsed

## process pileup file

# read in input file
if not DEBUG:
    out_file = open(out_file, "w")

with open(pileup_file, "rt") as pileup_file:
    for pileup_file_line in pileup_file:

        parsed = parse_pileup(pileup_file_line)

        proportion = parsed["proportion"]

        # print to output file
        if DEBUG:
            print("\t".join([parsed["chromosome"], str(parsed["position"]), parsed["reference"], parsed["top_base"], str(round(proportion, ROUND)),
                parsed["pileup"], parsed["quality"]]))
        else:
            print("\t".join([parsed["chromosome"], str(parsed["position"]), parsed["reference"], parsed["top_base"], str(round(proportion, ROUND)),
                parsed["pileup"], parsed["quality"]]), file=out_file)

if not DEBUG:
    out_file.close()
