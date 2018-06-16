
######## Snakemake header ########
import sys; sys.path.append("/home/tamburin/tools/miniconda3/lib/python3.6/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x05\x00\x00\x00inputq\x03csnakemake.io\nInputFiles\nq\x04)\x81q\x05X\x1c\x00\x00\x00pileup/4308_7-26-2016.pileupq\x06a}q\x07X\x06\x00\x00\x00_namesq\x08}q\tsbX\x06\x00\x00\x00outputq\ncsnakemake.io\nOutputFiles\nq\x0b)\x81q\x0cX\x1c\x00\x00\x00snp_calls/4308_7-26-2016.tsvq\ra}q\x0eh\x08}q\x0fsbX\x06\x00\x00\x00paramsq\x10csnakemake.io\nParams\nq\x11)\x81q\x12(K\x05G?\xe9\x99\x99\x99\x99\x99\x9aK\x14e}q\x13(h\x08}q\x14(X\x07\x00\x00\x00min_cvgq\x15K\x00N\x86q\x16X\x08\x00\x00\x00min_freqq\x17K\x01N\x86q\x18X\x08\x00\x00\x00min_qualq\x19K\x02N\x86q\x1auh\x15K\x05h\x17G?\xe9\x99\x99\x99\x99\x99\x9ah\x19K\x14ubX\t\x00\x00\x00wildcardsq\x1bcsnakemake.io\nWildcards\nq\x1c)\x81q\x1dX\x0e\x00\x00\x004308_7-26-2016q\x1ea}q\x1f(h\x08}q X\x06\x00\x00\x00sampleq!K\x00N\x86q"sX\x06\x00\x00\x00sampleq#h\x1eubX\x07\x00\x00\x00threadsq$K\x01X\t\x00\x00\x00resourcesq%csnakemake.io\nResources\nq&)\x81q\'(K\x01K\x01K K\x02e}q((h\x08}q)(X\x06\x00\x00\x00_coresq*K\x00N\x86q+X\x06\x00\x00\x00_nodesq,K\x01N\x86q-X\x03\x00\x00\x00memq.K\x02N\x86q/X\x04\x00\x00\x00timeq0K\x03N\x86q1uh*K\x01h,K\x01h.K h0K\x02ubX\x03\x00\x00\x00logq2csnakemake.io\nLog\nq3)\x81q4}q5h\x08}q6sbX\x06\x00\x00\x00configq7}q8(X\t\x00\x00\x00referenceq9XO\x00\x00\x00/home/tamburin/fiona/bacteremia/1.assemble_trimmed/filtered_assemblies/L3.fastaq:X\x07\x00\x00\x00samplesq;X0\x00\x00\x00/home/tamburin/fiona/crassphage/all_samples.listq<X\t\x00\x00\x00reads_dirq=X)\x00\x00\x00/home/tamburin/fiona/crassphage/readlinksq>X\n\x00\x00\x00paired_endq?X\x01\x00\x00\x00Yq@X\x07\x00\x00\x00min_cvgqAK\nuX\x04\x00\x00\x00ruleqBX\t\x00\x00\x00call_snpsqCub.'); from snakemake.logging import logger; logger.printshellcmds = False
######## Original script #########
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
        "proportion": 0.0, "chromosome": chromosome,
        "position": int(position), "reference": reference,
        "coverage": int(coverage), "pileup": pileup,
        "quality": quality, "top_base": "N"}

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
