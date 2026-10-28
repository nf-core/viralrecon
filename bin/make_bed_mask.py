#!/usr/bin/env python

import os
import sys
import re
import argparse
import gzip


def parse_args(args=None):
    Description = "Find indels positions in bed file"
    Epilog = "Example usage: python make_bed_mask.py <VCF_IN> <BED_IN> <BED_OUT>"
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("VCF_IN", help="Input vcf file.")
    parser.add_argument("BED_IN", help="Input bed file.")
    parser.add_argument("BED_OUT", help="Name of the output new bed file.")
    return parser.parse_args(args)


def find_indels_vcf(vcf_in):
    encoding = "utf-8"
    indels_pos_len = {}
    with gzip.open(vcf_in, "r") as f:
        for line in f:
            if "#" not in str(line, encoding):
                line = re.split("\t", str(line, encoding))
                var_pos = line[1]
                ref = line[3]
                alt = line[4]
                if len(ref) != len(alt):
                    indels_pos_len[var_pos] = len(ref)
    return indels_pos_len


def make_bed_mask(bed_in, bed_out, indels_pos_len):
    fout = open(bed_out, "w")
    indels_positions = []
    for pos in indels_pos_len:
        indels_positions.append(pos)
    with open(bed_in) as b:
        for line in b:
            line = re.split("\t", line)
            ref_genome = line[0]
            init_pos = line[1]
            end_pos = line[2]
            range_length = int(end_pos) - int(init_pos)
            oline = ref_genome + "\t" + init_pos + "\t" + end_pos
            test = True
            for position in indels_positions:
                indel_init_pos = position
                indel_whole_length = indels_pos_len[position]
                indel_end_pos = int(indel_init_pos) + int(indel_whole_length) - 1
                if int(init_pos) in range(int(indel_init_pos), int(indel_end_pos)) or int(end_pos) in range(
                    int(indel_init_pos), int(indel_end_pos)
                ):
                    test = False
                    break
                else:
                    oline = ref_genome + "\t" + init_pos + "\t" + end_pos
            if test:
                fout.write(oline + "\n")


def main(args=None):
    args = parse_args(args)
    indels_pos_len = find_indels_vcf(args.VCF_IN)
    make_bed_mask(args.BED_IN, args.BED_OUT, indels_pos_len)


if __name__ == "__main__":
    sys.exit(main())
