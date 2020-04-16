#!/usr/bin/env python
import os
import sys
import re
import argparse

def parse_args(args=None):
    Description = 'Reformat iVar tsv files to vcf files.'
    Epilog = """Example usage: python tsv_to_vcf.py <FILE_IN> <FILE_OUT>"""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('FILE_IN', help="Input tsv file.")
    parser.add_argument('FILE_OUT', help="Output file with .vfc.")
    return parser.parse_args(args)

def parse_tsv(FileIn,OutPrefix):
    file_basename=os.path.splitext(FileIn)[0]
    header_1='##fileformat=VCFv4.2\n##source=iVar\n##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n##FILTER=<ID=PASS,Description="Result of p-value <= 0.05">\n##FILTER=<ID=nopass,Description="Result of p-value > 0.05">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=REF_DP,Number=1,Type=Integer,Description="Depth of reference base">\n##FORMAT=<ID=REF_RV,Number=1,Type=Integer,Description="Depth of reference base on reverse reads">\n##FORMAT=<ID=REF_QUAL,Number=1,Type=Integer,Description="Mean quality of reference base">\n##FORMAT=<ID=ALT_DP,Number=1,Type=Integer,Description="Depth of alternate base">\n##FORMAT=<ID=ALT_RV,Number=1,Type=Integer,Description="Deapth of alternate base on reverse reads">\n##FORMAT=<ID=ALT_QUAL,Number=1,Type=String,Description="Mean quality of alternate base">\n##FORMAT=<ID=ALT_FREQ,Number=1,Type=String,Description="Frequency of alternate base">\n'
    header_2='#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+file_basename+'\n'

    o = open(OutPrefix,'w')
    o.write(header_1)
    o.write(header_2)
    with open(FileIn) as f:
        for line in f:
            if not re.match("REGION",line):
                line = re.split("\t", line)
                CHROM=line[0]
                POS=line[1]
                ID='.'
                REF=line[2]
                ALT=line[3]
                QUAL='.'
                pass_test=line[13]
                if pass_test == 'TRUE':
                    FILTER='PASS'
                else:
                    FILTER='nopass'
                INFO='DP='+line[11]
                FORMAT='GT:REF_DP:REF_RV:REF_QUAL:ALT_DP:ALT_RV:ALT_QUAL:ALT_FREQ'
                SAMPLE='1:'+line[4]+':'+line[5]+':'+line[6]+':'+line[7]+':'+line[8]+':'+line[9]+':'+line[10]
                line = CHROM+'\t'+POS+'\t'+ID+'\t'+REF+'\t'+ALT+'\t'+QUAL+'\t'+FILTER+'\t'+INFO+'\t'+FORMAT+'\t'+SAMPLE+'\n'
                o.write(line)
    o.close()

def main(args=None):
    args = parse_args(args)
    parse_tsv(args.FILE_IN,args.FILE_OUT)


if __name__ == '__main__':
    sys.exit(main())
