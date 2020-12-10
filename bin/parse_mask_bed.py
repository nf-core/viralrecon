#!/usr/bin/env python
import os
import sys
import re
import argparse
import gzip


def parse_args(args=None):
    Description = 'Find indels positions in bed file'
    Epilog = """Example usage: python parse_mask_bed.py <BED_IN> <BED_OUT>"""
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('VCF_IN', help="Input vcf file.")
    parser.add_argument('BED_IN', help="Input bed file.")
    parser.add_argument('BED_OUT', help="Name of the output new bed file.")
    return parser.parse_args(args)

def find_indels_vcf(VcfIn):
    encoding = 'utf-8'
    indels_pos_len={}
    with gzip.open(VcfIn,'r') as f:
        for line in f:
            if '#' not in str(line, encoding):
                line = re.split("\t", str(line, encoding))
                var_pos=line[1]
                ref=line[3]
                alt=line[4]
                if len(alt) > len(ref):
                    indels_pos_len[var_pos] = len(alt)
                elif len(ref) > len(alt):
                    indels_pos_len[var_pos] = len(ref)
    return(indels_pos_len)

def find_dels_vcf(VcfIn):
    encoding = 'utf-8'
    dels_pos_len={}
    with gzip.open(VcfIn,'r') as f:
        for line in f:
            if '#' not in str(line, encoding):
                line = re.split("\t", str(line, encoding))
                var_pos=line[1]
                ref=line[3]
                alt=line[4]
                if len(ref) > len(alt):
                    dels_pos_len[var_pos] = len(ref)
    return(dels_pos_len)

def parse_mask_bed(BedIn,BedOut,indels_pos_len):
    fout = open(BedOut,'w')
    indels_postions=[]
    for pos in indels_pos_len:
        indels_postions.append(pos)
    with open(BedIn) as b:
        for line in b:
            line = re.split("\t", line)
            ref_genome=line[0]
            init_pos=line[1]
            end_pos=line[2]
            coverage=line[3]
            range_length=int(end_pos)-int(init_pos)
            oline=ref_genome+'\t'+init_pos+'\t'+end_pos+'\t'+coverage
            if len(indels_postions) > 0:
                for position in indels_postions:
                    if int(position) in range(int(init_pos),int(end_pos)):
                        indels_postions.remove(position)
                        indel_whole_length=indels_pos_len[position]
                        if end_pos == position:
                            new_end_pos=int(end_pos)-1
                            if int(new_end_pos) > int(init_pos):
                                oline=ref_genome+'\t'+str(new_init_pos)+'\t'+end_pos+'\t'+coverage
                                fout.write(oline)
                                break
                            else:
                                None #If the indel has size enought to cover all range, remove that line
                                break
                        else:
                            end_pos1=int(position)-1
                            init_pos2=int(position)+indel_whole_length-1
                            if int(end_pos1) > int(init_pos):
                                oline=ref_genome+'\t'+init_pos+'\t'+str(end_pos1)+'\t'+coverage
                                fout.write(oline)
                            if int(end_pos) > int(init_pos2):
                                oline=ref_genome+'\t'+str(init_pos2)+'\t'+end_pos+'\t'+coverage
                                fout.write(oline)
                            break
                    else:
                        oline=ref_genome+'\t'+init_pos+'\t'+end_pos+'\t'+coverage
                        fout.write(oline)
                        break
            else:
                oline=ref_genome+'\t'+init_pos+'\t'+end_pos+'\t'+coverage
                fout.write(oline)

########More def functions
def main(args=None):
    args = parse_args(args)
    dels_pos_len=find_dels_vcf(args.VCF_IN)
    parse_mask_bed(args.BED_IN, args.BED_OUT,dels_pos_len)

if __name__ == '__main__':
    sys.exit(main())
