#!/usr/bin/env python
import os
import sys
import re
import errno
import argparse


def parse_args(args=None):
    Description = 'Collapse LEFT/RIGHT primers in amplicon BED to single intervals.'
    Epilog = """Example usage: python collapse_amplicon_bed.py <FILE_IN> <FILE_OUT>"""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('FILE_IN', help="Input BED file.")
    parser.add_argument('FILE_OUT', help="Output BED file.")
    parser.add_argument('-lp', '--left_primer_suffix', type=str, dest="LEFT_PRIMER_SUFFIX", default='_LEFT', help="Suffix for left primer in name column of BED file (default: '_LEFT').")
    parser.add_argument('-rp', '--right_primer_suffix', type=str, dest="RIGHT_PRIMER_SUFFIX", default='_RIGHT', help="Suffix for right primer in name column of BED file (default: '_RIGHT').")
    return parser.parse_args(args)


def make_dir(path):
    if not len(path) == 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


## See https://stackoverflow.com/a/480227
def uniqify(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]


def collapse_amplicon_bed(FileIn,FileOut,LeftPrimerSuffix,RightPrimerSuffix):
    StartPosList = []
    IntervalDict = {}
    fin = open(FileIn,'r')
    while True:
        line = fin.readline()
        if line:
            chrom,start,end,name,score,strand = line.strip().split('\t')
            amplicon = re.sub(r'(?:{}|{}).*'.format(LeftPrimerSuffix,RightPrimerSuffix),'',name)
            if amplicon not in IntervalDict:
                IntervalDict[amplicon] = []
            IntervalDict[amplicon].append((chrom,int(start),int(end),score))
            StartPosList.append((int(start),amplicon))
        else:
            fin.close()
            break

    fout = open(FileOut,'w')
    for amplicon in uniqify([x[1] for x in sorted(StartPosList)]):
        posList = [item for elem in IntervalDict[amplicon] for item in elem[1:3]]
        chrom = IntervalDict[amplicon][0][0]
        start = min(posList)
        end = max(posList)
        strand = '+'
        score = IntervalDict[amplicon][0][3]
        fout.write(f'{chrom}\t{start}\t{end}\t{amplicon}\t{score}\t{strand}\n')
    fout.close()


def main(args=None):
    args = parse_args(args)
    collapse_amplicon_bed(args.FILE_IN,args.FILE_OUT,args.LEFT_PRIMER_SUFFIX,args.RIGHT_PRIMER_SUFFIX)


if __name__ == '__main__':
    sys.exit(main())
