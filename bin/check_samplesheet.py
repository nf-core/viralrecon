#!/usr/bin/env python

import os
import sys
import argparse

def parse_args(args=None):
    Description = 'Reformat nf-core/viralrecon samplesheet file and check its contents.'
    Epilog = """Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('FILE_IN', help="Input samplesheet file.")
    parser.add_argument('FILE_OUT', help="Output samplesheet file.")

    return parser.parse_args(args)


def print_error(error,line):
    print("ERROR: Please check samplesheet -> {}\nLine: '{}'".format(error,line.strip()))
    sys.exit(1)

def check_samplesheet(FileIn,FileOut):
    HEADER = ['sample', 'fastq_1', 'fastq_2']

    ## Check header
    fin = open(FileIn,'r')
    header = fin.readline().strip().split(',')
    if header != HEADER:
        print("ERROR: Please check samplesheet header -> {} != {}".format(','.join(header),','.join(HEADER)))
        sys.exit(1)

    sampleRunDict = {}
    while True:
        line = fin.readline()
        if line:
            lspl = [x.strip() for x in line.strip().split(',')]

            ## Check valid number of columns per row
            if len(lspl) != len(header):
                print_error("Invalid number of columns (minimum = {})!".format(len(header)),line)

            numCols = len([x for x in lspl if x])
            if numCols < 1:
                print_error("Invalid number of populated columns (minimum = 1)!",line)

            ## Cheack sample name entries
            sample,fastQFiles = lspl[0],lspl[1:]
            if sample:
                if sample.find(' ') != -1:
                    print_error("Sample entry contains spaces!",line)
            else:
                print_error("Sample entry has not been specified!",line)

            ## Check FastQ file extension
            for fastq in fastQFiles:
                if fastq:
                    if fastq.find(' ') != -1:
                        print_error("FastQ file contains spaces!",line)
                    if fastq[-9:] != '.fastq.gz' and fastq[-6:] != '.fq.gz':
                        print_error("FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",line)

            ## Auto-detect paired-end/single-end/is_sra
            single_end = '0'; is_sra = '0'
            fastq_1,fastq_2 = fastQFiles
            if sample and fastq_1 and fastq_2:              ## Paired-end short reads
                pass
            elif sample and fastq_1 and not fastq_2:        ## Single-end short reads
                single_end = '1'
            elif sample and not fastq_1 and not fastq_2:    ## SRA accession
                if sample[:3] == 'SRR':
                    is_sra = '1'
                else:
                    print_error("Please provide a valid SRA run accession starting with 'SRR'!",line)
            else:
                print_error("Invalid combination of columns provided!",line)

            sampleInfoList = [single_end, is_sra, fastq_1, fastq_2]
            if sample not in sampleRunDict:
                sampleRunDict[sample] = []
            else:
                if sampleInfoList in sampleRunDict[sample]:
                    print_error("Samplesheet contains duplicate rows!",line)
            sampleRunDict[sample].append(sampleInfoList)

        else:
            fin.close()
            break

    ## Write to file
    fout = open(FileOut,'w')
    fout.write(','.join(['sample_id', 'single_end', 'is_sra', 'fastq_1', 'fastq_2']) + '\n')
    for sample in sorted(sampleRunDict.keys()):

        ## Check that multiple runs of the same sample are of the same datatype
        if not all(x[:2] == sampleRunDict[sample][0][:2] for x in sampleRunDict[sample]):
            print_error("Multiple runs of a sample must be of the same datatype","Sample: {}".format(sample))

        for idx,val in enumerate(sampleRunDict[sample]):
            sample_id = "{}_T{}".format(sample,idx+1)
            fout.write(','.join([sample_id] + val) + ',\n')
    fout.close()


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN,args.FILE_OUT)


if __name__ == '__main__':
    sys.exit(main())
