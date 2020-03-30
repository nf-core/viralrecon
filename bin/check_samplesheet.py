#!/usr/bin/env python

import os
import sys
import argparse

def parse_args(args=None):
    Description = 'Reformat nf-core/artic samplesheet file and check its contents.'
    Epilog = """Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('FILE_IN', help="Input samplesheet file.")
    parser.add_argument('FILE_OUT', help="Output samplesheet file.")

    return parser.parse_args(args)


def print_error(error,line):
    print("ERROR: Please check samplesheet -> {}\nLine: '{}'".format(error,line.strip()))


def check_samplesheet(FileIn,FileOut):
    HEADER = ['sample', 'run', 'short_fastq_1', 'short_fastq_2', 'long_fastq']

    ## CHECK HEADER
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

            ## CHECK VALID NUMBER OF COLUMNS PER SAMPLE
            if len(lspl) != len(header):
                print_error("Invalid number of columns (minimum = {})!".format(len(header)),line)
                sys.exit(1)

            numCols = len([x for x in lspl if x])
            if numCols < 3:
                print_error("Invalid number of populated columns (minimum = 3)!",line)
                sys.exit(1)

            ## CHECK SAMPLE ID ENTRIES
            sample,run,fastQFiles = lspl[0],lspl[1],lspl[2:]
            if sample:
                if sample.find(' ') != -1:
                    print_error("Sample entry contains spaces!",line)
                    sys.exit(1)
            else:
                print_error("Sample entry has not been specified!",line)
                sys.exit(1)

            ## CHECK RUN COLUMN IS INTEGER
            if not run.isdigit():
                print_error("Run id not an integer!",line)
                sys.exit(1)

            ## CHECK FASTQ FILE EXTENSION
            for fastq in fastQFiles:
                if fastq:
                    if fastq.find(' ') != -1:
                        print_error("FastQ file contains spaces!",line)
                        sys.exit(1)
                    if fastq[-9:] != '.fastq.gz' and fastq[-6:] != '.fq.gz':
                        print_error("FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",line)
                        sys.exit(1)

            ## AUTO-DETECT ILLUMINA/NANOPORE
            readDict = {}
            short_fastq_1,short_fastq_2,long_fastq = fastQFiles

            ## Paired-end short reads only
            if short_fastq_1 and short_fastq_2 and not long_fastq:
                readDict[sample+'_SR'] = ['0', '0', short_fastq_1, short_fastq_2]  ## [ SINGLE_END?, LONG_READS?, FASTQ_1, FASTQ_2 ]

            ## Paired-end short reads and long reads
            elif short_fastq_1 and short_fastq_2 and long_fastq:
                readDict[sample+'_SR'] = ['0', '0', short_fastq_1, short_fastq_2]
                readDict[sample+'_LR'] = ['0', '1', long_fastq, '']

            ## Single-end short reads only
            elif short_fastq_1 and not short_fastq_2 and not long_fastq:
                readDict[sample+'_SR'] = ['1', '0', short_fastq_1, '']

            ## Single-end short reads and long reads
            elif short_fastq_1 and not short_fastq_2 and long_fastq:
                readDict[sample+'_SR'] = ['1', '0', short_fastq_1, '']
                readDict[sample+'_LR'] = ['0', '1', long_fastq, '']

            elif not short_fastq_1 and not short_fastq_2 and long_fastq:    ## Long reads only
                readDict[sample+'_LR'] = ['0', '1', long_fastq, '']

            else:
                print_error("'short_fastq_2' cannot be specified without 'short_fastq_1'!",line)
                sys.exit(1)

            ## CREATE SAMPLE MAPPING DICT = {SAMPLE_ID: {RUN_ID:[ SINGLE_END, LONG_READS, FASTQ_1, FASTQ_2 ]}
            run = int(run)
            for rsample in readDict.keys():
                if rsample not in sampleRunDict:
                    sampleRunDict[rsample] = {}
                if run not in sampleRunDict[rsample]:
                    sampleRunDict[rsample][run] = readDict[rsample]
                else:
                    print_error("Duplicate run IDs found!",line)
                    sys.exit(1)

        else:
            fin.close()
            break

    ## WRITE TO FILE
    fout = open(FileOut,'w')
    fout.write(','.join(['sample_id', 'single_end', 'long_reads', 'fastq_1', 'fastq_2']) + '\n')
    for sample in sorted(sampleRunDict.keys()):

        ## CHECK THAT RUN IDS ARE IN FORMAT 1..<NUM_RUNS>
        run_ids = set(sampleRunDict[sample].keys())
        if len(run_ids) != max(run_ids):
            print_error("Run IDs must start with 1..<num_runs>!","Sample: {}, Run IDs: {}".format(sample,list(run_ids)))
            sys.exit(1)

        ## CHECK THAT MULTIPLE RUNS ARE FROM THE SAME DATATYPE
        if not all(x[:2] == list(sampleRunDict[sample].values())[0][:2] for x in list(sampleRunDict[sample].values())):
            print_error("Multiple runs of a sample must be of the same datatype","Sample: {}, Run IDs: {}".format(sample,list(run_ids)))
            sys.exit(1)

        for run in sorted(sampleRunDict[sample].keys()):
            sample_id = "{}_T{}".format(sample,run)
            fout.write(','.join([sample_id] + sampleRunDict[sample][run]) + ',\n')
    fout.close()


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN,args.FILE_OUT)


if __name__ == '__main__':
    sys.exit(main())
