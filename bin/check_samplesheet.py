#!/usr/bin/env python

import os
import sys
import csv
import requests
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


def get_sra_runinfo(sra_id):
    url = "https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=" + sra_id
    try:
        r = requests.get(url)
    except requests.exceptions.RequestException as e:
        raise SystemExit(e)
    if r.status_code != 200:
        print_error("ERROR: Connection to SRA failed\nError code '{}'".format(r.status_code))

    runInfoDict = {}
    csv_content = csv.DictReader(r.content.decode('utf-8').splitlines(), delimiter=',')
    for row in csv_content:
        runInfoDict[row['Run']] = row

    return runInfoDict


def check_samplesheet(FileIn,FileOut):
    HEADER = ['sample', 'fastq_1', 'fastq_2']

    ## Check header
    fin = open(FileIn,'r')
    header = fin.readline().strip().split(',')
    if header != HEADER:
        print("ERROR: Please check samplesheet header -> {} != {}".format(','.join(header),','.join(HEADER)))
        sys.exit(1)

    sraRunInfoDict = {}
    sraWarningList = []
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
            is_sra = False
            sampleInfoList = []                             ## [sample, single_end, is_sra, fastq_1, fastq_2]
            fastq_1,fastq_2 = fastQFiles
            if sample and fastq_1 and fastq_2:                              ## Paired-end short reads
                sampleInfoList = [[sample, '0', '0', fastq_1, fastq_2]]
            elif sample and fastq_1 and not fastq_2:                        ## Single-end short reads
                sampleInfoList = [[sample, '1', '0', fastq_1, fastq_2]]
            elif sample and not fastq_1 and not fastq_2:    ## SRA accession
                if sample[:2] in ['SR', 'PR']:
                    is_sra = True
                    runInfoDict = get_sra_runinfo(sample)
                    if len(runInfoDict) != 0:
                        for sra_id in runInfoDict.keys():
                            platform = runInfoDict[sra_id]['Platform']
                            library = runInfoDict[sra_id]['LibraryLayout']
                            if runInfoDict[sra_id]['Platform'].lower() in ['illumina']:
                                if library.lower() == 'single':
                                    sampleInfoList = [[sra_id, '1', '1', '', '']]
                                    single_end = '1'; fastq_1 = sra_id + "_1.fastq.gz"
                                elif library.lower() == 'paired':
                                    sampleInfoList = [[sra_id, '0', '1', '', '']]
                                else:
                                    sraWarningList.append("WARNING: Library layout '{}' != 'SINGLE' or 'PAIRED' for SRA id ({},{})!".format(layout,sample,sra_id))
                            else:
                                sraWarningList.append("WARNING: Illumina platform currently supported. SRA id ({},{}) was sequenced on the '{}' platform!".format(sample,sra_id,platform))
                    else:
                        sraWarningList.append("WARNING: No data available for SRA id {}!".format(sample))
                else:
                    print_error("Please provide a valid SRA id starting with 'SR' or 'PR'!",line)
            else:
                print_error("Invalid combination of columns provided!",line)

            for sampleList in sampleInfoList:
                sample_id = sampleList[0]; infoList = sampleList[1:]
                if sample_id not in sampleRunDict:
                    sampleRunDict[sample_id] = []
                else:
                    if infoList in sampleRunDict[sample_id]:
                        if is_sra:
                            sraWarningList.append("WARNING: Duplicate SRA ids observed for ({},{})!".format(sample,sample_id))
                        else:
                            print_error("Samplesheet contains duplicate rows!",line)
                sampleRunDict[sample_id].append(infoList)

        else:
            fin.close()
            break

    for i in sraWarningList:
        print i

    ## Write to file
    fout = open(FileOut,'w')
    fout.write(','.join(['sample_id', 'single_end', 'is_sra', 'fastq_1', 'fastq_2']) + '\n')
    for sample in sorted(sampleRunDict.keys()):

        ## Check that multiple runs of the same sample are of the same datatype
        if not all(x[:2] == sampleRunDict[sample][0][:2] for x in sampleRunDict[sample]):
            print_error("Multiple runs of a sample must be of the same datatype","Sample: {}".format(sample))

        for idx,val in enumerate(sampleRunDict[sample]):
            sample_id = sample
            if val[1] == '0':
                sample_id = "{}_T{}".format(sample_id,idx+1)
            fout.write(','.join([sample_id] + val) + ',\n')
    fout.close()


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN,args.FILE_OUT)


if __name__ == '__main__':
    sys.exit(main())
