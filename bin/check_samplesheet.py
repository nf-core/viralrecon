#!/usr/bin/env python

import os
import sys
import csv
import requests
import argparse

def parse_args(args=None):
    Description = 'Reformat nf-core/viralrecon samplesheet file and check its contents.'
    Epilog = """Example usage: python check_samplesheet.py <FILE_IN> <OUT_PREFIX>"""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('FILE_IN', help="Input samplesheet file.")
    parser.add_argument('OUT_PREFIX', help="Output file prefix.")
    parser.add_argument('-ss', '--skip_sra', dest="SKIP_SRA", help="Skip steps involving SRA identifiers.",action='store_true')
    parser.add_argument('-is', '--ignore_sra_errors', dest="IGNORE_SRA_ERRORS", help="Ignore SRA validation errors as opposed to failing.",action='store_true')

    return parser.parse_args(args)

def make_dir(path):
    if not len(path) == 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

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
        print("ERROR: Connection to SRA failed\nError code '{}'".format(r.status_code))
        sys.exit(1)

    runInfoDict = {}
    csv_content = csv.DictReader(r.content.decode('utf-8').splitlines(), delimiter=',')
    for row in csv_content:
        runInfoDict[row['Run']] = row

    return runInfoDict


def check_samplesheet(FileIn,OutPrefix,ignoreSRAErrors=False,skipSRA=False):
    HEADER = ['sample', 'fastq_1', 'fastq_2']

    ## Check header
    fin = open(FileIn,'r')
    header = fin.readline().strip().split(',')
    if header != HEADER:
        print("ERROR: Please check samplesheet header -> {} != {}".format(','.join(header),','.join(HEADER)))
        sys.exit(1)

    sampleRunDict = {}
    sraRunInfoHeader = []
    sraRunInfoDict = {}
    sraWarningList = []
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
            sampleInfoList = []                                             ## [sample, single_end, is_sra, fastq_1, fastq_2]
            fastq_1,fastq_2 = fastQFiles
            if sample and fastq_1 and fastq_2:                              ## Paired-end short reads
                sampleInfoList = [[sample, '0', '0', fastq_1, fastq_2]]

            elif sample and fastq_1 and not fastq_2:                        ## Single-end short reads
                sampleInfoList = [[sample, '1', '0', fastq_1, fastq_2]]

            elif sample and not fastq_1 and not fastq_2:                    ## SRA accession
                if not skipSRA:
                    if sample[:2] in ['SR', 'PR']:
                        is_sra = True
                        runInfoDict = get_sra_runinfo(sample)
                        if len(runInfoDict) != 0:
                            for sra_id in runInfoDict.keys():
                                platform = runInfoDict[sra_id]['Platform']
                                library = runInfoDict[sra_id]['LibraryLayout']
                                if runInfoDict[sra_id]['Platform'].lower() in ['illumina']:
                                    if library.lower() == 'single':
                                        sampleInfoList.append([sra_id, '1', '1', '', ''])
                                    elif library.lower() == 'paired':
                                        sampleInfoList.append([sra_id, '0', '1', '', ''])
                                    else:
                                        errorStr = "Library layout '{}' != 'SINGLE' or 'PAIRED' for SRR id:'{}'! User provided id:'{}'.".format(layout,sra_id,sample)
                                        if not ignoreSRAErrors:
                                            print_error(errorStr,line)
                                        sraWarningList.append(errorStr)
                                else:
                                    errorStr = "Only Illumina platform is currently supported. SRR id:'{}' is from the '{}' platform! User provided id:'{}'.".format(sra_id,platform,sample)
                                    if not ignoreSRAErrors:
                                        print_error(errorStr,line)
                                    sraWarningList.append(errorStr)

                                sraRunInfoHeader = list(set(sraRunInfoHeader) | set(runInfoDict[sra_id].keys()))
                                sraRunInfoDict[sra_id] = runInfoDict[sra_id]
                        else:
                            errorStr = "No data available for SRA id:{}!".format(sample)
                            if not ignoreSRAErrors:
                                print_error(errorStr,line)
                            sraWarningList.append(errorStr)
                    else:
                        print_error("Please provide a valid SRA id starting with 'SR' or 'PR'!",line)
            else:
                print_error("Invalid combination of columns provided!",line)

            for sampleList in sampleInfoList:
                sample_id = sampleList[0]; infoList = sampleList[1:]
                if sample_id not in sampleRunDict:
                    sampleRunDict[sample_id] = [infoList]
                else:
                    if is_sra:
                        errorStr = "Duplicate SRR id:'{}'! User provided id:'{}'.".format(sample_id,sample)
                        sraWarningList.append(errorStr)
                    else:
                        if infoList in sampleRunDict[sample_id]:
                            print_error("Samplesheet contains duplicate rows!",line)
                        else:
                            sampleRunDict[sample_id].append(infoList)

        else:
            fin.close()
            break

    ## Write validated samplesheet with appropriate columns
    OutDir = os.path.dirname(OutPrefix)
    make_dir(OutDir)
    fout = open(os.path.join(OutDir,'{}.csv'.format(OutPrefix)),'w')
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

    ## Write SRA warnings to file
    if len(sraWarningList) != 0:
        fout = open(os.path.join(OutDir,'{}.sra_warnings.txt'.format(OutPrefix)),'w')
        for line in sorted(sraWarningList):
            fout.write("WARNING: {}\n".format(line))
        fout.close()

    ## Write SRA runInfo to file
    if len(sraRunInfoDict) != 0:
        fout = open(os.path.join(OutDir,'{}.sra_runinfo.txt'.format(OutPrefix)),'w')
        fout.write('\t'.join(sorted(sraRunInfoHeader)) + '\n')
        for sra_id in sorted(sraRunInfoDict.keys()):
            rowList = []
            for col in sraRunInfoHeader:
                if col in sorted(sraRunInfoDict[sra_id]):
                    rowList.append(sraRunInfoDict[sra_id][col])
                else:
                    rowList.append('NA')
            fout.write('\t'.join(rowList) + '\n')
        fout.close()


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN,args.OUT_PREFIX,args.IGNORE_SRA_ERRORS,args.SKIP_SRA)


if __name__ == '__main__':
    sys.exit(main())
