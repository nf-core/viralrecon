#!/usr/bin/env python

# based on https://www.biostars.org/p/77526/

#  conda install -c bioconda parallel-fastq-dump
# conda install -c bioconda sra-tools
# SRA_id = 'SRR11092056'
# SRA_id = 'SRR11140744'
# SRA_id = 'SRR11177792'
# parallel-fastq-dump --sra-id SRR11092056 --split-files --threads 4

import os
import sys
import argparse
import requests
import csv

def parse_args(args=None):
    Description = 'Gets metainformation of SRA samples'
    # Epilog = """Example usage: python get_SRA_metainfo.py SRA_ID <FILE_OUT>"""
    Epilog = """Example usage: python get_SRA_metainfo.py SRA_ID"""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('SRA_ID', help="Sample ID to obtain metainfo.")
    parser.add_argument('--is_single', required=False, action='store_true',
                        default=False, help='Returns true if SRA is single otherwise false')

    # parser.add_argument('FILE_OUT', help="Output samplesheet file.")
    return parser.parse_args(args)

def print_error(error):
    print("ERROR: {}\n".format(error))
    sys.exit(1)

def get_sra_metainfo(sra_id, is_single_bool):
    url = "https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=" + sra_id

    if sra_id[:3] != 'SRR':
        print_error("Please provide a valid SRA run accession starting with 'SRR'!" )

    try:
        r = requests.get(url)
    except requests.exceptions.RequestException as e:
        raise SystemExit(e)

    if r.status_code != 200:
        print_error("ERROR: Connection to sra failed\nError code '{}'".format(r.status_code))

    csv_content = csv.DictReader(r.content.decode('utf-8').splitlines(), delimiter=',')
    empty_file = True
    sampleRunDict = {}

    for row in csv_content:
        empty_file = False

        single_end = '0'; is_sra = '1'; fastq_1=""; fastq_2=""

        if row['Platform'] != "ILLUMINA":
            print_error("ERROR: Pipeline currently only suports ILLUMINA reads\nThe requested SRA corresponds to '{}' reads".format(row['Platform']))

        if row['LibraryLayout'] == "SINGLE":
            single_end = '1'; fastq_1 = sra_id + "_1.fastq.gz"
            if is_single_bool: print ("true"); exit(0)

        elif row['LibraryLayout'] == "PAIRED":
            single_end = '0'; fastq_1 = sra_id + "_1.fastq.gz"; fastq_2 = sra_id + "_2.fastq.gz"
            if is_single_bool: print("false"); exit(0)
        else:
            row['LibraryLayout']
            print_error("ERROR: Library layout unknown '{}'".format(row['LibraryLayout']))

        # print (row['Platform'])
        # print (row['download_path']) #DOWNLOAD PATH
        # print (row['LibraryLayout'])

        sampleInfoList = [single_end, is_sra, fastq_1, fastq_2]
        sampleRunDict[sra_id] = []
        # if sample not in sampleRunDict:
        #     sampleRunDict[sample] = []
        # else:
        #     if sampleInfoList in sampleRunDict[sample]:
        #         print_error("Samplesheet contains duplicate rows!", line)
        sampleRunDict[sra_id].append(sampleInfoList)

    if empty_file:
        print_error("No data available for the provided SRA_id: {}".format(sra_id))

    # Write to file
    file_name =  sra_id + "_metainfo.csv"
    fout = open(file_name, 'w')
    fout.write(','.join(['sample_id', 'single_end', 'is_sra', 'fastq_1', 'fastq_2']) + '\n')
    for sample in sorted(sampleRunDict.keys()):

        ## Check that multiple runs of the same sample are of the same datatype
        if not all(x[:2] == sampleRunDict[sample][0][:2] for x in sampleRunDict[sample]):
            print_error("Multiple runs of a sample must be of the same datatype", "Sample: {}".format(sample))

        for idx, val in enumerate(sampleRunDict[sample]):
            sample_id = "{}_T{}".format(sample, idx + 1)
            fout.write(','.join([sample_id] + val) + ',\n')
    fout.close()

def is_single(sra_id):

    url = "https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=" + sra_id

    if sra_id[:3] != 'SRR':
        print_error("Please provide a valid SRA run accession starting with 'SRR'!")

    try:
        r = requests.get(url)
    except requests.exceptions.RequestException as e:
        raise SystemExit(e)

    if r.status_code != 200:
        print_error("ERROR: Connection to sra failed\nError code '{}'".format(r.status_code))

    csv_content = csv.DictReader(r.content.decode('utf-8').splitlines(), delimiter=',')
    empty_file = True
    sampleRunDict = {}

    for row in csv_content:
        empty_file = False

        single_end = '0';
        is_sra = '1';
        fastq_1 = "";
        fastq_2 = ""

        if row['Platform'] != "ILLUMINA":
            print_error(
                "ERROR: Pipeline currently only suports ILLUMINA reads\nThe requested SRA corresponds to '{}' reads".format(
                    row['Platform']))

        if row['LibraryLayout'] == "SINGLE":
            print ("SINGLE")
            break
        elif row['LibraryLayout'] == "PAIRED":
            print ("PAIRED")
            break
        else:
            print (row['LibraryLayout'])

    if empty_file:
        print_error("No data available for the provided SRA_id: {}".format(sra_id))

def main(args=None):
    args = parse_args(args)

    get_sra_metainfo(args.SRA_ID, args.is_single)

if __name__ == '__main__':
    sys.exit(main())
