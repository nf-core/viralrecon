#!/usr/bin/env python

# based on https://www.biostars.org/p/77526/

import os
import sys
import argparse
import requests
import csv

# def parse_args(args=None):
#     Description = 'Get metainformation of SRA samples, specifically the platform and the library layout (paired/single)'
#     Epilog = """Example usage: python get_SRA_metainfo.py SRA_ID <FILE_OUT>"""
#
#     parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
#     parser.add_argument('SRA_ID', help="Sample ID to.")
#     # parser.add_argument('FILE_OUT', help="Output samplesheet file.")
#
#     return parser.parse_args(args)

#  conda install -c bioconda parallel-fastq-dump
# conda install -c bioconda sra-tools
SRA_id = 'SRR11092056'
# SRA_id = 'SRR11177792'
# parallel-fastq-dump --sra-id SRR11092056 --split-files --threads 4

url = "https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=" + SRA_id

with requests.Session() as request:
    r = requests.post(url)
    csv = csv.DictReader(r.content.decode('utf-8').splitlines(), delimiter=',')
    for row in csv:
        print row
        print row['Platform']
        print row['download_path']
        print row['LibraryLayout']

# if __name__ == '__main__':
#     sys.exit(main())