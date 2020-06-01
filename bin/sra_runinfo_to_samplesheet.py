#!/usr/bin/env python

import os
import sys
import errno
import argparse

def parse_args(args=None):
    Description = "Create valid nf-core/viralrecon samplesheet file from output of 'fetch_sra_runinfo.py' script."
    Epilog = """Example usage: python sra_runinfo_to_samplesheet.py <FILE_IN> <FILE_OUT>"""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('FILE_IN', help="Input metadata file created from 'fetch_sra_runinfo.py' script.")
    parser.add_argument('FILE_OUT', help="Output file.")
    return parser.parse_args(args)


def make_dir(path):
    if not len(path) == 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


def sra_runinfo_to_samplesheet(FileIn,FileOut):

    sampleRunDict = {}
    fin = open(FileIn,'r')
    header = fin.readline().strip().split('\t')
    while True:
        line = fin.readline()
        if line:
            line_dict = dict(zip(header,line.strip().split('\t')))
            run_id = line_dict['run_accession']
            exp_id = line_dict['experiment_accession']
            library = line_dict['library_layout']
            fastq_files = line_dict['fastq_ftp']
            fastq_md5 = line_dict['fastq_md5']

            db_id = exp_id
            sample_info = []          ## [single_end, is_sra, is_ftp, fastq_1, fastq_2, md5_1, md5_2]
            if library == 'SINGLE':
                if fastq_files:
                    sample_info = ['1', '1', '1', fastq_files , '', fastq_md5, '']
                else:
                    db_id = run_id
                    sample_info = ['1', '1', '0', '', '', '', '']
            elif library == 'PAIRED':
                if fastq_files:
                    fq_files = fastq_files.split(';')[-2:]
                    if fq_files[0].find('_1.fastq.gz') != -1 and fq_files[1].find('_2.fastq.gz') != -1:
                        sample_info = ['0', '1', '1'] + fq_files + fastq_md5.split(';')[-2:]
                    else:
                        print("Invalid FastQ files found for database id:'{}'!.".format(run_id))
                else:
                    db_id = run_id
                    sample_info = ['0', '1', '0', '', '', '', '']

            if sample_info:
                if db_id not in sampleRunDict:
                    sampleRunDict[db_id] = [sample_info]
                else:
                    if sample_info in sampleRunDict[db_id]:
                        print("Input run info file contains duplicate rows!\nLine: '{}'".format(line))
                    else:
                        sampleRunDict[db_id].append(sample_info)
        else:
            break
            fin.close()

    ## Write samplesheet with appropriate columns
    if len(sampleRunDict) != 0:
        OutDir = os.path.dirname(FileOut)
        make_dir(OutDir)
        fout = open(FileOut,'w')
        fout.write(','.join(['sample_id', 'single_end', 'is_sra', 'is_ftp', 'fastq_1', 'fastq_2', 'md5_1', 'md5_2']) + '\n')
        for db_id in sorted(sampleRunDict.keys()):
            for idx,val in enumerate(sampleRunDict[db_id]):
                fout.write(','.join(["{}_T{}".format(db_id,idx+1)] + val) + '\n')
        fout.close()


def main(args=None):
    args = parse_args(args)
    sra_runinfo_to_samplesheet(args.FILE_IN,args.FILE_OUT)


if __name__ == '__main__':
    sys.exit(main())
