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
    parser.add_argument('-ga', '--get_aspera_links', dest="GET_ASPERA_LINK", help="Output Aspera links instead of FTP.",action='store_true')
    parser.add_argument('-is', '--ignore_sra_errors', dest="IGNORE_SRA_ERRORS", help="Ignore SRA validation errors as opposed to failing.",action='store_true')
    parser.add_argument('-ss', '--skip_sra', dest="SKIP_SRA", help="Skip steps involving SRA identifiers.",action='store_true')

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


def fetch_url(url,encoding='utf-8'):
    try:
        r = requests.get(url)
    except requests.exceptions.RequestException as e:
        raise SystemExit(e)
    if r.status_code != 200:
        print("ERROR: Connection failed\nError code '{}'".format(r.status_code))
        sys.exit(1)

    return r.content.decode(encoding).splitlines()


## Get SRA ids for GSE, GSM ids.
    ## Can resolve SRA ids via ENA url
    ## Can resolve GSM (not GSE) ids via SRA url
    ## Cant resolve GEO ids via ENA url
    ## Cant resolve PRJ ids via ENA url
def id_to_sra(db_id):
    db_ids = []
    if db_id[:3] in ['ERR', 'GSM', 'PRJ', 'SAM', 'SRR']:
        db_ids = [db_id]
    elif db_id[:3] in ['GSE']:
        url = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={}&targ=gsm&view=data&form=text'.format(db_id)
        db_ids = [x.split('=')[1].strip() for x in fetch_url(url) if x.find('GSM') != -1]

    sra_ids = []
    for id in db_ids:
        delimiter = ','
        exp_col_name = 'Experiment'
        url = 'https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=' + id
        if id[:3] in ['ERR']:
            delimiter = '\t'
            exp_col_name = 'experiment_accession'
            url = 'http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession={}&result=read_run'.format(id)
        csv_dict = csv.DictReader(fetch_url(url), delimiter=delimiter)
        for row in csv_dict:
            sra_ids.append(row[exp_col_name])

    return set(sra_ids)

## SRA supported ids = ['PRJNA63463', 'SAMN00765663', 'SRA023522', 'SRP003255', 'SRR390278', 'SRS282569', 'SRX111814']
## ENA supported ids = ['ERA2421642', 'ERP120836', 'ERR674736', 'ERS4399631', 'ERX629702', 'PRJEB7743', 'SAMEA3121481']
## GEO supported ids = ['GSE18729', 'GSM465244']
## DBJ unsupported   = ['DRA000001', 'DRP000001', 'DRR000001', 'DRS000001', 'DRX000001', 'PRJDA38027', 'SAMD00016353']
def get_runinfo(db_id):
    db_ids = []
    if db_id[:3] in ['SRA', 'SRP', 'SRS', 'SRX', 'ERA', 'ERP', 'ERS', 'ERX']:
        db_ids = [db_id]
    elif db_id[:4] in ['PRJE', 'SAME']:
        db_ids = [db_id]
    elif db_id[:3] in ['SRR', 'ERR', 'GSE', 'GSM']:
        db_ids = id_to_sra(db_id)
    elif db_id[:4] in ['PRJN', 'SAMN']:
        db_ids = id_to_sra(db_id)

    runInfoDict = {}
    for id in db_ids:
        url = 'http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession={}&result=read_run'.format(id)
        csv_dict = csv.DictReader(fetch_url(url), delimiter='\t')
        for row in csv_dict:
            sampleid = row['experiment_accession']
            runid = row['run_accession']
            if not sampleid in runInfoDict:
                runInfoDict[sampleid] = {}
            runInfoDict[sampleid][runid] = row

    return runInfoDict


def check_samplesheet(FileIn,OutPrefix,getAsperaLinks=False,ignoreSRAErrors=False,skipSRA=False):
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
            sampleInfoList = []                                             ## [sample, single_end, is_sra, is_ftp, fastq_1, fastq_2, md5_1, md5_2]
            fastq_1,fastq_2 = fastQFiles
            if sample and fastq_1 and fastq_2:                              ## Paired-end short reads
                sampleInfoList = [[sample, '0', '0', '0', fastq_1, fastq_2, '', '']]

            elif sample and fastq_1 and not fastq_2:                        ## Single-end short reads
                sampleInfoList = [[sample, '1', '0', '0', fastq_1, fastq_2, '', '']]

            elif sample and not fastq_1 and not fastq_2:                    ## SRA accession
                if not skipSRA:
                    if sample[:2] in ['SR', 'SA', 'PR', 'ER', 'GS']:
                        is_sra = True
                        runInfoDict = get_runinfo(sample)
                        if len(runInfoDict) != 0:
                            for exp in runInfoDict.keys():
                                for run in runInfoDict[exp].keys():
                                    platform = runInfoDict[exp][run]['instrument_platform']
                                    library = runInfoDict[exp][run]['library_layout']
                                    fastq_files = runInfoDict[exp][run]['fastq_ftp']
                                    fastq_md5 = runInfoDict[exp][run]['fastq_md5']
                                    if getAsperaLinks:
                                        fastq_files = runInfoDict[exp][run]['fastq_aspera']
                                    if platform.lower() in ['illumina']:
                                        if library.lower() == 'single':
                                            if fastq_files != '':
                                                sampleInfoList.append([exp, '1', '1', '1', fastq_files , '', fastq_md5, ''])
                                            else:
                                                sampleInfoList.append([run, '1', '1', '0', '', '', '', ''])
                                        elif library.lower() == 'paired':
                                            if fastq_files != '':
                                                sampleInfoList.append([exp, '0', '1', '1'] + fastq_files.split(';') + fastq_md5.split(';'))
                                            else:
                                                sampleInfoList.append([run, '0', '1', '0', '', '', '', ''])
                                        else:
                                            errorStr = "Library layout '{}' != 'SINGLE' or 'PAIRED' for database id:'{}'! User provided id:'{}'.".format(layout,run,sample)
                                            if not ignoreSRAErrors:
                                                print_error(errorStr,line)
                                            sraWarningList.append(errorStr)
                                    else:
                                        errorStr = "Only Illumina platform is currently supported. Database id:'{}' is from the '{}' platform! User provided id:'{}'.".format(run,platform,sample)
                                        if not ignoreSRAErrors:
                                            print_error(errorStr,line)
                                        sraWarningList.append(errorStr)

                                    sraRunInfoHeader = list(set(sraRunInfoHeader) | set(runInfoDict[exp][run].keys()))
                                    sraRunInfoDict[run] = runInfoDict[exp][run]
                        else:
                            errorStr = "No data available for database id:{}!".format(sample)
                            if not ignoreSRAErrors:
                                print_error(errorStr,line)
                            sraWarningList.append(errorStr)
                    else:
                        print_error("Please provide a valid SRA id starting with 'SR', 'SA', 'PR', 'ER', 'GS'!",line)
            else:
                print_error("Invalid combination of columns provided!",line)

            for sampleList in sampleInfoList:
                sample_id = sampleList[0]; infoList = sampleList[1:]
                if sample_id not in sampleRunDict:
                    sampleRunDict[sample_id] = [infoList]
                else:
                    if infoList in sampleRunDict[sample_id]:
                        if is_sra:
                            errorStr = "Duplicate database id:'{}'! User provided id:'{}'.".format(sample_id,sample)
                            sraWarningList.append(errorStr)
                        else:
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
    fout.write(','.join(['sample_id', 'single_end', 'is_sra', 'is_ftp', 'fastq_1', 'fastq_2', 'md5_1', 'md5_2']) + '\n')
    for sample in sorted(sampleRunDict.keys()):

        ## Check that multiple runs of the same sample are of the same datatype
        if not all(x[:2] == sampleRunDict[sample][0][:2] for x in sampleRunDict[sample]):
            print_error("Multiple runs of a sample must be of the same datatype","Sample: {}".format(sample))

        for idx,val in enumerate(sampleRunDict[sample]):
            fout.write(','.join(["{}_T{}".format(sample,idx+1)] + val) + '\n')
    fout.close()

    ## Write SRA warnings to file
    if len(sraWarningList) != 0:
        fout = open(os.path.join(OutDir,'{}.sra_warnings.txt'.format(OutPrefix)),'w')
        for line in sorted(set(sraWarningList)):
            fout.write("WARNING: {}\n".format(line))
        fout.close()

    ## Write SRA runInfo to file
    if len(sraRunInfoDict) != 0:
        fout = open(os.path.join(OutDir,'{}.sra_runinfo.txt'.format(OutPrefix)),'w')
        fout.write('"{}"'.format('"\t "'.join(map(str, sorted(sraRunInfoHeader)))) + '\n')
        for db_id in sorted(sraRunInfoDict.keys()):
            rowList = []
            for col in sorted(sraRunInfoHeader):
                if col in sraRunInfoDict[db_id]:
                    rowList.append(sraRunInfoDict[db_id][col])
                else:
                    rowList.append('NA')
            fout.write('"{}"'.format('"\t "'.join(map(str, rowList))) + '\n')
        fout.close()


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN,args.OUT_PREFIX,args.GET_ASPERA_LINK,args.IGNORE_SRA_ERRORS,args.SKIP_SRA)


if __name__ == '__main__':
    sys.exit(main())
