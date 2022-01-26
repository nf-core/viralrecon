#!/usr/bin/env python

from matplotlib import table
import pandas as pd
import sys
import numpy as np
import re
import argparse
import glob, os
import io

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

#Current working directory
cwd = os.getcwd()
os.chdir(cwd)

#create SampleTables folder in the folder where the script is running
if not os.path.isdir("SampleTables"):
    sample_table_path = os.path.join(cwd, "SampleTables")
    os.mkdir(sample_table_path)

def parser_args(args=None):
    Description = 'Create long/wide tables fo ivar/bcftools'
    Epilog = """Example usage: python parser_ivar_bcftools.py <sample> <snpsift> <pangolin> <software> """
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('-sample', help="Input sample table files path." )
    parser.add_argument('-snpsift', help="Input snpsift txt files path.",required=True)
    parser.add_argument('-pangolin', help="Input pangolin csv files path.",required=True)
    parser.add_argument('-software', help="Input bcftools of ivar.",required=True)

    return parser.parse_args(args)


def create_long(in_table,snp_table,pango_table,counter,software):
    if not in_table:
        intable_pathname = cwd
    else:
        intable_pathname = os.path.join(cwd, in_table)

    os.chdir(intable_pathname)
    table_list = []
    for file in glob.glob("*_norm.table"):
        table_list.append(file)
    table_list.sort()

    snp_pathname = os.path.join(cwd, snp_table)
    os.chdir(snp_pathname)
    snp_list = []
    for file in glob.glob("*_norm.snpsift.txt"):
        snp_list.append(file)
    snp_list.sort()

    pango_pathname = os.path.join(cwd, pango_table)
    os.chdir(pango_pathname)
    pangolin_list = []
    for file in glob.glob("*.csv"):
        pangolin_list.append(file)
    pangolin_list.sort()

     ### format of sample table
    sample_table=pd.read_table(str(intable_pathname) + "/"  + table_list[counter], header='infer')
    sample_table = sample_table.dropna(how = 'all', axis =1)

    if software=='bcftools':
        sample_table.rename(columns={sample_table.columns[5]: "DP",sample_table.columns[6]: "AD"}, inplace=True)

        new_column = sample_table
        new_column[['REF_DP','ALT_DP']] = sample_table['AD'].str.split(',', expand=True)
        sample_table = pd.merge(sample_table,new_column,how = 'left')
        sample_table[["ALT_DP", "DP"]] = sample_table[["ALT_DP", "DP"]].apply(pd.to_numeric)
        sample_table['AF']=sample_table['ALT_DP']/sample_table['DP']
        sample_table['AF'] = sample_table['AF'].round(2)
        sample_table = sample_table.loc[:, ~sample_table.columns.str.contains('AD')]

    elif software=='ivar':
        sample_table.rename(columns={sample_table.columns[0]: "CHROM",sample_table.columns[1]: "POS",sample_table.columns[2]: "REF",sample_table.columns[3]: "ALT",sample_table.columns[4]: "FILTER",sample_table.columns[5]: "DP",sample_table.columns[6]: "REF_DP", sample_table.columns[7]: "ALT_DP"}, inplace=True)
        sample_table[["ALT_DP", "DP"]] = sample_table[["ALT_DP", "DP"]].apply(pd.to_numeric)
        sample_table['AF']=sample_table['ALT_DP']/sample_table['DP']
        sample_table['AF'] = sample_table['AF'].round(2)

    ### format of snpsift table
    snpsift_table = pd.read_csv(str(snp_pathname) + '/' + snp_list[counter], sep="\t", header = "infer")
    snpsift_table = snpsift_table.loc[:, ~snpsift_table.columns.str.contains('^Unnamed')]
    colnames_snpsift = list(snpsift_table.columns)
    colnames_snpsift = [i.replace('ANN[*].', '') for i in colnames_snpsift]
    for i in range(len(colnames_snpsift)):
       snpsift_table.rename(columns = {snpsift_table.columns[i]:colnames_snpsift[i]}, inplace = True)
    snpsift_table =  snpsift_table.loc[:, ['CHROM','POS','REF','ALT','GENE','EFFECT','HGVS_C','HGVS_P']]

    snpsift_table_copy = snpsift_table.copy()

    for i in  range(len(snpsift_table_copy)):
        for j in range(3,8):
         snpsift_table_copy.iloc[i,j]= str(snpsift_table.iloc[i,j]).split(",")[0]

    #format of lineages
    pangolin_table = pd.read_csv(str(pango_pathname) +  '/' + pangolin_list[counter], sep=",", header = "infer")
    lineages = pangolin_table.loc[:,['taxon','lineage']]

    #table long one sample
    tl_onesample = pd.DataFrame(data =sample_table)
    if software=='bcftools':
        tl_onesample["Sample"] = lineages.iloc[0,0]
    elif software=='ivar':
        tl_onesample["Sample"] = lineages.iloc[0,0].split('_')[1].split('.')[0]
    tl_onesample["software"] = software
    tl_onesample["Lineage"] = lineages.iloc[0,1]

    return tl_onesample,snpsift_table_copy,counter

def mergetables(sample_intable,snp_intable,path,counter):

    if not path:
        merge_pathname = cwd
    else:
        merge_pathname = os.path.join(cwd, path)

    os.chdir(merge_pathname)
    table_list = []
    for file in glob.glob("*.table"):
        table_list.append(file[0:6])
    table_list.sort()

    left = snp_intable
    right = sample_intable
    merged_table_long = pd.merge(left,right,how = 'outer')

    merged_table_long.to_csv(str(merge_pathname) + '/SampleTables/'+ table_list[counter] + '.csv', header='infer', index=None, sep=' ', mode='a')
    return(merged_table_long)

def loop(path_in):
    if not path_in:
        in_pathname = cwd
    else:
        in_pathname = os.path.join(cwd, path_in)

    os.chdir(in_pathname)
    table_list = []
    for file in glob.glob("*.table"):
        table_list.append(file)
    table_list.sort()

    length = len(table_list)

    return length


def concatenatetable(path_in_concat):

    if not path_in_concat:
        concat_pathname = os.path.join(cwd, 'SampleTables')
    else:
        concat_pathname = os.path.join(cwd,path_in_concat,'SampleTables')

    os.chdir(concat_pathname)

    extension = 'csv'
    all_filenames = [i for i in glob.glob('*.{}'.format(extension))]

    #combine all files in the list
    #combined_csv = pd.concat([pd.read_csv(f) for f in all_filenames ])
    list_of_dataframes = []
    for i in all_filenames:
        list_of_dataframes.append(pd.read_csv(i,sep=" "))

    merged_df = pd.concat(list_of_dataframes)
    merged_df.to_csv("final_long_table.csv", index=False, encoding='utf-8-sig')
    return merged_df

def main(args=None):
    args = parser_args(args)
    length = loop(args.sample)

    for i in range(0,length):
        tl_onesample_out,snpsift_table_final,counter_it = create_long(args.sample,args.snpsift, args.pangolin,i,args.software)
        merged = mergetables(tl_onesample_out,snpsift_table_final,args.sample,counter_it)

    merged_df= concatenatetable(args.sample)


if __name__ == '__main__':
    sys.exit(main())
