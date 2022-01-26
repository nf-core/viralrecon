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


def create_long(snp_file,snpsift_file,pangolin_file,software):

    ### format of sample table
    snp_table=pd.read_table(snp_file, header='infer')
    snp_table = snp_table.dropna(how = 'all', axis =1)

    if software=='bcftools':
        snp_table.rename(columns={snp_table.columns[5]: "DP",snp_table.columns[6]: "AD"}, inplace=True)
        new_column = snp_table
        new_column[['REF_DP','ALT_DP']] = snp_table['AD'].str.split(',', expand=True)
        snp_table = pd.merge(snp_table,new_column,how = 'left')
        snp_table[["ALT_DP", "DP"]] = snp_table[["ALT_DP", "DP"]].apply(pd.to_numeric)
        snp_table['AF']=snp_table['ALT_DP']/snp_table['DP']
        snp_table['AF'] = snp_table['AF'].round(2)
        snp_table = snp_table.loc[:, ~snp_table.columns.str.contains('AD')]

    elif software=='ivar':
        snp_table.rename(columns={snp_table.columns[0]: "CHROM",snp_table.columns[1]: "POS",snp_table.columns[2]: "REF",snp_table.columns[3]: "ALT",snp_table.columns[4]: "FILTER",snp_table.columns[5]: "DP",snp_table.columns[6]: "REF_DP", snp_table.columns[7]: "ALT_DP"}, inplace=True)
        snp_table[["ALT_DP", "DP"]] = snp_table[["ALT_DP", "DP"]].apply(pd.to_numeric)
        snp_table['AF']=snp_table['ALT_DP']/snp_table['DP']
        snp_table['AF'] = snp_table['AF'].round(2)

    ### format of snpsift table
    snpsift_table = pd.read_csv(snpsift_file, sep="\t", header = "infer")
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
    pangolin_table = pd.read_csv(pangolin_file, sep=",", header = "infer")
    lineages = pangolin_table.loc[:,['taxon','lineage']]

    #table long one sample
    tl_onesample = pd.DataFrame(data =snp_table)
    if software=='bcftools':
        tl_onesample["Sample"] = lineages.iloc[0,0]
    elif software=='ivar':
        tl_onesample["Sample"] = lineages.iloc[0,0].split('_')[1].split('.')[0]
    tl_onesample["software"] = software
    tl_onesample["Lineage"] = lineages.iloc[0,1]

    return tl_onesample,snpsift_table_copy

def mergetables(sample_intable,snp_intable,path,counter):

#     if not path:
#         merge_pathname = cwd
#     else:
#         merge_pathname = os.path.join(cwd, path)

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

def same_len(list):
    return len(set(list)) == 1

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

    # List vcf table files
    table_list = []
    for file in glob.glob(args.sample + "/*_norm.table"):
        table_list.append(file)
    table_list.sort()

    #List snpsift files
    snpsift_list = []
    for file in glob.glob(args.snpsift + "/*_norm.snpsift.txt"):
        snpsift_list.append(file)
    snpsift_list.sort()

    # List pangolin files
    pangolin_list = []
    for file in glob.glob(args.pangolin + "/*.csv"):
        pangolin_list.append(file)
    pangolin_list.sort()

    if not same_len([len(table_list),len(snpsift_list),len(pangolin_list)]):
        print("not same number of files for variants, snpsift and pangolin results ")
        exit()

    sample_names = [os.path.basename(filename).split("_norm.table")[0] for filename in table_list]
    for table,snpsift,pangolin in zip(table_list,snpsift_list,pangolin_list):
        tl_onesample_out,snpsift_table_final = create_long(table,snpsift,pangolin,args.software)
        #merged = mergetables(tl_onesample_out,snpsift_table_final,args.sample,counter_it)

#      merged_df= concatenatetable(args.sample)

if __name__ == '__main__':
    sys.exit(main())
