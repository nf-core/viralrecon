############This code parsers .vcf and .snpsift files into long table ##############

# Import libraries
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


def parser_args(args=None):
    Description = 'Create long/wide tables fo ivar/bcftools'
    Epilog = """Example usage: python parser_ivar_bcftools.py <sample> <snpsift> <pangolin> <software> """
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('-sample', help="Input sample table files path.")
    parser.add_argument('-snpsift', help="Input snpsift txt files path.")
    parser.add_argument('-pangolin', help="Input pangolin csv files path.")
    parser.add_argument('-software', help="Input bcftools of ivar.")
   
    return parser.parse_args(args)

 
def files(sample,snpsift,pangolin):
    # function to obtain all required files
    path = sample
    sample_dir = os.listdir(sample)
    snpsift_dir = os.listdir(snpsift)
    pangolin_dir = os.listdir(pangolin)

    norm_tables = sorted([i for i in sample_dir if i.endswith("_norm.table")])
    norm_snpsift = sorted([i for i in snpsift_dir if i.endswith("_norm.snpsift.txt")])
    csv_pango  = sorted([i for i in pangolin_dir if i.endswith(".csv")])
    return path,norm_tables,norm_snpsift,csv_pango
    

def create_long_ivar(sample_dir,in_table,snp_table,pango_table):
    # creating long table for software ivar 
    software='ivar'
     ### format of sample table
    sample_table=pd.read_table(str(sample_dir) + "/"  + in_table, header='infer') 
    sample_table = sample_table.dropna(how = 'all', axis =1)
    sample_table.rename(columns={sample_table.columns[0]: "CHROM",sample_table.columns[1]: "POS",sample_table.columns[2]: "REF",sample_table.columns[3]: "ALT",sample_table.columns[4]: "FILTER",sample_table.columns[5]: "DP",sample_table.columns[6]: "REF_DP", sample_table.columns[7]: "ALT_DP"}, inplace=True)
    sample_table[["ALT_DP", "DP"]] = sample_table[["ALT_DP", "DP"]].apply(pd.to_numeric)
    sample_table['AF']=sample_table['ALT_DP']/sample_table['DP']
    sample_table['AF'] = sample_table['AF'].round(2)
   
    #snpsift function formats snpsift table
    snpsift_table_copy = f_snpsift(sample_dir,snp_table)

    #format of lineages table long one sample 
    tl_onesample= f_pangolin(sample_table,software,pango_table,sample_dir)

    return tl_onesample,snpsift_table_copy


def create_long_bcftools(sample_dir,in_table,snp_table,pango_table):
    software = 'bcftools'
     ### format of sample table
    sample_table=pd.read_table(str(sample_dir) + "/"  + in_table, header='infer') 
    sample_table = sample_table.dropna(how = 'all', axis =1)
    sample_table.rename(columns={sample_table.columns[5]: "DP",sample_table.columns[6]: "AD"}, inplace=True)
    new_column = sample_table
    new_column[['REF_DP','ALT_DP']] = sample_table['AD'].str.split(',', expand=True)
    sample_table = pd.merge(sample_table,new_column,how = 'left')
    sample_table[["ALT_DP", "DP"]] = sample_table[["ALT_DP", "DP"]].apply(pd.to_numeric)
    sample_table['AF']=sample_table['ALT_DP']/sample_table['DP']
    sample_table['AF'] = sample_table['AF'].round(2)
    sample_table = sample_table.loc[:, ~sample_table.columns.str.contains('AD')]
    
    #snpsift function formats snpsift table
    snpsift_table_copy = f_snpsift(sample_dir,snp_table)

    #format of lineages table long one sample 
    tl_onesample= f_pangolin(sample_table,software,pango_table,sample_dir)

    return tl_onesample,snpsift_table_copy


def f_snpsift(sample_dir,snp_table):

     ### format of snpsift table
    snpsift_table = pd.read_csv(str(sample_dir) + '/snpeff/' + snp_table, sep="\t", header = "infer")
    snpsift_table = snpsift_table.loc[:, ~snpsift_table.columns.str.contains('^Unnamed')]
    snpsift_table.rename(columns={snpsift_table.columns[0]: "CHROM",snpsift_table.columns[1]: "POS",snpsift_table.columns[2]: "REF",snpsift_table.columns[3]: "ALT",snpsift_table.columns[4]: "GENE",snpsift_table.columns[5]: "EFFECT",snpsift_table.columns[6]: "HGVS_C", snpsift_table.columns[7]: "HGVS_P"}, inplace=True)
    
    #removing , and . from columns
    snpsift_table_copy = snpsift_table.copy()
    for i in  range(len(snpsift_table_copy)):
        for j in range(3,8):
         snpsift_table_copy.iloc[i,j]= str(snpsift_table.iloc[i,j]).split(",")[0]
    return snpsift_table_copy

def f_pangolin(sample_table,software,pango_table,sample_dir):
    #format of lineages table long one sample 
    tl_onesample = pd.DataFrame(data =sample_table)

    if software == 'ivar':
        pangolin_table = pd.read_csv(str(sample_dir) +  '/pangolin/' + 'lineage_report.csv', sep=",", header = "infer")
        pangolin_table['taxon'] = pangolin_table['taxon'].astype(str).str[10:16]
        pangolin_table = pangolin_table[pangolin_table['taxon'] == pango_table[0:6] ]
    elif software=='bcftools':
        pangolin_table = pd.read_csv(str(sample_dir) +  '/pangolin/' + pango_table, sep=",", header = "infer")
        
    lineages = pangolin_table.loc[:,['taxon','lineage']]
    tl_onesample["sample"] = lineages.iloc[0,0]
    tl_onesample["software"] = software
    tl_onesample["lineage"] = lineages.iloc[0,1]
    
    return tl_onesample

def mergetables(sample_intable,snp_intable,sample_dir,norm_tables):
    #creating one table from snpsift and .vcf 
    left = snp_intable 
    right = sample_intable
    merged_table_long = pd.merge(left,right,how = 'outer')
    
    merged_table_long.to_csv(sample_dir + '/SampleTables/'+ norm_tables + '.csv', header='infer', index=None, sep=' ', mode='a')
    return(merged_table_long)

    
def concatenatetable(path_in_concat):
    #Creating final long table
    os.chdir(path_in_concat + '/SampleTables/')
    extension = 'csv'
    all_filenames = [i for i in glob.glob('*.{}'.format(extension))]

    #combine all files in the list
    list_of_dataframes = []
    for i in all_filenames:
        list_of_dataframes.append(pd.read_csv(i,sep=" "))
    merged_df = pd.concat(list_of_dataframes)
    merged_df.to_csv("final_long_table.csv", index=False, encoding='utf-8-sig')

    return merged_df



def main(args=None):

    
    args = parser_args(args)

    sample_dir,norm_tables,norm_snpsift,csv_pango = files(args.sample,args.snpsift,args.pangolin )
    
    for i in range(0,len(norm_tables)):

        if args.software =='ivar':
            tl_onesample_out,snpsift_table_final = create_long_ivar(sample_dir,norm_tables[i],norm_snpsift[i], csv_pango[i])
            merged = mergetables(tl_onesample_out,snpsift_table_final,sample_dir,norm_tables[i])

        elif args.software =='bcftools':
            tl_onesample_out,snpsift_table_final = create_long_bcftools(sample_dir,norm_tables[i],norm_snpsift[i], csv_pango[i])
            merged = mergetables(tl_onesample_out,snpsift_table_final,sample_dir,norm_tables[i])
   
    merged_df= concatenatetable(args.sample)
    

if __name__ == '__main__':
    sys.exit(main())