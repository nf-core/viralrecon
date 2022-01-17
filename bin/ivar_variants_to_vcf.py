#!/usr/bin/env python
import os
import sys
import re
import errno
import argparse
import numpy as np
from scipy.stats import fisher_exact


def parse_args(args=None):
    Description = "Convert iVar variants tsv file to vcf format."
    Epilog = """Example usage: python ivar_variants_to_vcf.py <FILE_IN> <FILE_OUT>"""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input tsv file.")
    parser.add_argument("FILE_OUT", help="Full path to output vcf file.")
    parser.add_argument(
        "-po",
        "--pass_only",
        dest="PASS_ONLY",
        help="Only output variants that PASS all filters.",
        action="store_true",
    )
    parser.add_argument(
        "-af",
        "--allele_freq_thresh",
        type=float,
        dest="ALLELE_FREQ_THRESH",
        default=0,
        help="Only output variants where allele frequency greater than this number (default: 0).",
    )
    parser.add_argument(
        "-nsb",
        "--not_strand_bias",
        dest="NOT_STRAND_BIAS",
        default=False,
        help="Does not take into account strand bias, use this option when not using amplicons for sequencing",
        action="store_true"
    )
    parser.add_argument(
        "-nmc",
        "--not_merge_codons",
        dest="NOT_MERGE_CODONS",
        help="Only output variants without taking into accout if the positions are consecutive and belong to the same codon.",
        action="store_true"
    )

    return parser.parse_args(args)

def checkConsecutive(mylist):
    '''
    Input dict_lines['POS'] position list.
    If len(list) = 3 and consecutive positions, returns "triple"
    If len(list) = 2 and consecutive positions, returns "double"
    If not consecutive returns False
    '''
    my_list = list(map(int, mylist))
    if sorted(my_list) == list(range(min(my_list), max(my_list)+1)):
        return len(my_list)
    else:
        if len(my_list) > 1:
            my_list.pop()
            if sorted(my_list) == list(range(min(my_list), max(my_list)+1)):
                return len(my_list)
        else: 
            return False
        return False

def is_same_codon(seq1,seq2):
    '''
    Returns position where seq1 != seq2 
    '''
    if seq1 =="NA":
        return False

    ind_diff = [i for i in range(len(seq1)) if seq1[i] != seq2[i]]
    if len(ind_diff) > 1:
            print("There has been an issue, more than one difference between the seqs.")
            return False
    else:        
        return ind_diff[0]
        
def renameVars(dict_lines,dummy):
    CHROM = dict_lines["CHROM"][0]
    POS = dict_lines["POS"][0]
    ID = dict_lines["ID"][0]
    if dummy =="double":
        REF = str(dict_lines["REF"][0]) + str(dict_lines["REF"][1]) 
        ALT = str(dict_lines["ALT"][0]) + str(dict_lines["ALT"][1]) 
    elif dummy =="triple":
        REF = str(dict_lines["REF"][0]) + str(dict_lines["REF"][1]) + str(dict_lines["REF"][2])
        ALT = str(dict_lines["ALT"][0]) + str(dict_lines["ALT"][1]) + str(dict_lines["ALT"][2])
    ## TODO Check how much differences we found among DPs in the three positions of a codon.
    REF_DP = dict_lines["REF_DP"][0]
    REF_RV = dict_lines["REF_RV"][0]
    ALT_DP = dict_lines["ALT_DP"][0]
    ALT_RV = dict_lines["ALT_RV"][0]
    QUAL = dict_lines["QUAL"][0]
    REF_CODON = REF
    ALT_CODON = ALT
    FILTER =dict_lines["FILTER"][0]  
    # INFO DP depends on the decision in the todo above. SB is left with the first one (en principio)
    INFO = dict_lines["INFO"][0]  
    FORMAT = dict_lines["FORMAT"][0]
    # sample depends on the decision in the todo above.
    SAMPLE = dict_lines["SAMPLE"][0]
    return CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,SAMPLE

def make_dir(path):
    if not len(path) == 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


def ivar_variants_to_vcf(FileIn, FileOut, passOnly=False, minAF=0,NotStrandBias=False,NotMergeCodon=False):
    filename = os.path.splitext(FileIn)[0]
    header = (
        "##fileformat=VCFv4.2\n"
        "##source=iVar\n"
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n'
        '##FILTER=<ID=PASS,Description="Result of p-value <= 0.05">\n'
        '##FILTER=<ID=FAIL,Description="Result of p-value > 0.05">\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '##FORMAT=<ID=REF_DP,Number=1,Type=Integer,Description="Depth of reference base">\n'
        '##FORMAT=<ID=REF_RV,Number=1,Type=Integer,Description="Depth of reference base on reverse reads">\n'
        '##FORMAT=<ID=REF_QUAL,Number=1,Type=Integer,Description="Mean quality of reference base">\n'
        '##FORMAT=<ID=ALT_DP,Number=1,Type=Integer,Description="Depth of alternate base">\n'
        '##FORMAT=<ID=ALT_RV,Number=1,Type=Integer,Description="Depth of alternate base on reverse reads">\n'
        '##FORMAT=<ID=ALT_QUAL,Number=1,Type=String,Description="Mean quality of alternate base">\n'
        '##FORMAT=<ID=ALT_FREQ,Number=1,Type=String,Description="Frequency of alternate base">\n'
    )
    header += (
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tINFO[SB]\tFORMAT\t" + filename + "\n"
    )
    
    varList = []
    varCountDict = {"SNP": 0, "INS": 0, "DEL": 0}
    dict_lines = {'CHROM':[],'POS':[],'ID':[],'REF':[],'ALT':[],'REF_DP':[],'REF_RV':[],'ALT_DP':[],'ALT_RV':[],'QUAL':[],'REF_CODON':[],'ALT_CODON':[],'FILTER': [],'INFO':[],'FORMAT':[],'SAMPLE':[]}
    writeLine=False
    write2Line = False
    OutDir = os.path.dirname(FileOut)
    make_dir(OutDir)
    fout = open(FileOut, "w")
    fout.write(header)

    with open(FileIn) as f:
        
        for line in f:

            if not re.match("REGION", line):
                
                line = re.split("\t", line)
                CHROM = line[0]
                POS = line[1]
                ID = "."
                REF = line[2]
                ALT = line[3]
                REF_DP = int(line[4])
                REF_RV = int(line[5])
                REF_FW = REF_DP - REF_RV
                ALT_RV = int(line[8])
                ALT_DP = int(line[7])
                ALT_FW = ALT_DP - ALT_RV
                table = np.array([[REF_FW, REF_RV], [ALT_FW, ALT_RV]])
                oddsr, p = fisher_exact(table, alternative='greater')

                var_type = "SNP"
                if ALT[0] == "+":
                    ALT = REF + ALT[1:]
                    var_type = "INS"
                elif ALT[0] == "-":
                    REF += ALT[1:]
                    ALT = line[2]
                    var_type = "DEL"

                QUAL = "."
                pass_test = line[13]
                REF_CODON = line[15]
                ALT_CODON = line[17]
              
                
                
                

                if NotStrandBias:

                    if pass_test =="TRUE":
                        FILTER = "PASS"
                    else:
                        FILTER = "FAIL"
                    INFO = "DP=" + line[11]
                else:

                    if p<0.05 and pass_test =="TRUE":
                        FILTER = "SB"
                    elif p>0.05 and pass_test =="TRUE":
                        FILTER = "PASS"
                    elif  p<=0.05 and pass_test == "FALSE":
                        FILTER = "SB,other"
                    else:
                        FILTER = "FAIL"
                    INFO = "DP=" + line[11]+":SB_pvalue="+str(round(p,5))

                
                FORMAT = "GT:REF_DP:REF_RV:REF_QUAL:ALT_DP:ALT_RV:ALT_QUAL:ALT_FREQ"
                SAMPLE = (
                    "1:"
                    + line[4]
                    + ":"
                    + line[5]
                    + ":"
                    + line[6]
                    + ":"
                    + line[7]
                    + ":"
                    + line[8]
                    + ":"
                    + line[9]
                    + ":"
                    + line[10]
                )
                param_list = [CHROM,POS,ID,REF,ALT,REF_DP,REF_RV,ALT_DP,ALT_RV,QUAL,REF_CODON,ALT_CODON,FILTER,INFO,FORMAT,SAMPLE]
                
                if NotMergeCodon:
                    writeLine = True
                    oline = (CHROM+ "\t"+ POS+ "\t"+ ID+ "\t"+ REF+ "\t"+ ALT+ "\t"+ QUAL+ "\t"+ FILTER+ "\t"+ INFO+ "\t"+ FORMAT+ "\t"+ SAMPLE+ "\n" )

                else:
                    #For the first two lines of the files, write without checking if the 3 positions are consecutive
                    if len(dict_lines["POS"]) ==0 or len(dict_lines["POS"]) ==1  :
                        for i,j in enumerate(dict_lines):
                            dict_lines.setdefault(j, []).append(param_list[i])
                        writeLine=False


                    elif  len(dict_lines["POS"]) == 2:

                        for i,j in enumerate(dict_lines):
                            dict_lines.setdefault(j, []).append(param_list[i])

                        #¿Consecutive?
                        if checkConsecutive(dict_lines["POS"]) == 2:
                            if is_same_codon(dict_lines["REF_CODON"][0],dict_lines["ALT_CODON"][0]) != 2:
                                writeLine = True
                                dummy = "double"
                                CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,SAMPLE = renameVars(dict_lines,dummy)
                                oline = (CHROM+ "\t"+ POS+ "\t"+ ID+ "\t"+ REF+ "\t"+ ALT+ "\t"+ QUAL+ "\t"+ FILTER+ "\t"+ INFO+"\t"+   FORMAT+ "\t"+ SAMPLE+ "\n" )
                                for i,j in enumerate(dict_lines):
                                    dict_lines[list(dict_lines.keys())[i]].pop(0)
                                    dict_lines[list(dict_lines.keys())[i]].pop(0)

                        elif checkConsecutive(dict_lines["POS"]) == 3:

                            if is_same_codon(dict_lines["REF_CODON"][0],dict_lines["ALT_CODON"][0]) == 0:
                                writeLine = True
                                dummy = "triple"
                                CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,SAMPLE = renameVars(dict_lines,dummy)
                                oline = (CHROM+ "\t"+ POS+ "\t"+ ID+ "\t"+ REF+ "\t"+ ALT+ "\t"+ QUAL+ "\t"+ FILTER+ "\t"+ INFO+ "\t"+  FORMAT+ "\t"+ SAMPLE+ "\n" )
                                for i,j in enumerate(dict_lines):
                                    dict_lines[list(dict_lines.keys())[i]].pop(0)
                                    dict_lines[list(dict_lines.keys())[i]].pop(0)
                                dict_lines = {'CHROM':[],'POS':[],'ID':[],'REF':[],'ALT':[],'REF_DP':[],'REF_RV':[],'ALT_DP':[],    'ALT_RV':[],'QUAL':[],'REF_CODON':[],'ALT_CODON':[],'FILTER':[],'INFO':[],'FORMAT':[],'SAMPLE':[]}
                            
                            elif is_same_codon(dict_lines["REF_CODON"][0],dict_lines["ALT_CODON"][0]) == 1:
                                writeLine = True
                                dummy = "double"
                                CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,SAMPLE = renameVars(dict_lines,dummy)
                                oline = (CHROM+ "\t"+ POS+ "\t"+ ID+ "\t"+ REF+ "\t"+ ALT+ "\t"+ QUAL+ "\t"+ FILTER+ "\t"+ INFO+ "\t"+  FORMAT+ "\t"+ SAMPLE+ "\n" ) 
                                for i,j in enumerate(dict_lines):
                                    dict_lines[list(dict_lines.keys())[i]].pop(0)
                                    dict_lines[list(dict_lines.keys())[i]].pop(0)
                            elif is_same_codon(dict_lines["REF_CODON"][0],dict_lines["ALT_CODON"][0]) == 2:
                                writeLine = True
                                oline =(dict_lines["CHROM"][0] + "\t"+ dict_lines["POS"][0]+ "\t"+ dict_lines["ID"][0]+ "\t"+ dict_lines    ["REF"][0]+ "\t"+ dict_lines["ALT"][0]+ "\t"+ dict_lines["QUAL"][0]+ "\t"+ dict_lines["FILTER"][0]+ "\t"    + dict_lines["INFO"][0]+ "\t"+ dict_lines["FORMAT"][0]+ "\t"+ dict_lines["SAMPLE"][0]+ "\n")

                                for i,j in enumerate(dict_lines):
                                    dict_lines[list(dict_lines.keys())[i]].pop(0)
                        elif checkConsecutive(dict_lines["POS"]) == False:
                            writeLine = True

                            oline =(dict_lines["CHROM"][0] + "\t"+ dict_lines["POS"][0]+ "\t"+ dict_lines["ID"][0]+ "\t"+ dict_lines    ["REF"][0]+ "\t"+ dict_lines["ALT"][0]+ "\t"+ dict_lines["QUAL"][0]+ "\t"+ dict_lines["FILTER"][0]+ "\t"+   dict_lines["INFO"][0]+ "\t"+ dict_lines["FORMAT"][0]+ "\t"+ dict_lines["SAMPLE"][0]+ "\n")

                            for i,j in enumerate(dict_lines):
                                dict_lines[list(dict_lines.keys())[i]].pop(0)
                    else:
                        print("Something went terribly wrong!!" + str(len(dict_lines["POS"])))        

                
                      
                if passOnly and FILTER != "PASS":
                    writeLine = False
                if float(line[10]) < minAF:
                    writeLine = False
                if (CHROM, POS, REF, ALT) in varList:
                    writeLine = False
                else:
                    varList.append((CHROM, POS, REF, ALT))
                if writeLine:
                    varCountDict[var_type] += 1
                    fout.write(oline)

            

    

    ## Print variant counts to pass to MultiQC
    varCountList = [(k, str(v)) for k, v in sorted(varCountDict.items())]
    print("\t".join(["sample"] + [x[0] for x in varCountList]))
    print("\t".join([filename] + [x[1] for x in varCountList]))


    if  len(dict_lines["POS"]) == 2:
        if checkConsecutive(dict_lines["POS"]) == 2:
            if is_same_codon(dict_lines["REF_CODON"][0],dict_lines["ALT_CODON"][0]) != 2:
                writeLine = True
                dummy = "double"
                CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,SAMPLE = renameVars(dict_lines,dummy)
                oline = (CHROM+ "\t"+ POS+ "\t"+ ID+ "\t"+ REF+ "\t"+ ALT+ "\t"+ QUAL+ "\t"+ FILTER+ "\t"+ INFO+"\t"+   FORMAT+ "\t"+ SAMPLE+ "\n" )
                fout.write(oline)
        else:
            oline =(dict_lines["CHROM"][0] + "\t"+ dict_lines["POS"][0]+ "\t"+ dict_lines["ID"][0]+ "\t"+ dict_lines    ["REF"][0]+ "\t"+ dict_lines["ALT"][0]+ "\t"+ dict_lines["QUAL"][0]+ "\t"+ dict_lines["FILTER"][0]+ "\t"+   dict_lines["INFO"][0]+ "\t"+ dict_lines["FORMAT"][0]+ "\t"+ dict_lines["SAMPLE"][0]+ "\n")
            oline1 =(dict_lines["CHROM"][1] + "\t"+ dict_lines["POS"][1]+ "\t"+ dict_lines["ID"][1]+ "\t"+ dict_lines    ["REF"][1]+ "\t"+ dict_lines["ALT"][1]+ "\t"+ dict_lines["QUAL"][1]+ "\t"+ dict_lines["FILTER"][1]+ "\t"+   dict_lines["INFO"][1]+ "\t"+ dict_lines["FORMAT"][1]+ "\t"+ dict_lines["SAMPLE"][1]+ "\n")
            fout.write(oline)
            fout.write(oline1)
    elif  len(dict_lines["POS"]) == 1:
        oline =(dict_lines["CHROM"][0] + "\t"+ dict_lines["POS"][0]+ "\t"+ dict_lines["ID"][0]+ "\t"+ dict_lines    ["REF"][0]+ "\t"+ dict_lines["ALT"][0]+ "\t"+ dict_lines["QUAL"][0]+ "\t"+ dict_lines["FILTER"][0]+ "\t"+   dict_lines["INFO"][0]+ "\t"+ dict_lines["FORMAT"][0]+ "\t"+ dict_lines["SAMPLE"][0]+ "\n")
        fout.write(oline)


def main(args=None):
    args = parse_args(args)
    ivar_variants_to_vcf(
        args.FILE_IN, args.FILE_OUT, args.PASS_ONLY, args.ALLELE_FREQ_THRESH, args.NOT_STRAND_BIAS, args.NOT_MERGE_CODONS
    )


if __name__ == "__main__":
    sys.exit(main())