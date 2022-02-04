#!/usr/bin/env python

import os
import sys
import re
import errno
import argparse
import numpy as np
from scipy.stats import fisher_exact


def parse_args(args=None):
    Description = "Convert iVar variants TSV file to VCF format."
    Epilog = """Example usage: python ivar_variants_to_vcf.py <file_in> <file_out>"""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("file_in", help="Input iVar TSV file.")
    parser.add_argument("file_out", help="Full path to output VCF file.")
    parser.add_argument(
        "-po",
        "--pass_only",
        help="Only output variants that PASS filters.",
        action="store_true",
    )
    parser.add_argument(
        "-af",
        "--allele_freq_threshold",
        type=float,
        default=0,
        help="Only output variants where allele frequency is greater than this number (default: 0).",
    )
    parser.add_argument(
        "-is",
        "--ignore_strand_bias",
        default=False,
        help="Does not take strand bias into account, use this option when not using amplicon sequencing.",
        action="store_true"
    )
    parser.add_argument(
        "-ic",
        "--ignore_merge_codons",
        help="Output variants without taking into account if consecutive positions belong to the same codon.",
        action="store_true"
    )

    return parser.parse_args(args)


def check_consecutive(mylist):
    '''
    Description:
        This function checks if a list of three or two  numbers are consecutive and returns how many items are consecutive.
    input:
        my_list - A list of integers
    return:
        Number of items consecutive in the list - [False, 1, 2]
    '''
    my_list = list(map(int, mylist))

    ## Check if the list contains consecutive numbers
    if sorted(my_list) == list(range(min(my_list), max(my_list)+1)):
        return len(my_list)
    else:
        ## If not, and the list is > 1, remove the last item and reevaluate.
        if len(my_list) > 1:
            my_list.pop()
            if sorted(my_list) == list(range(min(my_list), max(my_list)+1)):
                return len(my_list)
        else:
            return False
        return False


def codon_position(seq1,seq2):
    '''
    Description:
        Function to compare two codon nucleotide sequences (size 3) and retuns the position where it differs.
    Input:
        seq1 - list size 3 [A,T,C,G]
        seq2 - list size 3 [A,T,C,G]
    Returns:
        Returns position where seq1 != seq2
    '''
    if seq1 == "NA":
        return False

    ind_diff = [i for i in range(len(seq1)) if seq1[i] != seq2[i]]
    if len(ind_diff) > 1:
            print("There has been an issue, more than one difference between the seqs.")
            return False
    else:
        return ind_diff[0]


def rename_vars(dict_lines,num_collapse):
    '''
    Description:
        The function set the vars acordingly to the lines to collapse do to consecutive variants.
    Input:
        dict_lines - Dict with var lines.
        num_collapse - number of lines to collapse [2,3]
    Returns::
        Vars fixed.
    '''
    CHROM = dict_lines["CHROM"][0]
    POS = dict_lines["POS"][0]
    ID = dict_lines["ID"][0]
    # If two consecutive collapse 2 lines into one.
    if int(num_collapse) == 2:
        REF = str(dict_lines["REF"][0]) + str(dict_lines["REF"][1])
        ALT = str(dict_lines["ALT"][0]) + str(dict_lines["ALT"][1])
    # If three consecutive collapse 3 lines into one.
    elif int(num_collapse) == 3:
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
    # INFO DP depends on the decision in the todo above. SB is left with the first one.
    INFO = dict_lines["INFO"][0]
    FORMAT = dict_lines["FORMAT"][0]
    # sample depends on the decision in the todo above.
    SAMPLE = dict_lines["SAMPLE"][0]
    return CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,SAMPLE


def make_dir(path):
    '''
    Description:
        Create directory if it doesn't exist.
    Input:
        path - path where the directory will be created.
    Returns:
        None
    '''
    if not len(path) == 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


def ivar_variants_to_vcf(file_in, file_out, pass_only=False, min_allele_frequency=0, ignore_strand_bias=False, ignore_merge_codons=False):
    '''
    Description:
        Main function to convert iVar variants TSV to VCF.
    Input:
        file_in             : iVar variants TSV file
        file_out            : VCF output file
        pass_only           : Only keep variants that PASS filter [True, False]
        min_allele_freq     : Minimum allele frequency to keep a variant [0]
        ignore_strand_bias  : Do not apply strand-bias filter [True, False]
        ignore_merge_codons : Do not take into account consecutive positions belong to the same codon.
    Returns:
        None
    '''
    ## Create output directory
    filename = os.path.splitext(file_in)[0]
    out_dir = os.path.dirname(file_out)
    make_dir(out_dir)

    ## Define VCF header
    header_source = [
        "##fileformat=VCFv4.2",
        "##source=iVar"
    ]
    header_info = [
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">'
    ]
    header_filter = [
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##FILTER=<ID=ft,Description="Fisher\'s exact test of variant frequency compared to mean error rate, p-value > 0.05">'
    ]
    header_format = [
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=REF_DP,Number=1,Type=Integer,Description="Depth of reference base">',
        '##FORMAT=<ID=REF_RV,Number=1,Type=Integer,Description="Depth of reference base on reverse reads">',
        '##FORMAT=<ID=REF_QUAL,Number=1,Type=Integer,Description="Mean quality of reference base">',
        '##FORMAT=<ID=ALT_DP,Number=1,Type=Integer,Description="Depth of alternate base">',
        '##FORMAT=<ID=ALT_RV,Number=1,Type=Integer,Description="Depth of alternate base on reverse reads">',
        '##FORMAT=<ID=ALT_QUAL,Number=1,Type=Integer,Description="Mean quality of alternate base">',
        '##FORMAT=<ID=ALT_FREQ,Number=1,Type=Float,Description="Frequency of alternate base">',
    ]
    header_cols = [
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{filename}"
    ]
    if not ignore_strand_bias:
        header_info += [
            '##INFO=<ID=SB_PV,Number=1,Type=Float,Description="Strand-bias fisher-test p-value">'
        ]
        header_filter += [
            '##FILTER=<ID=sb,Description="Strand-bias fisher-test p-value < 0.05">'
        ]
    header = header_source + header_info + header_filter + header_format + header_cols

    ## Initialise variables
    var_list = []
    var_count_dict = {"SNP": 0, "INS": 0, "DEL": 0}
    dict_lines = {'CHROM':[], 'POS':[], 'ID':[], 'REF':[], 'ALT':[], 'REF_DP':[], 'REF_RV':[], 'ALT_DP':[], 'ALT_RV':[], 'QUAL':[], 'REF_CODON':[], 'ALT_CODON':[], 'FILTER': [], 'INFO':[], 'FORMAT':[], 'SAMPLE':[]}
    write_line = False
    fout = open(file_out, "w")
    fout.write('\n'.join(header) + '\n')
    with open(file_in, 'r') as fin:
        for line in fin:
            if not re.match("REGION", line):
                line = re.split("\t", line)

                ## Assign intial fields to variables
                CHROM = line[0]
                POS = line[1]
                ID = "."
                REF = line[2]
                ALT = line[3]

                ## REF/ALF depths
                REF_DP = int(line[4])
                REF_RV = int(line[5])
                REF_FW = REF_DP - REF_RV
                ALT_RV = int(line[8])
                ALT_DP = int(line[7])
                ALT_FW = ALT_DP - ALT_RV

                ## Perform a fisher_exact test for strand bias detection
                table = np.array([[REF_FW, REF_RV], [ALT_FW, ALT_RV]])
                oddsr, pvalue = fisher_exact(table, alternative='greater')

                ## Determine variant type
                var_type = "SNP"
                if ALT[0] == "+":
                    ALT = REF + ALT[1:]
                    var_type = "INS"
                elif ALT[0] == "-":
                    REF += ALT[1:]
                    ALT = line[2]
                    var_type = "DEL"

                QUAL = "."

                ## Determine FILTER field
                INFO = f"DP={line[11]}"
                pass_test = line[13]
                if ignore_strand_bias:
                    if pass_test == "TRUE":
                        FILTER = "PASS"
                    else:
                        FILTER = "ft"
                else:
                    ## Add SB in the FILTER field if strand-bias p-value is significant
                    if pvalue < 0.05 and pass_test == "TRUE":
                        FILTER = "sb"
                    elif pvalue > 0.05 and pass_test == "TRUE":
                        FILTER = "PASS"
                    elif pvalue <= 0.05 and pass_test == "FALSE":
                        FILTER = "ft;sb"
                    else:
                        FILTER = "ft"
                    INFO += f":SB_PV={str(round(pvalue, 5))}"

                FORMAT = "GT:REF_DP:REF_RV:REF_QUAL:ALT_DP:ALT_RV:ALT_QUAL:ALT_FREQ"
                SAMPLE = f'1:{":".join(line[4:11])}'

                REF_CODON = line[15]
                ALT_CODON = line[17]
                param_list = [CHROM, POS, ID, REF, ALT, REF_DP, REF_RV, ALT_DP, ALT_RV, QUAL, REF_CODON, ALT_CODON, FILTER, INFO, FORMAT, SAMPLE]

                if ignore_merge_codons or var_type != "SNP":
                    write_line = True
                    oline = (CHROM + "\t" + POS + "\t" + ID + "\t" + REF + "\t" + ALT + "\t" + QUAL + "\t" + FILTER + "\t" + INFO + "\t" + FORMAT + "\t" + SAMPLE + "\n")

                else:
                    ## dict_lines contains all the informative fields for 3 positions in the vcf.
                    # dict_lines has a maximum size of three.

                    ## Always fill dict_lines until size 2.
                    if len(dict_lines["POS"]) == 0 or len(dict_lines["POS"]) == 1:
                        for i,j in enumerate(dict_lines):
                            dict_lines.setdefault(j, []).append(param_list[i])
                        write_line=False

                    # If queue has size 2, we include the third line
                    elif len(dict_lines["POS"]) == 2:
                        for i,j in enumerate(dict_lines):
                            dict_lines.setdefault(j, []).append(param_list[i])
                        # Are two positions in the dict consecutive?
                        if check_consecutive(dict_lines["POS"]) == 2:
                            ## If the first position is not on the third position of the codon they are in the same codon.
                            if codon_position(dict_lines["REF_CODON"][0],dict_lines["ALT_CODON"][0]) != 2:
                                write_line = True
                                num_collapse = "2"
                                CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE = rename_vars(dict_lines, num_collapse)
                                oline = (CHROM + "\t" + POS + "\t" + ID + "\t" + REF + "\t" + ALT + "\t" + QUAL + "\t" + FILTER + "\t" + INFO + "\t" + FORMAT + "\t" + SAMPLE + "\n")
                                ## We removed the first two items in dict_lines with have been just processed.
                                for i,j in enumerate(dict_lines):
                                    dict_lines[list(dict_lines.keys())[i]].pop(0)
                                    dict_lines[list(dict_lines.keys())[i]].pop(0)
                            else:
                                write_line = True
                                oline =(dict_lines["CHROM"][0] + "\t" + dict_lines["POS"][0] + "\t" + dict_lines["ID"][0] + "\t" + dict_lines["REF"][0] + "\t" + dict_lines["ALT"][0] + "\t" + dict_lines["QUAL"][0] + "\t" + dict_lines["FILTER"][0] + "\t" + dict_lines["INFO"][0] + "\t" + dict_lines["FORMAT"][0] + "\t" + dict_lines["SAMPLE"][0] + "\n")
                                for i,j in enumerate(dict_lines):
                                    dict_lines[list(dict_lines.keys())[i]].pop(0)

                        # Are the three positions in the dict consecutive?
                        elif check_consecutive(dict_lines["POS"]) == 3:
                            ## we check the first position in which codon position is to process it acordingly.
                            # If first position is in the first codon position all three positions belong to the same codon.
                            if codon_position(dict_lines["REF_CODON"][0], dict_lines["ALT_CODON"][0]) == 0:
                                write_line = True
                                num_collapse = 3
                                CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE = rename_vars(dict_lines, num_collapse)
                                oline = (CHROM + "\t" + POS + "\t" + ID + "\t" + REF + "\t" + ALT + "\t" + QUAL + "\t" + FILTER + "\t" + INFO + "\t" + FORMAT + "\t" + SAMPLE + "\n")
                                for i,j in enumerate(dict_lines):
                                    dict_lines[list(dict_lines.keys())[i]].pop(0)
                                    dict_lines[list(dict_lines.keys())[i]].pop(0)
                                # we empty the dict_lines
                                dict_lines = {'CHROM':[], 'POS':[], 'ID':[], 'REF':[], 'ALT':[], 'REF_DP':[], 'REF_RV':[], 'ALT_DP':[], 'ALT_RV':[], 'QUAL':[], 'REF_CODON':[], 'ALT_CODON':[], 'FILTER':[], 'INFO':[], 'FORMAT':[], 'SAMPLE':[]}
                            # If first position is in the second codon position, we have the two first positions belonging to the same codon and the last one independent.
                            elif codon_position(dict_lines["REF_CODON"][0], dict_lines["ALT_CODON"][0]) == 1:
                                write_line = True
                                num_collapse = 2
                                CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE = rename_vars(dict_lines, num_collapse)
                                oline = (CHROM + "\t" + POS + "\t" + ID + "\t" + REF + "\t" + ALT + "\t" + QUAL + "\t" + FILTER + "\t" + INFO + "\t" + FORMAT + "\t" + SAMPLE + "\n")
                                for i,j in enumerate(dict_lines):
                                    dict_lines[list(dict_lines.keys())[i]].pop(0)
                                    dict_lines[list(dict_lines.keys())[i]].pop(0)
                            ## Finally if we have the first position in the last codon position, we write first position and left the remaining two to be evaluated in the next iteration.
                            elif codon_position(dict_lines["REF_CODON"][0], dict_lines["ALT_CODON"][0]) == 2:
                                write_line = True
                                oline =(dict_lines["CHROM"][0] + "\t" + dict_lines["POS"][0] + "\t" + dict_lines["ID"][0] + "\t" + dict_lines["REF"][0] + "\t" + dict_lines["ALT"][0] + "\t" + dict_lines["QUAL"][0] + "\t" + dict_lines["FILTER"][0] + "\t" + dict_lines["INFO"][0] + "\t" + dict_lines["FORMAT"][0] + "\t" + dict_lines["SAMPLE"][0] + "\n")
                                for i,j in enumerate(dict_lines):
                                    dict_lines[list(dict_lines.keys())[i]].pop(0)

                        elif check_consecutive(dict_lines["POS"]) == False:
                            write_line = True
                            oline =(dict_lines["CHROM"][0] + "\t" + dict_lines["POS"][0] + "\t" + dict_lines["ID"][0] + "\t" + dict_lines["REF"][0] + "\t" + dict_lines["ALT"][0] + "\t" + dict_lines["QUAL"][0] + "\t" + dict_lines["FILTER"][0] + "\t" + dict_lines["INFO"][0] + "\t" + dict_lines["FORMAT"][0] + "\t" + dict_lines["SAMPLE"][0] + "\n")
                            for i,j in enumerate(dict_lines):
                                dict_lines[list(dict_lines.keys())[i]].pop(0)
                    else:
                        print("Something went terribly wrong!!" + str(len(dict_lines["POS"])))

                ## Determine whether to output variant
                if pass_only and FILTER != "PASS":
                    write_line = False
                if float(line[10]) < min_allele_frequency:
                    write_line = False
                if (CHROM, POS, REF, ALT) in var_list:
                    write_line = False
                else:
                    var_list.append((CHROM, POS, REF, ALT))

                ## Write to file
                if write_line:
                    var_count_dict[var_type] += 1
                    fout.write(oline)

    ## Print variant counts to pass to MultiQC
    var_count_list = [(k, str(v)) for k, v in sorted(var_count_dict.items())]
    print("\t".join(["sample"] + [x[0] for x in var_count_list]))
    print("\t".join([filename] + [x[1] for x in var_count_list]))

    ## Handle last 3 lines.
    if  len(dict_lines["POS"]) == 2:
        if check_consecutive(dict_lines["POS"]) == 2:
            if codon_position(dict_lines["REF_CODON"][0],dict_lines["ALT_CODON"][0]) != 2:
                write_line = True
                num_collapse = 2
                CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE = rename_vars(dict_lines, num_collapse)
                oline = (CHROM + "\t" + POS + "\t" + ID + "\t" + REF + "\t" + ALT + "\t" + QUAL + "\t" + FILTER + "\t" + INFO + "\t" + FORMAT + "\t" + SAMPLE + "\n")
                fout.write(oline)
        else:
            oline = (dict_lines["CHROM"][0] + "\t" + dict_lines["POS"][0] + "\t" + dict_lines["ID"][0] + "\t" + dict_lines["REF"][0] + "\t" + dict_lines["ALT"][0] + "\t" + dict_lines["QUAL"][0] + "\t" + dict_lines["FILTER"][0] + "\t" + dict_lines["INFO"][0] + "\t" + dict_lines["FORMAT"][0] + "\t" + dict_lines["SAMPLE"][0] + "\n")
            oline1 = (dict_lines["CHROM"][1] + "\t" + dict_lines["POS"][1] + "\t" + dict_lines["ID"][1] + "\t" + dict_lines["REF"][1] + "\t" + dict_lines["ALT"][1] + "\t" + dict_lines["QUAL"][1] + "\t" + dict_lines["FILTER"][1] + "\t" + dict_lines["INFO"][1] + "\t" + dict_lines["FORMAT"][1] + "\t" + dict_lines["SAMPLE"][1] + "\n")
            fout.write(oline)
            fout.write(oline1)
    elif len(dict_lines["POS"]) == 1:
        oline =(dict_lines["CHROM"][0] + "\t" + dict_lines["POS"][0] + "\t" + dict_lines["ID"][0] + "\t" + dict_lines["REF"][0] + "\t" + dict_lines["ALT"][0] + "\t" + dict_lines["QUAL"][0] + "\t" + dict_lines["FILTER"][0] + "\t" + dict_lines["INFO"][0] + "\t" + dict_lines["FORMAT"][0] + "\t" + dict_lines["SAMPLE"][0] + "\n")
        fout.write(oline)
    fout.close()


def main(args=None):
    args = parse_args(args)
    ivar_variants_to_vcf(
        args.file_in,
        args.file_out,
        args.pass_only,
        args.allele_freq_threshold,
        args.ignore_strand_bias,
        args.ignore_merge_codons,
    )


if __name__ == "__main__":
    sys.exit(main())
