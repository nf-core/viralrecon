#!/usr/bin/env python

from email.charset import QP
import os
import sys
import re
import errno
import argparse
from collections import OrderedDict
from collections import deque

import numpy as np
from Bio import SeqIO
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
        action="store_true",
    )
    parser.add_argument(
        "-ic",
        "--ignore_merge_codons",
        help="Output variants without taking into account if consecutive positions belong to the same codon.",
        action="store_true",
    )
    parser.add_argument(
        "-f",
        "--fasta",
        type=str,
        default=None,
        help="Fasta file used in mapping and variant calling for vcf header reference genome lenght info.",
    )
    return parser.parse_args(args)


def make_dir(path):
    """
    Description:
        Create directory if it doesn't exist.
    Input:
        path - path where the directory will be created.
    Returns:
        None
    """
    if not len(path) == 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


def parse_ivar_line(line):
    """
    Description:
        Parse ivar line to get needed variables for vcf format.
    input:
        line - ivar tsv line
    return:
        CHROM, POS, ID, REF, ALT, QUAL, INFO, FORMAT, REF_CODON, ALT_CODON, pass_test, var_type
    """

    line = line.strip("\n").split("\t")
    ## Assign intial fields to variables
    CHROM = line[0]
    POS = line[1]
    ID = "."
    REF = line[2]
    ALT = line[3]

    ## REF/ALF depths and quals
    try:
        REF_DP = int(line[4])
    except ValueError:
        print(line)
        print(line[4])
        exit(-1)
    REF_RV = int(line[5])
    REF_FW = REF_DP - REF_RV
    REF_QUAL = int(line[6])
    ALT_RV = int(line[8])
    ALT_DP = int(line[7])
    ALT_FW = ALT_DP - ALT_RV
    ALT_QUAL = int(line[9])
    ALT_FREQ = float(line[10])
    FORMAT = [REF_DP, REF_RV, REF_QUAL, ALT_DP, ALT_RV, ALT_QUAL, ALT_FREQ]

    ## Codon annotation
    REF_CODON = line[15]
    ALT_CODON = line[17]

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
    INFO = f"DP={int(float(line[11]))}"
    pass_test = line[13]

    return (
        CHROM,
        POS,
        ID,
        REF,
        ALT,
        QUAL,
        INFO,
        FORMAT,
        REF_CODON,
        ALT_CODON,
        pass_test,
        var_type,
    )


######################
## FILTER FUNCTIONS ##
######################


def ivar_filter(pass_test):
    """
    Description:
        process ivar filter into vcf filter format.
    input:
        pass_test - ivar fisher exact test [ True, False ]
    return:
        Whether it passes the filter or not. [False, "ft"]
    """
    if pass_test == "TRUE":
        return False
    else:
        return "ft"


def strand_bias_filter(format):
    """
    Description:
        Calculate strand-bias fisher test.
    input:
        format - format variables
    return:
        Whether it passes the filter or not. [False, "sb"]
    """
    # format=[REF_DP, REF_RV, REF_QUAL, ALT_DP, ALT_RV, ALT_QUAL, ALT_FREQ]
    # table:
    ##  REF_FW  REF_RV
    ##  ALT_FW  ALT_RV
    table = np.array([[format[0] - format[1], format[1]], [format[3] - format[4], format[4]]])
    oddsr, pvalue = fisher_exact(table, alternative="greater")

    # h0: both strands are equally represented.
    # If test is significant h0 is refused so there is an strand bias.
    if pvalue < 0.05:
        return "sb"
    else:
        return False


def write_vcf_header(ref, ignore_strand_bias, file_out, filename):
    """
    Description:
        Write vcf header for VCFv4.2
    input:
        ref - (optional), ref in fasta format
        ignore_strand_bias - if no strand-bias is calculated [True, False]
        file_out - output file_in
        filename - name of the output file
    return:
        Nothing.
    """
    ## Define VCF header
    header_source = ["##fileformat=VCFv4.2", "##source=iVar"]
    if ref:
        header_contig = []
        for record in SeqIO.parse(ref, "fasta"):
            header_contig += ["##contig=<ID=" + record.id + ",length=" + str(len(record.seq)) + ">"]

        header_source += header_contig

    header_info = ['##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">']
    header_filter = [
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##FILTER=<ID=ft,Description="Fisher\'s exact test of variant frequency compared to mean error rate, p-value > 0.05">',
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
    header_cols = [f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{filename}"]
    if not ignore_strand_bias:
        header_filter += ['##FILTER=<ID=sb,Description="Strand-bias fisher-test p-value < 0.05">']

    header = header_source + header_info + header_filter + header_format + header_cols
    fout = open(file_out, "w")
    fout.write("\n".join(header) + "\n")
    fout.close()


def write_vcf_line(chrom, pos, id, ref, alt, filter, qual, info, format, file_out):
    """
    Description:
        Format variables into vcf line format and write line to file.
    input:
        chrom, pos, id, ref, alt, filter, qual, info, format - vcf variables
        file_out                                             - file output
    return:
        Nothing.
    """
    sample = f'1:{":".join(str(x) for x in format)}'
    format = "GT:REF_DP:REF_RV:REF_QUAL:ALT_DP:ALT_RV:ALT_QUAL:ALT_FREQ"

    oline = (
        chrom
        + "\t"
        + pos
        + "\t"
        + id
        + "\t"
        + ref
        + "\t"
        + alt
        + "\t"
        + qual
        + "\t"
        + filter
        + "\t"
        + info
        + "\t"
        + format
        + "\t"
        + sample
        + "\n"
    )
    fout = open(file_out, "a")
    fout.write(oline)
    fout.close()


############################
## MERGE CODONS FUNCTIONS ##
############################


def check_consecutive(mylist):
    """
    Description:
        This function checks a list of  numbers and returns how many items are consecutive.
    input:
        my_list - A list of integers
    return:
        Number of items consecutive in the list - [False, 2, 3,..]
    """
    # getting first index of tuple for consecutive checking
    my_list = list(map(int, [i[0] for i in mylist]))
    ## Check if the list contains consecutive numbers
    if len(my_list) == 1:
        return False
    elif sorted(my_list) == list(range(min(my_list), max(my_list) + 1)):
        return len(my_list)
    else:
        ## If not, and the list is > 1, remove the last item and reevaluate.
        if len(my_list) > 2:
            my_list.pop()
            if sorted(my_list) == list(range(min(my_list), max(my_list) + 1)):
                return len(my_list)
        else:
            return False
        return False


def get_diff_position(seq1, seq2):
    """
    Description:
        Function to compare two codon nucleotide sequences (size 3) and retuns the position where it differs.
    Input:
        seq1 - string size 3 [A,T,C,G]. Ex. "ATC"
        seq2 - string size 3 [A,T,C,G]. Ex. "ACC"
    Returns:
        Returns position where seq1 != seq2
    """
    # If codon is NA treat as not same codon
    if seq1 == "NA":
        return 2

    ind_diff = [i for i in range(len(seq1)) if seq1[i] != seq2[i]]
    if len(ind_diff) > 1:
        print("There has been an issue, more than one difference between the seqs.")
        return False
    else:
        return ind_diff[0]


def check_merge_codons(q_pos, fe_codon_ref, fe_codon_alt):
    """
    Description:
        Logic for determine if variant lines need to be collapsed into one determining
        if they are consecutive and belong to the same codon.
    Input:
        qpos         - list of positions. Ex. [4441, 4442, 4443]
        fe_codon_ref - first position codon annotation for ref. Ex. "ATG"
        fe_codon_alt - first position codon annotation for alt. Ex. "AGG"
    Returns:
        Returns num_collapse. Number of lines that need to be collapsed into one.
    """
    # Are two positions in the queue consecutive?
    # q_pos = [4441, 4442, 5067]
    num_collapse = 0
    if check_consecutive(list(q_pos)) == 2:
        ## If the first position is not on the third position of the codon they are in the same codon.
        if get_diff_position(fe_codon_ref, fe_codon_alt) != 2:
            num_collapse = 2
        else:
            num_collapse = 1
    # Are the three positions in the queue consecutive?
    # q_pos = [4441, 4442, 4443]
    elif check_consecutive(list(q_pos)) == 3:
        ## we check the first position in which codon position is to process it acordingly.
        # If first position is in the first codon position all three positions belong to the same codon.
        if get_diff_position(fe_codon_ref, fe_codon_alt) == 0:
            num_collapse = 3
        # If first position is in the second codon position, we have the two first positions belonging to the same codon and the last one independent.
        elif get_diff_position(fe_codon_ref, fe_codon_alt) == 1:
            num_collapse = 2
        ## Finally if we have the first position in the last codon position, we write first position and left the remaining two to be evaluated in the next iteration.
        elif get_diff_position(fe_codon_ref, fe_codon_alt) == 2:
            num_collapse = 1
    # If no consecutive process only one line.
    elif check_consecutive(list(q_pos)) == False:
        num_collapse = 1

    return num_collapse


def process_variants(variants, num_collapse):
    """
    Description:
        The function set the variables acordingly to the lines to collapse do to consecutive variants.
    Input:
        variants - Dict with var lines.
        num_collapse - number of lines to collapse [2,3]
    Returns::
        Vars fixed: chrom, pos, id, ref, alt, qual, filter, info, format
    """
    # Collapsed variant parameters equal to first variant
    key_list = ["chrom", "pos", "id", "qual", "filter", "info", "format"]
    chrom, pos, id, qual, filter, info, format = [variants[next(iter(variants))][key] for key in key_list]

    # If no consecutive, process one variant line
    # If two consecutive, process two variant lines into one
    # If three consecutive process three variant lines and write one
    ref = ""
    alt = ""
    iter_variants = iter(variants)
    for _ in range(num_collapse):  # fixed notation
        var = next(iter_variants)
        ref += variants[var]["ref"]
        alt += variants[var]["alt"]

    return chrom, pos, id, ref, alt, qual, filter, info, format


def main(args=None):
    # Process args
    args = parse_args(args)

    # Initialize vars
    filename = os.path.splitext(args.file_in)[0]
    out_dir = os.path.dirname(args.file_out)
    var_list = []  # store variants
    var_count_dict = {"SNP": 0, "INS": 0, "DEL": 0}  # variant counts
    variants = OrderedDict()  # variant dict (merge codon)
    q_pos = deque([], maxlen=3)  # pos fifo queue (merge codon)
    last_pos = ""

    # Create output directory
    make_dir(out_dir)

    ##############################
    ## Write vcf header to file ##
    ##############################
    write_vcf_header(args.fasta, args.ignore_strand_bias, args.file_out, filename)

    #################################
    ## Read and process input file ##
    #################################
    with open(args.file_in, "r") as fin:
        for line in fin:
            if "REGION" not in line:
                ################
                ## Parse line ##
                ################
                ## format=
                # [REF_DP, REF_RV, REF_QUAL, ALT_DP, ALT_RV, ALT_QUAL, ALT_FREQ]
                write_line = True
                (
                    chrom,
                    pos,
                    id,
                    ref,
                    alt,
                    qual,
                    info,
                    format,
                    ref_codon,
                    alt_codon,
                    pass_test,
                    var_type,
                ) = parse_ivar_line(line)

                ## If pos is duplicated due to annotation skip lines
                if pos == last_pos:
                    continue

                last_pos = pos
                #####################
                ## Process filters ##
                #####################
                ## ivar fisher test
                filter = ""
                if ivar_filter(pass_test):
                    filter = ivar_filter(pass_test)
                ## strand-bias fisher test
                if not args.ignore_strand_bias:
                    if strand_bias_filter(format):
                        if filter:
                            filter += ";" + strand_bias_filter(format)
                        else:
                            filter = strand_bias_filter(format)

                if not filter:
                    filter = "PASS"

                #####################
                ## Filter variants ##
                #####################
                if args.pass_only and filter != "PASS":
                    write_line = False
                ### AF filtering. ALT_DP/(ALT_DP+REF_DP)
                if float(format[3] / (format[0] + format[3])) < args.allele_freq_threshold:
                    write_line = False
                ### Duplication filter
                if (chrom, pos, ref, alt) in var_list:
                    write_line = False
                else:
                    var_list.append((chrom, pos, ref, alt))

                ############################################################
                ##                MERGE_CODONS                            ##
                ## Merge consecutive variants belonging to the same codon ##
                ############################################################
                if not args.ignore_merge_codons and var_type == "SNP":
                    ## re-fill queue and dict accordingly
                    q_pos.append((pos, var_type))  # adding type information
                    variants[(chrom, pos, ref, alt)] = {
                        "chrom": chrom,
                        "pos": pos,
                        "id": id,
                        "ref": ref,
                        "alt": alt,
                        "qual": qual,
                        "filter": filter,
                        "info": info,
                        "format": format,
                        "ref_codon": ref_codon,
                        "alt_codon": alt_codon,
                    }

                    if len(q_pos) == q_pos.maxlen:
                        fe_codon_ref = variants[next(iter(variants))]["ref_codon"]
                        fe_codon_alt = variants[next(iter(variants))]["alt_codon"]
                        num_collapse = check_merge_codons(q_pos, fe_codon_ref, fe_codon_alt)
                        (
                            chrom,
                            pos,
                            id,
                            ref,
                            alt,
                            qual,
                            filter,
                            info,
                            format,
                        ) = process_variants(variants, num_collapse)

                        ## Empty variants dict and queue accordingly
                        for _ in range(num_collapse):
                            variants.popitem(last=False)
                            q_pos.popleft()
                    else:
                        write_line = False

                ##############################
                ## Write output to vcf file ##
                ##############################
                if write_line:
                    var_count_dict[var_type] += 1
                    write_vcf_line(
                        chrom,
                        pos,
                        id,
                        ref,
                        alt,
                        filter,
                        qual,
                        info,
                        format,
                        args.file_out,
                    )

    if not args.ignore_merge_codons:
        #######################
        ## handle last lines ##
        #######################
        while len(q_pos) > 0:
            try:
                fe_codon_ref = variants[next(iter(variants))]["ref_codon"]
                fe_codon_alt = variants[next(iter(variants))]["alt_codon"]
            except StopIteration:
                break
            else:
                num_collapse = check_merge_codons(q_pos, fe_codon_ref, fe_codon_alt)
                (chrom, pos, id, ref, alt, qual, filter, info, format) = process_variants(variants, num_collapse)

                var_count_dict[q_pos[0][1]] += 1
                write_vcf_line(chrom, pos, id, ref, alt, filter, qual, info, format, args.file_out)
                ## Empty variants dict and queue accordingly
                for _ in range(num_collapse):
                    variants.popitem(last=False)
                    q_pos.popleft()

    #############################################
    ##  variant counts to pass to MultiQC      ##
    #############################################
    var_count_list = [(k, str(v)) for k, v in sorted(var_count_dict.items())]

    # format output table a little more cleanly
    # row_spacing = len(filename)

    row = create_f_string(30, "<")  # an arbitraily long value to fit most sample names
    row += create_f_string(10) * len(var_count_list)  # A spacing of ten looks pretty

    headers = ["sample"]
    headers.extend([x[0] for x in var_count_list])
    data = [filename]
    data.extend([x[1] for x in var_count_list])
    print(row.format(*headers))
    print(row.format(*data))


def create_f_string(str_size, placement="^"):
    row_size = "{: " + placement + str(str_size) + "}"
    return row_size


if __name__ == "__main__":
    sys.exit(main())
