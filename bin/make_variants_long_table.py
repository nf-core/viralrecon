#!/usr/bin/env python

import os
import sys
import glob
import errno
import shutil
import logging
import argparse
import pandas as pd
from matplotlib import table


logger = logging.getLogger()


pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", None)


def parser_args(args=None):
    Description = "Create long/wide tables containing variant information."
    Epilog = """Example usage: python make_variants_long_table.py --bcftools_query_dir ./bcftools_query/ --snpsift_dir ./snpsift/ --pangolin_dir ./pangolin/"""
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument(
        "-bd",
        "--bcftools_query_dir",
        type=str,
        default="./bcftools_query",
        help="Directory containing output of BCFTools query for each sample (default: './bcftools_query').",
    )
    parser.add_argument(
        "-sd",
        "--snpsift_dir",
        type=str,
        default="./snpsift",
        help="Directory containing output of SnpSift for each sample (default: './snpsift').",
    )
    parser.add_argument(
        "-pd",
        "--pangolin_dir",
        type=str,
        default="./pangolin",
        help="Directory containing output of Pangolin for each sample (default: './pangolin').",
    )
    parser.add_argument(
        "-bs",
        "--bcftools_file_suffix",
        type=str,
        default=".bcftools_query.txt",
        help="Suffix to trim off BCFTools query file name to obtain sample name (default: '.bcftools_query.txt').",
    )
    parser.add_argument(
        "-ss",
        "--snpsift_file_suffix",
        type=str,
        default=".snpsift.txt",
        help="Suffix to trim off SnpSift file name to obtain sample name (default: '.snpsift.txt').",
    )
    parser.add_argument(
        "-ps",
        "--pangolin_file_suffix",
        type=str,
        default=".pangolin.csv",
        help="Suffix to trim off Pangolin file name to obtain sample name (default: '.pangolin.csv').",
    )
    parser.add_argument(
        "-of",
        "--output_file",
        type=str,
        default="variants_long_table.csv",
        help="Full path to output file (default: 'variants_long_table.csv').",
    )
    parser.add_argument(
        "-vc", "--variant_caller", type=str, default="ivar", help="Tool used to call the variants (default: 'ivar')."
    )
    return parser.parse_args(args)


def make_dir(path):
    if not len(path) == 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


def get_file_dict(file_dir, file_suffix):
    files = glob.glob(os.path.join(file_dir, f"*{file_suffix}"))
    samples = [os.path.basename(x).removesuffix(f"{file_suffix}") for x in files]

    return dict(zip(samples, files))


def three_letter_aa_to_one(hgvs_three):
    aa_dict = {
        "Ala": "A",
        "Arg": "R",
        "Asn": "N",
        "Asp": "D",
        "Cys": "C",
        "Gln": "Q",
        "Glu": "E",
        "Gly": "G",
        "His": "H",
        "Ile": "I",
        "Leu": "L",
        "Lys": "K",
        "Met": "M",
        "Phe": "F",
        "Pro": "P",
        "Pyl": "O",
        "Ser": "S",
        "Sec": "U",
        "Thr": "T",
        "Trp": "W",
        "Tyr": "Y",
        "Val": "V",
        "Asx": "B",
        "Glx": "Z",
        "Xaa": "X",
        "Xle": "J",
        "Ter": "*",
    }
    hgvs_one = hgvs_three
    for key in aa_dict:
        if key in hgvs_one:
            hgvs_one = hgvs_one.replace(str(key), str(aa_dict[key]))

    return hgvs_one


## Returns a pandas dataframe in the format:
#         CHROM   POS  REF  ALT FILTER   DP  REF_DP  ALT_DP    AF
# 0  MN908947.3   241    C    T   PASS  642     375     266  0.41
# 1  MN908947.3  1875    C    T   PASS   99      63      34  0.34
def ivar_bcftools_query_to_table(bcftools_query_file):
    table = pd.read_table(bcftools_query_file, header="infer")
    table = table.dropna(how="all", axis=1)
    old_colnames = list(table.columns)
    new_colnames = [x.split("]")[-1].split(":")[-1] for x in old_colnames]
    table.rename(columns=dict(zip(old_colnames, new_colnames)), inplace=True)

    if not table.empty:
        table[["ALT_DP", "DP"]] = table[["ALT_DP", "DP"]].apply(pd.to_numeric)
        table["AF"] = table["ALT_DP"] / table["DP"]
        table["AF"] = table["AF"].round(2)

    return table


## Returns a pandas dataframe in the format:
#         CHROM    POS REF ALT FILTER  DP REF_DP  ALT_DP    AF
# 0  MN908947.3    241   C   T      .  24      8      16  0.67
# 1  MN908947.3   3037   C   T      .  17      5      12  0.71
def bcftools_bcftools_query_to_table(bcftools_query_file):
    table = pd.read_table(bcftools_query_file, header="infer")
    table = table.dropna(how="all", axis=1)
    old_colnames = list(table.columns)
    new_colnames = [x.split("]")[-1].split(":")[-1] for x in old_colnames]
    table.rename(columns=dict(zip(old_colnames, new_colnames)), inplace=True)

    if not table.empty:
        table[["REF_DP", "ALT_DP"]] = table["AD"].str.split(",", expand=True)
        table[["ALT_DP", "DP"]] = table[["ALT_DP", "DP"]].apply(pd.to_numeric)
        table["AF"] = table["ALT_DP"] / table["DP"]
        table["AF"] = table["AF"].round(2)
        table.drop("AD", axis=1, inplace=True)

    return table


## Returns a pandas dataframe in the format:
#         CHROM   POS REF ALT FILTER  DP  REF_DP  ALT_DP    AF
# 0  MN908947.3   241   C   T   PASS  30       1      29  0.97
# 1  MN908947.3  1163   A   T   PASS  28       0      28  1.00
def nanopolish_bcftools_query_to_table(bcftools_query_file):
    table = pd.read_table(bcftools_query_file, header="infer")
    table = table.dropna(how="all", axis=1)
    old_colnames = list(table.columns)
    new_colnames = [x.split("]")[-1].split(":")[-1] for x in old_colnames]
    table.rename(columns=dict(zip(old_colnames, new_colnames)), inplace=True)

    ## Split out ref/alt depths from StrandSupport column
    if not table.empty:
        table_cp = table.copy()
        table_cp[["FORW_REF_DP", "REV_REF_DP", "FORW_ALT_DP", "REV_ALT_DP"]] = table_cp["StrandSupport"].str.split(
            ",", expand=True
        )
        table_cp[["FORW_REF_DP", "REV_REF_DP", "FORW_ALT_DP", "REV_ALT_DP"]] = table_cp[
            ["FORW_REF_DP", "REV_REF_DP", "FORW_ALT_DP", "REV_ALT_DP"]
        ].apply(pd.to_numeric)

        table["DP"] = table_cp[["FORW_REF_DP", "REV_REF_DP", "FORW_ALT_DP", "REV_ALT_DP"]].sum(axis=1)
        table["REF_DP"] = table_cp[["FORW_REF_DP", "REV_REF_DP"]].sum(axis=1)
        table["ALT_DP"] = table_cp[["FORW_ALT_DP", "REV_ALT_DP"]].sum(axis=1)
        table["AF"] = table["ALT_DP"] / table["DP"]
        table["AF"] = table["AF"].round(2)
        table.drop("StrandSupport", axis=1, inplace=True)

    return table


## Returns a pandas dataframe in the format:
#         CHROM    POS REF ALT FILTER  DP REF_DP  ALT_DP    AF
# 0  MN908947.3    241   C   T   PASS  21      0      21  1.00
# 1  MN908947.3   3037   C   T   PASS  28      0      25  0.89
def medaka_bcftools_query_to_table(bcftools_query_file):
    table = pd.read_table(bcftools_query_file, header="infer")
    table = table.dropna(how="all", axis=1)
    old_colnames = list(table.columns)
    new_colnames = [x.split("]")[-1].split(":")[-1] for x in old_colnames]
    table.rename(columns=dict(zip(old_colnames, new_colnames)), inplace=True)

    if not table.empty:
        table[["REF_DP", "ALT_DP"]] = table["AC"].str.split(",", expand=True)
        table[["ALT_DP", "DP"]] = table[["ALT_DP", "DP"]].apply(pd.to_numeric)
        table["AF"] = table["ALT_DP"] / table["DP"]
        table["AF"] = table["AF"].round(2)
        table.drop("AC", axis=1, inplace=True)

    return table


def get_pangolin_lineage(pangolin_file):
    table = pd.read_csv(pangolin_file, sep=",", header="infer")

    return table["lineage"][0]


def snpsift_to_table(snpsift_file):
    table = pd.read_table(snpsift_file, sep="\t", header="infer")
    table = table.loc[:, ~table.columns.str.contains("^Unnamed")]
    old_colnames = list(table.columns)
    new_colnames = [x.replace("ANN[*].", "") for x in old_colnames]
    table.rename(columns=dict(zip(old_colnames, new_colnames)), inplace=True)
    table = table.loc[:, ["CHROM", "POS", "REF", "ALT", "GENE", "EFFECT", "HGVS_C", "HGVS_P"]]

    ## Split by comma and get first value in cols = ['ALT','GENE','EFFECT','HGVS_C','HGVS_P']
    for i in range(len(table)):
        for j in range(3, 8):
            table.iloc[i, j] = str(table.iloc[i, j]).split(",")[0]

    ## Amino acid substitution
    aa = []
    for index, item in table["HGVS_P"].iteritems():
        hgvs_p = three_letter_aa_to_one(str(item))
        aa.append(hgvs_p)
    table["HGVS_P_1LETTER"] = pd.Series(aa)

    return table


def main(args=None):
    args = parser_args(args)

    ## Create output directory if it doesn't exist
    out_dir = os.path.dirname(args.output_file)
    make_dir(out_dir)

    ## Check correct variant caller has been provided
    variant_callers = ["ivar", "bcftools", "nanopolish", "medaka"]
    if args.variant_caller not in variant_callers:
        logger.error(
            f"Invalid option '--variant caller {args.variant_caller}'. Valid options: " + ", ".join(variant_callers)
        )
        sys.exit(1)

    ## Find files and create a dictionary {'sample': '/path/to/file'}
    bcftools_files = get_file_dict(args.bcftools_query_dir, args.bcftools_file_suffix)
    snpsift_files = get_file_dict(args.snpsift_dir, args.snpsift_file_suffix)
    pangolin_files = get_file_dict(args.pangolin_dir, args.pangolin_file_suffix)

    ## Check all files are provided for each sample
    if set(bcftools_files) != set(snpsift_files):
        logger.error(
            f"Number of BCFTools ({len(bcftools_files)}) and SnpSift ({len(snpsift_files)}) files do not match!"
        )
        sys.exit(1)
    else:
        if pangolin_files:
            if set(bcftools_files) != set(pangolin_files):
                logger.error(
                    f"Number of BCFTools ({len(bcftools_files)}) and Pangolin ({len(pangolin_files)}) files do not match!"
                )
                sys.exit(1)

    ## Create per-sample table and write to file
    sample_tables = []
    for sample in sorted(bcftools_files):
        ## Read in BCFTools query file
        bcftools_table = None
        if args.variant_caller == "ivar":
            bcftools_table = ivar_bcftools_query_to_table(bcftools_files[sample])
        elif args.variant_caller == "bcftools":
            bcftools_table = bcftools_bcftools_query_to_table(bcftools_files[sample])
        elif args.variant_caller == "nanopolish":
            bcftools_table = nanopolish_bcftools_query_to_table(bcftools_files[sample])
        elif args.variant_caller == "medaka":
            bcftools_table = medaka_bcftools_query_to_table(bcftools_files[sample])

        if not bcftools_table.empty:
            ## Read in SnpSift file
            snpsift_table = snpsift_to_table(snpsift_files[sample])

            merged_table = pd.DataFrame(data=bcftools_table)
            merged_table.insert(0, "SAMPLE", sample)
            merged_table = pd.merge(merged_table, snpsift_table, how="outer")
            merged_table["CALLER"] = args.variant_caller

            ## Read in Pangolin lineage file
            if pangolin_files:
                merged_table["LINEAGE"] = get_pangolin_lineage(pangolin_files[sample])

            sample_tables.append(merged_table)

    ## Merge table across samples
    if sample_tables:
        merged_tables = pd.concat(sample_tables)
        merged_tables.to_csv(args.output_file, index=False, encoding="utf-8-sig")


if __name__ == "__main__":
    sys.exit(main())
