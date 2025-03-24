#!/usr/bin/env python

import os
import sys
import errno
import argparse
import warnings
from itertools import product

import numpy as np
from Bio import SeqIO
from scipy.stats import fisher_exact
import pandas as pd


def parse_args(args=None):
    Description = "Convert iVar variants TSV file to VCF format."
    Epilog = """Example usage: python ivar_variants_to_vcf.py <file_in> <file_out>"""

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("file_in", help="Input iVar TSV file.")
    parser.add_argument("file_out", help="Full path to output VCF file.")
    parser.add_argument(
        "-po",
        "--pass_only",
        action="store_true",
        help="Only output variants that PASS filters.",
    )
    parser.add_argument(
        "-af",
        "--allele_freq_threshold",
        type=float,
        default=0,
        help="Only output variants where allele frequency is greater than this number (default: 0).",
    )
    parser.add_argument(
        "-bq",
        "--bad_quality_threshold",
        type=int,
        default=20,
        help="Only output variants where ALT_QUAL is greater than this number (default: 20).",
    )
    parser.add_argument(
        "-mt",
        "--merge_af_threshold",
        type=float,
        default=0.25,
        help="Only merge variants within a range of Allele Frequency. Useful to distinguish haplotypes.",
    )
    parser.add_argument(
        "-cf",
        "--consensus_af",
        type=float,
        default=0.75,
        help="Allele Frenquecy threshold used to include variants in consensus sequence.",
    )
    parser.add_argument(
        "-is",
        "--ignore_strand_bias",
        action="store_true",
        help="Does not take strand bias into account, use this option when using amplicon sequencing.",
    )
    parser.add_argument(
        "-ic",
        "--ignore_merge_codons",
        action="store_true",
        help="Output variants without taking into account if consecutive positions belong to the same codon.",
    )
    parser.add_argument(
        "-f",
        "--fasta",
        type=str,
        default=None,
        help="Fasta file used in mapping and variant calling for vcf header reference genome lenght info.",
    )
    return parser.parse_args(args)


class IvarVariants:
    def __init__(
        self,
        file_in=None,
        file_out=None,
        pass_only=None,
        freq_threshold=None,
        bad_qual_threshold=None,
        merge_af_threshold=None,
        consensus_af=None,
        ignore_stbias=None,
        ignore_merge=None,
        fasta=None,
    ):
        self.file_in = file_in
        self.file_out = file_out
        if os.path.exists(self.file_out):
            self.filename = str(os.path.splitext(self.file_in)[0])
        else:
            self.filename = str(self.file_in)
        self.pass_only = pass_only
        try:
            self.freq_threshold = float(freq_threshold)
        except Exception:
            print("Invalid allele frequency threshold. Setting it to 0")
            self.freq_threshold = 0
        try:
            self.bad_qual_threshold = float(bad_qual_threshold)
        except Exception:
            print("Invalid bad quality threshold. Setting it to 20")
            self.bad_qual_threshold = 20
        try:
            self.merge_af_threshold = float(merge_af_threshold)
        except Exception:
            print("Invalid merge af threshold. Setting it to 0.25")
            self.merge_af_threshold = 0.25
        try:
            self.consensus_af = float(consensus_af)
        except Exception:
            print("Invalid consensus af threshold. Setting it to 0.75")
            self.consensus_af = 0.75
        self.ignore_stbias = ignore_stbias
        self.ignore_merge = ignore_merge
        if not self.file_out:
            exit("Output file not provided. Aborting...")
        if fasta:
            self.ref_fasta = (
                fasta if os.path.exists(str(fasta)) else exit("Invalid fasta path")
            )
        else:
            self.ref_fasta = None

        if self.file_in:
            try:
                self.raw_ivar_df = pd.read_csv(self.file_in, sep="\t")
            except Exception:
                exit(f"Could not read input file: {file_in}")
        else:
            exit("Input file not provided. Aborting...")
        if self.raw_ivar_df.empty:
            exit("Input tsv was empty")

    def strand_bias_filter(self, row):
        """Calculate strand-bias fisher test.

        Args:
            row (pd.Series)- a row from ivar.tsv dataframe as series
        Returns:
            str(): Whether it passes the filter ("") or not. ("sb")
        """
        data_matrix = np.array(
            [
                [row["REF_DP"] - row["REF_RV"], row["REF_RV"]],
                [row["ALT_DP"] - row["ALT_RV"], row["ALT_RV"]],
            ]
        )
        _, pvalue = fisher_exact(data_matrix, alternative="greater")
        if pvalue < 0.05:
            return "sb"
        else:
            return ""

    def apply_filters(self, row):
        """Apply all the filters to a row and return its output merged

        Args:
            row (pd.Series)- a row from ivar.tsv dataframe as series

        Returns:
            vcf_filter (str): Results from each filter joined by ';'
        """
        ivar_filter = "" if row["PASS"] == True else "ft"
        try:
            float(row["ALT_QUAL"])
        except ValueError:
            bad_quality_filter = "bq"
        bad_quality_filter = "" if row["ALT_QUAL"] >= self.bad_qual_threshold else "bq"
        if not self.ignore_stbias:
            stb_filter = self.strand_bias_filter(row)
        else:
            stb_filter = ""
        all_filters = [ivar_filter, stb_filter, bad_quality_filter]
        vcf_filter = ";".join([x for x in all_filters if x])
        if not vcf_filter:
            vcf_filter = "PASS"
        return vcf_filter

    def initiate_vcf_df(self):
        """Read the input ivar.tsv file, process the data depending on the parameters
        selected when running the script and create a pandas dataframe with it

        Returns:
            vcf_df: Pandas dataframe after processing the file
        """
        ivar_df = self.raw_ivar_df.copy()
        ivar_df = ivar_df.dropna(thresh=2)
        if self.pass_only:
            ivar_df = ivar_df[ivar_df["PASS"] is True]
        ivar_df = ivar_df[ivar_df["ALT_FREQ"] >= self.freq_threshold]
        vcf_dict = {}
        vcf_dict["REGION"] = ivar_df["REGION"]
        vcf_dict["POS"] = ivar_df["POS"]
        vcf_dict["ID"] = ["."] * len(ivar_df)
        # Dealing with insertions and deletions
        vcf_dict["REF"] = np.where(
            ivar_df["ALT"].str[0] == "-",
            ivar_df["REF"] + ivar_df["ALT"].str[1:],
            ivar_df["REF"],
        )
        vcf_dict["ALT"] = np.where(
            ivar_df["ALT"].str[0] == "+",
            ivar_df["REF"] + ivar_df["ALT"].str[1:],
            np.where(ivar_df["ALT"].str[0] == "-", ivar_df["REF"], ivar_df["ALT"]),
        )
        vcf_dict["QUAL"] = ["."] * len(ivar_df)
        vcf_dict["FILTER"] = ivar_df.apply(self.apply_filters, axis=1)
        vcf_dict["INFO"] = np.select(
            [ivar_df["ALT"].str[0] == "+", ivar_df["ALT"].str[0] == "-"],
            ["TYPE=INS", "TYPE=DEL"],
            default="TYPE=SNP",
        )
        format_cols = [
            "GT",
            "DP",
            "REF_DP",
            "REF_RV",
            "REF_QUAL",
            "ALT_DP",
            "ALT_RV",
            "ALT_QUAL",
            "ALT_FREQ",
        ]
        vcf_dict["FORMAT"] = ":".join(format_cols)
        # simple workaround by setting the genotype to a constant 1
        ivar_df["DP"] = ivar_df["TOTAL_DP"]
        ivar_df["GT"] = 1
        vcf_dict["FILENAME"] = ivar_df[format_cols].astype(str).apply(":".join, axis=1)
        if not self.ignore_merge:
            # These columns are needed to merge codons. They will be deleted later
            vcf_dict["REF_CODON"] = ivar_df["REF_CODON"]
            vcf_dict["ALT_CODON"] = ivar_df["ALT_CODON"]
        vcf_df = pd.DataFrame.from_dict(vcf_dict)
        return vcf_df

    def find_consecutive(self, vcf_df):
        """Find and extract the consecutive variants in the vcf dataframe

        Args:
            vcf_df (pd.DataFrame): Pandas df from initiate_vcf_df()

        Returns:
            consecutive_df(pd.DataFrame): dataframe only with variants in
            consecutive positions, including those with duplicates
        """
        consecutive_mask = vcf_df["POS"].diff() <= 1
        consecutive_mask = consecutive_mask | consecutive_mask.shift(-1)
        consecutive_df = vcf_df[consecutive_mask]
        return consecutive_df

    def split_non_consecutive(self, consecutive_df):
        """Split rows that are not consecutive to form groups of consecutive variants

        Args:
            consecutive_df (pd.DataFrame): Consecutive variants from find_consecutive()

        Returns:
            split_rows_list (list): List of dataframes with consecutive rows
        """
        # Numpy raises a warning from a pandas function that is wrongly set as deprecated
        with warnings.catch_warnings():
            warnings.simplefilter(action="ignore", category=FutureWarning)
            split_rows_list = np.split(
                consecutive_df, np.where(np.diff(consecutive_df["POS"]) > 1)[0] + 1
            )
        return split_rows_list

    def find_duplicates(self, row_set):
        """Check wether there are duplicated variants in a group of consecutive rows

        Args:
            row_set(pd.DataFrame): Rows from split_non_consecutive() & get_same_codon()

        Returns:
            Bool: True if there are duplicates in the dataframe, False if not
        """
        consdups_mask = row_set["POS"].diff() == 0
        consdups_mask = consdups_mask | consdups_mask.shift(-1)
        if row_set[consdups_mask].empty:
            return False
        else:
            return True

    def get_same_codon(self, splitted_groups):
        """Receive a list of dataframes and exclude dataframes with only consecutive
        rows that share codon sequence, excluding those purely formed by duplicates.

        Args:
            splitted_groups (list(pd.DataFrame)): Groups from split_non_consecutive()

        Returns:
            same_codon_rows(list(pd.DataFrame)): list with separated dataframes
        """
        same_codon_dfs = []
        for consec_df in splitted_groups:
            clean_codon_rows = [
                self.find_consecutive(group)
                for _, group in consec_df.groupby("REF_CODON")
                if not self.find_consecutive(group).empty
                and not self.find_consecutive(group)["POS"].nunique() == 1
            ]
            same_codon_dfs.append(clean_codon_rows)
        same_codon_consecutive = [df for group in same_codon_dfs for df in group]
        return same_codon_consecutive

    def split_by_codon(self, same_codon_rows):
        """Split each dataframe into a dictionary with rows belonging to the same codon

        Args:
            same_codon_rows (list(pd.DataFrame)): Groups from get_same_codon()

        Returns:
            split_rows_dict(dict(index:pd.DataFrame)): Dictionary containing dataframe
            with rows belonging to the same codon as values and the index of the first
            row as keys, which is the position where the rows are going to be merged.
        """
        split_rows_dict = {}
        for rowset in same_codon_rows:
            last_pos = 0
            rows_groups = []
            first_index = None
            for row in rowset.itertuples():
                alt_pos = [x for x in range(3) if row.REF_CODON[x] != row.ALT_CODON[x]]
                if len(alt_pos) > 1:
                    print("Conflicting variants in position %s. Skipped" % row.POS)
                    continue
                alt_pos = alt_pos[0]
                first_index = row.Index if first_index is None else first_index
                if alt_pos < last_pos:
                    split_rows_dict[first_index] = pd.DataFrame(rows_groups)
                    rows_groups = []
                    first_index = row.Index
                rows_groups.append(row)
                last_pos = alt_pos
            split_rows_dict[first_index] = pd.DataFrame(rows_groups)
        return split_rows_dict

    def merge_rows(self, consec_rows):
        """Merge certain columns from a set of rows into the position of the first one

        Args:
            consec_rows (pd.DataFrame): Clean dataframe only with rows to be merged

        Returns:
            row_to_merge (list): List with the values for each cell in the merged row
        """
        merged_index = consec_rows.head(1).index[0]
        consec_rows.at[merged_index, "REF"] = "".join(consec_rows["REF"])
        consec_rows.at[merged_index, "ALT"] = "".join(consec_rows["ALT"])
        freqs_to_merge = ",".join(consec_rows["FILENAME"].str.split(":").str[8])
        depths_to_merge = ",".join(consec_rows["FILENAME"].str.split(":").str[5])
        stats_to_merge = ":".join([freqs_to_merge, depths_to_merge])
        consec_rows.at[merged_index, "FILENAME"] += f":{stats_to_merge}"
        consec_rows.at[merged_index, "FORMAT"] += ":MERGED_AF:MERGED_DP"
        lowest_af = min(freqs_to_merge.split(","))
        filecol_list = consec_rows.at[merged_index, "FILENAME"].split(":")
        filecol_list[8] = lowest_af
        consec_rows.at[merged_index, "FILENAME"] = ":".join(filecol_list)
        merged_row = consec_rows.loc[merged_index].values.tolist()
        return merged_row

    def merge_hap_vcf(self, clean_rows_list):
        """Merge all the given rows in a single one for each dataframe of consecutive
        rows in a given list

        Args:
            clean_rows_list (list(pd.DataFrame())): Dataframes with consecutive rows

        Returns:
            rows_to_merge (list(list())): List of merged rows as a list of values
        """
        rows_to_merge = []
        for rowbatch in clean_rows_list:
            if rowbatch.empty:
                continue
            if len(rowbatch) == 1:
                rows_to_merge.append(rowbatch.values.tolist()[0])
            else:
                if self.find_consecutive(rowbatch).empty:
                    continue
                merged_row = self.merge_rows(rowbatch)
                rows_to_merge.append(merged_row)
        return rows_to_merge

    def handle_dup_rows(self, row_set):
        """Split dataframe with multiple variants in the same position and create a list
        of rows for each possible resulting codon in the position of the first variant

        Args:
            row_set (pd.DataFrame): Consecutive variants df from split_rows_dict()

        Returns:
            merged_rowlist (list(list)): List of merged rows from merge_rows()
        """
        # Numpy raises a warning from a pandas function that is wrongly set as deprecated
        # https://github.com/numpy/numpy/issues/24889
        # https://github.com/apache/arrow/issues/36412
        with warnings.catch_warnings():
            warnings.simplefilter(action="ignore", category=FutureWarning)
            split_dups = np.split(
                row_set, np.where(np.diff(row_set["POS"]) >= 1)[0] + 1
            )
        split_indexes = [rowlist.index.to_list() for rowlist in split_dups]
        index_product_list = [indexlist for indexlist in product(*split_indexes)]
        merged_rowlist = []
        for indexlist in index_product_list:
            consec_rows = row_set.loc[indexlist, :]
            clean_rows_list = self.get_haplotypes(consec_rows.copy())
            cleaned_ref_rows_list = self.remove_edge_ref(clean_rows_list)
            batch_rowlist = self.merge_hap_vcf(cleaned_ref_rows_list.copy())
            merged_rowlist.extend(batch_rowlist)
        return merged_rowlist

    def get_ref_rowset(self, row_set):
        """Create a new row for each variant row in the dataframe emulating the
        reference for that position

        Args:
            row_set (pd.DataFrame()): A certain group of consecutive variants

        Returns:
            merged_ref_rows: Same Df but with duplicated rows that emulate reference
            for each variant position.
        """
        ref_row_set = row_set.copy()
        ref_row_set["ALT"] = ref_row_set["REF"]
        ref_row_set["ALT_CODON"] = ref_row_set["REF_CODON"]
        filecol = ref_row_set["FILENAME"].values.tolist()
        ref_filecol = []
        for row in filecol:
            # values = GT:DP:REF_DP:REF_RV:REF_QUAL:ALT_DP:ALT_RV:ALT_QUAL:ALT_FREQ
            split_vals = row.split(":")
            ref_dp = split_vals[2]
            total_dp = split_vals[1]
            split_vals[8] = str(round(int(ref_dp) / int(total_dp), 3))
            split_vals[5] = ref_dp
            ref_filecol.append(":".join(split_vals))
        ref_row_set["FILENAME"] = ref_filecol
        merged_ref_rows = (
            pd.concat([row_set, ref_row_set]).sort_values("POS").reset_index(drop=True)
        )
        return merged_ref_rows

    def merge_rule_check(self, alt_dictlist):
        """Evaluate a list of possible codons and decide if they can be merged together
        based on certain conditions regarding Allele Frequency

        Args:
            alt_dictlist (list(dict()): A list of dictionaries with consecutive locus

        Returns:
            consec_series (list(dict()): Same list but only with validated locus
        """
        consec_series = []

        for altdict in alt_dictlist:
            if len(altdict) <= 1:
                consec_series.append(altdict)
                continue
            af_list = [x["AF"] for x in altdict.values()]
            key_list = list(altdict.keys())
            distances = []
            for i in range(len(altdict) - 1):
                key1 = key_list[i]
                key2 = key_list[i + 1]
                distance = abs(altdict[key2]["AF"] - altdict[key1]["AF"])
                distances.append(distance)

            # First combination is always alt,alt,alt
            # If all af > consensus_af keep just this combination and return
            if all(af > self.consensus_af for af in af_list):
                consec_series = [altdict]
                return consec_series

            # If all af between 04 and 0.6 keep just this combination and return.
            # Note: All other combinations are going to be valid considering distance
            # but we've empirically checked that only ALT combination is found in real bams.
            if all(af > 0.4 or af < 0.6 for af in af_list):
                consec_series = [altdict]
                return consec_series

            # If any af == 0 combination not posible. Continue
            if any(af == 0 for af in af_list):
                continue

            # If all distances < distance threshold, haplotype is valid
            if all(dist < self.merge_af_threshold for dist in distances):
                consec_series.append(altdict)
                continue

            # If any distance > distance threshold, check if needed nucleotide af is > af_consensus
            correct = False
            single_add = {}
            for idx, dist in enumerate(distances):
                if dist > self.merge_af_threshold:
                    for i in range(idx, idx + 1):
                        if af_list[i] > self.consensus_af:
                            single_add[0] = altdict[i]
                            if single_add not in consec_series:
                                consec_series.append(single_add)
                            correct = True
                        else:
                            correct = False
                            break
                else:
                    correct = True

                if not correct:
                    break

            if correct:
                consec_series.append(altdict)
        return consec_series

    def get_haplotypes(self, consec_rows):
        """Create a list of all possible combinations of REF and ALT consecutive codons
        following certain conditions of similarity using Allele Frequency values

        Args:
            consec_rows (list(pd.DataFrame)): List of dataframes with consecutive rows

        Returns:
            clean_rows_list (list(pd.DataFrame)): Filtered list with viable combinations
        """
        # Create all posible combination between ref and alt (different haplotypes) for consecutive variants
        merged_ref_rows = self.get_ref_rowset(consec_rows.copy())
        merged_ref_rows["AF"] = merged_ref_rows["FILENAME"].str.split(":").str[8]
        alt_rows = merged_ref_rows[
            merged_ref_rows["REF_CODON"] != merged_ref_rows["ALT_CODON"]
        ].reset_index(drop=True)
        ref_rows = merged_ref_rows[
            merged_ref_rows["REF_CODON"] == merged_ref_rows["ALT_CODON"]
        ].reset_index(drop=True)
        ref_dict = {
            x: {"AF": float(y), "set": "ref"}
            for x, y in ref_rows["AF"].to_dict().items()
        }
        alt_dict = {
            x: {"AF": float(y), "set": "alt"}
            for x, y in alt_rows["AF"].to_dict().items()
        }
        hap_combinations = list(
            product(*[[(k, alt_dict[k]), (k, ref_dict[k])] for k in alt_dict.keys()])
        )
        hap_comb_dict = [dict(comb) for comb in hap_combinations]

        # Keep together only those haplotypes that fulfill certain similarity rules
        valid_hap = self.merge_rule_check(hap_comb_dict)

        if self.consensus_merge:
            valid_hap_fil = []

            # Process each dictionary in the list
            for hap in valid_hap:
                # Create a dictionary of items that meet the threshold
                valid = all(value['AF'] > self.consensus_af for value in hap.values())
                if valid:
                    valid_hap_fil.append(hap)
                else:
                    # Create a list of dictionaries for items that do not meet the threshold
                    separate_dicts = [{key: value} for key, value in hap.items() if value['AF'] <= self.consensus_af]

                    # Extend filtered_list with separate_dicts
                    for sep in separate_dicts:
                        if sep not in valid_hap_fil:
                            valid_hap_fil.append(sep)
        else:
            valid_hap_fil = valid_hap

        clean_rows_list = []
        for hap in valid_hap_fil:
            consec_rowsdict = {}
            for key, vals in hap.items():
                if vals["set"] == "alt":
                    consec_rowsdict[key] = alt_rows.loc[int(key)]
                else:
                    consec_rowsdict[key] = ref_rows.loc[int(key)]
            consec_df = pd.DataFrame.from_dict(consec_rowsdict, orient="index")

            clean_loc_df = consec_df.drop("AF", axis=1)
            if not clean_loc_df.empty:
                clean_rows_list.append(clean_loc_df)
        return clean_rows_list

    def remove_edge_ref(self, clean_rows_list):
        """Remove reference nucleotides from both edges of the variant codon

        Args:
            clean_rows_list (List(pd.DataFrame)): List of variants to be merged
        Returns:
            cleaned_ref_rows_list (List(pd.DataFrame)): List of rows without edge refs
        """

        cleaned_ref_rows_list = []
        for df in clean_rows_list:
            ref_col = df["REF"]
            alt_col = df["ALT"]
            idx_matches = [idx for idx in alt_col.index if alt_col[idx] == ref_col[idx]]
            if not idx_matches:
                cleaned_ref_rows_list.append(df)
                continue
            if len(df) == 3 and idx_matches == [1]:
                cleaned_ref_rows_list.append(df)
                continue

            df = df.drop(idx_matches, axis=0)
            if not df.empty:
                cleaned_ref_rows_list.append(df)

        return cleaned_ref_rows_list

    def process_vcf_df(self, vcf_df):
        """Merge rows with consecutive SNPs that passed all filters and without NAs

        Args:
            vcf_df - dataframe from self.initiate_vcf_df()
        Returns:
            processed_vcf_df: dataframe with consecutive variants merged
        """

        def include_rows(vcf_df, last_index, rows_to_merge):
            indexes_to_merge = [
                x for x in range(last_index, last_index + len(rows_to_merge))
            ]
            for index, row in zip(indexes_to_merge, rows_to_merge):
                try:
                    vcf_df.loc[index] = row
                except ValueError:
                    print(f"Invalid row found: {str(row)}. Skipped")
            return vcf_df

        # Clean ivar tsv duplicates due to duplicated annotation
        vcf_df = vcf_df.drop_duplicates(subset=["REGION", "POS", "REF", "ALT"])
        vcf_df_fil = vcf_df[vcf_df["INFO"] == "TYPE=SNP"]
        vcf_df_fil = vcf_df_fil[vcf_df_fil["FILTER"] == "PASS"].dropna()
        consecutive_df = self.find_consecutive(vcf_df_fil)
        if consecutive_df.empty:
            return vcf_df
        consecutive_split = self.split_non_consecutive(consecutive_df)
        codon_split = self.get_same_codon(consecutive_split)
        codon_split_dict = self.split_by_codon(codon_split)
        for _, row_set in sorted(codon_split_dict.items()):
            last_index = vcf_df.tail(1).index[0] + 1
            vcf_df = vcf_df.drop(row_set.Index)
            # Redundant "Index" column is generated by Itertuples in split_by_codon()
            row_set = row_set.drop(["Index"], axis=1)
            if self.find_duplicates(row_set):
                haps_merged = self.handle_dup_rows(row_set.copy())
            else:
                hap_rows = self.get_haplotypes(row_set.copy())
                hap_rows_noref = self.remove_edge_ref(hap_rows)
                haps_merged = self.merge_hap_vcf(hap_rows_noref.copy())
            vcf_df = include_rows(vcf_df, last_index, haps_merged)
        vcf_df = vcf_df[vcf_df["REF"] != vcf_df["ALT"]]
        vcf_df = vcf_df.sort_index().reset_index(drop=True)
        processed_vcf_df = vcf_df.drop_duplicates().sort_values("POS")
        return processed_vcf_df

    def get_vcf_header(self):
        """Create the vcf header for VCFv4.2

        Returns:
            header: String containing all the vcf header lines separated by newline.
        """
        # Define VCF header
        header_source = ["##fileformat=VCFv4.2", "##source=iVar"]
        if self.ref_fasta:
            header_contig = []
            for record in SeqIO.parse(self.ref_fasta, "fasta"):
                header_contig += [
                    "##contig=<ID="
                    + record.id
                    + ",length="
                    + str(len(record.seq))
                    + ">"
                ]

            header_source += header_contig

        header_info = [
            '##INFO=<ID=TYPE,Number=1,Type=String,Description="Either SNP (Single Nucleotide Polymorphism), DEL (deletion) or INS (Insertion)">',
        ]
        header_filter = [
            '##FILTER=<ID=PASS,Description="All filters passed">',
            '##FILTER=<ID=ft,Description="Fisher\'s exact test of variant frequency compared to mean error rate, p-value > 0.05">',
            '##FILTER=<ID=bq,Description="Bad quality variant: ALT_QUAL lower than 20">',
        ]
        if not self.ignore_stbias:
            header_filter.append(
                '##FILTER=<ID=sb,Description="Strand bias filter not passed">'
            )
        header_format = [
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
            '##FORMAT=<ID=REF_DP,Number=1,Type=Integer,Description="Depth of reference base">',
            '##FORMAT=<ID=REF_RV,Number=1,Type=Integer,Description="Depth of reference base on reverse reads">',
            '##FORMAT=<ID=REF_QUAL,Number=1,Type=Integer,Description="Mean quality of reference base">',
            '##FORMAT=<ID=ALT_DP,Number=1,Type=Integer,Description="Depth of alternate base">',
            '##FORMAT=<ID=ALT_RV,Number=1,Type=Integer,Description="Depth of alternate base on reverse reads">',
            '##FORMAT=<ID=ALT_QUAL,Number=1,Type=Integer,Description="Mean quality of alternate base">',
            '##FORMAT=<ID=ALT_FREQ,Number=1,Type=Float,Description="Frequency of alternate base">',
        ]
        if not self.ignore_merge:
            header_format.append(
                '##FORMAT=<ID=MERGED_AF,Number=A,Type=Float,Description="Frequency of each merged variant comma separated">'
            )
            header_format.append(
                '##FORMAT=<ID=MERGED_DP,Number=A,Type=Float,Description="Total Depth of each merged variant comma separated">'
            )
        header = header_source + header_info + header_filter + header_format
        return header

    def write_vcf(self):
        """Process ivar.tsv, merge the vcf header and table and write them into a file"""
        vcf_header = "\n".join(self.get_vcf_header())
        vcf_table = self.initiate_vcf_df()

        def export_vcf(vcf_table, consensus=True):
            self.consensus_merge = consensus
            if not self.ignore_merge:
                processed_vcf = self.process_vcf_df(vcf_table)
            else:
                processed_vcf = vcf_table

            try:
                processed_vcf = processed_vcf.drop(["REF_CODON", "ALT_CODON"], axis=1)
            except KeyError:
                pass
            # Workaround because itertuples cannot handle special characters in column names
            processed_vcf = processed_vcf.rename(
                columns={"REGION": "#CHROM", "FILENAME": self.filename}
            )
            if consensus:
                filepath = self.file_out
            else:
                basename = os.path.splitext(os.path.basename(self.file_out))[0]
                filename = str(basename) + "_all_hap.vcf"
                filepath = os.path.join(os.path.dirname(self.file_out), filename)
            with open(filepath, "w") as file_out:
                file_out.write(vcf_header + "\n")
            processed_vcf.to_csv(filepath, sep="\t", index=False, header=True, mode="a")

        export_vcf(vcf_table, consensus=True)
        export_vcf(vcf_table, consensus=False)
        return


def make_dir(path):
    """
    Description:
        Create directory if it doesn't exist.
    Args:
        path - path where the directory will be created.
    """
    if not len(path) == 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise
    return


def main(args=None):
    args = parse_args(args)

    ivar_to_vcf = IvarVariants(
        file_in=args.file_in,
        file_out=args.file_out,
        pass_only=args.pass_only,
        freq_threshold=args.allele_freq_threshold,
        bad_qual_threshold=args.bad_quality_threshold,
        merge_af_threshold=args.merge_af_threshold,
        consensus_af=args.consensus_af,
        ignore_stbias=args.ignore_strand_bias,
        ignore_merge=args.ignore_merge_codons,
        fasta=args.fasta,
    )
    ivar_to_vcf.write_vcf()
    return


if __name__ == "__main__":
    sys.exit(main())
