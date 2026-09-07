#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import json
import logging
import argparse
import re
import pandas as pd

logger = logging.getLogger()
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", None)


def parser_args(args=None):
    Description = "Parse Sierra-local JSON reports and corresponding codfreq file and generate tables with mutations and resistance information."
    Epilog = """Example usage:
    python resistance_tables.py --sierralocal_file sample_resistance.json --codfreq_file sample.codfreq --output_mutation_file sample_mutation_table.csv --output_resistance_file sample_resistance_table.csv
    """
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)

    parser.add_argument("-sf", "--sierralocal_file", type=str, required=True, help="JSON file containing sierra-local report.")
    parser.add_argument("-cf", "--codfreq_file", type=str, required=True, help="Path to codfreq file.")
    parser.add_argument("-s", "--sample_name", type=str, required=True, help="Name of the sample")
    parser.add_argument( "-om", "--output_mutation_file", required=True, type=str, help="Full path to output mutation CSV file.")
    parser.add_argument("-os", "--output_mutation_short", required=True, type=str,help="Full path to output mutation shortenned CSV file.")
    parser.add_argument("-or", "--output_resistance_file", required=True, type=str, help="Full path to output resistance CSV file.")
    parser.add_argument("-maf", "--min_allele_frequency", default=0.9, type=float, help="Minimum allele frequency threshold used in iVar consensus calling.")
    return parser.parse_args(args)

def parse_codfreq(codfreq_path):
    return pd.read_csv(codfreq_path)

def build_mutation_row(sample_name, gene_name, mut_text, original_mut_text, mut, resistance_comments):
    return {
        "Sample_name": sample_name,
        "Gene_name": gene_name,
        "Mutations": mut_text,
        "Mutations_type": mut.get("primaryType", "NA"),
        # Use the original combined mutation text, such as M184IM,
        # because sierra-local stores the comment under that mutation name.
        "Mutations_comments": resistance_comments.get(original_mut_text, ""),
        "isInsertion": mut.get("isInsertion", False),
        "isDeletion": mut.get("isDeletion", False),
        "isApobecMutation": mut.get("isApobecMutation", False),
        "isApobecDRM": mut.get("isApobecDRM", False),
        "isUnusual": mut.get("isUnusual", False),
        "isSDRM": mut.get("isSDRM", False),
        "hasStop": mut.get("hasStop", False),
        "Mutation_AF": "NA",
        "Coverage": "NA",
        "INDEL>5%": "NA",
    }

def parse_sierra_json(sample_name, json_path):
    """Parse one Sierra-local JSON and return a pandas DataFrame with all mutations."""
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)[0]

    rows = []

    # === 1. Extract comments of resistance associated with mutations ===
    resistance_comments = {}
    if "drugResistance" in data:
        for entry in data["drugResistance"]:
            for drugscore in entry.get("drugScores", []):
                for partial in drugscore.get("partialScores", []):
                    for mut in partial.get("mutations", []):
                        mut_name = mut.get("text", "")
                        comments = mut.get("comments", [])
                        if comments:
                            # Concatenate all comments if there is more than one
                            comment_text = " ".join([c.get("text", "") for c in comments if c.get("text")])
                            # Save only if it does not already exist or if it is longer (avoids overwriting with duplicates)
                            if mut_name not in resistance_comments or len(comment_text) > len(resistance_comments[mut_name]):
                                resistance_comments[mut_name] = comment_text

    # === 2. Detect warnings ===
    warnings_dict = {}
    for warning in data.get("validationResults", []):
        msg = warning.get("message", "")
        if "sequence had" in msg:
            # Replace "3\u2032-end" with "3'-end"
            msg = msg.replace("3\u2032-end", "3'-end")
            # Detect affected protein (RT, PR, or IN)
            for gene_name in ["RT", "PR", "IN"]:
                if f"('{gene_name}'," in msg:
                    warnings_dict[gene_name] = msg

    if "alignedGeneSequences" not in data:
        logger.warning(f"No 'alignedGeneSequences' field found in {json_path}")
        return pd.DataFrame()

    # === 3. Extract mutations ===
    for gene_entry in data["alignedGeneSequences"]:
        gene_name = gene_entry.get("gene", {}).get("name", "NA")
        lastAA = gene_entry['lastAA']

        for mut in gene_entry.get("mutations", []):
            consensus = mut.get("consensus", "")
            text = mut.get("text", "")
            # Keep the original combined mutation text so split mutations
            # can retrieve the resistance comment linked to the combined form. For instance: M184IM
            original_mut_text = text
            aas = mut.get("AAs", "")
            pos = mut.get("position", "")
            match = re.match(rf"{re.escape(consensus)}{str(pos)}(.+)", text)
            if match:
                mutant = match.group(1) # "KR" for "K70KR"
            else:
                mutant = ""  # fallback
            # If there are more than one possible amino acids, create a row for each.
            # Deletions are reported as e.g. P90del and should remain a single row.
            if len(mutant) > 1 and not mut.get("isDeletion"):
                for aa in mutant:
                    if pos == lastAA and aa == "X" and not mut.get("isDeletion"):
                        continue  # skip, false positive at end of sequence
                    mut_text = f"{consensus}{pos}{aa}"
                    rows.append(build_mutation_row(sample_name, gene_name, mut_text, original_mut_text, mut, resistance_comments))
            else:
                if pos == lastAA and aas == "X" and not mut.get("isDeletion"):
                    continue  # skip, false positive at end of sequence
                # If there is only one possible amino acid, keep a single row
                mut_text = mut.get("text", "NA")
                rows.append(build_mutation_row(sample_name, gene_name, mut_text, original_mut_text, mut, resistance_comments))

    df = pd.DataFrame(rows)

    return df

def is_insertion_codon(codon):
    return isinstance(codon, str) and len(codon) > 3

def is_deletion_codon(codon):
    return isinstance(codon, str) and "-" in codon

def extract_mutation_position(mutation):
    match = re.match(r"^[A-Z*_-](\d+)", str(mutation))
    if not match:
        raise ValueError(f"Could not extract mutation position from {mutation}")
    return match.group(1)

def integrate_codfreq_info(df_json, codfreq_df):
    """Update df_json with Mutation_AF and Coverage from codfreq_df."""
    if codfreq_df.empty:
        return df_json

    updated_rows = []

    for _, row in df_json.iterrows():
        gene = row["Gene_name"]
        pos = extract_mutation_position(row["Mutations"])
        aa = row["Mutations"][-1]  # Extract mutated amino acid
        is_indel = row["isInsertion"] or row["isDeletion"]

        # Search that position in codfreq
        candidates = codfreq_df[
            (codfreq_df["gene"] == gene) &
            (codfreq_df["position"] == int(pos))
        ]
        if candidates.empty:
            continue

        total = candidates["total"].iloc[0]

        if is_indel:
            if aa == "_" and row["isInsertion"] and not candidates.empty:
                # Use the insertion with the highest count
                insertions = candidates[candidates["codon"].apply(is_insertion_codon)]
                coverage = insertions["count"].max() if not insertions.empty else 0
            elif row["isInsertion"] and not candidates.empty:
                # Normal codon close to insertion -> look for the codon corresponding to the AA
                subset = candidates[candidates["aa_codon"] == aa]
                coverage = subset["count"].sum() if not subset.empty else 0
                row["isInsertion"] = False  # Deactivate isInsertion
            elif row["isDeletion"] and not candidates.empty:
                # Use the deletion codon with the highest count at this position.
                # TODO: Validate this scenario on more deletion positive samples.
                deletions = candidates[candidates["codon"].apply(is_deletion_codon)]
                coverage = deletions["count"].max() if not deletions.empty else 0
            else:
                raise RuntimeError(f"NEW CASE SCENARIO not yet implemented for gene={gene}, pos={pos}, aa={aa}")
        else:
            # Normal codon
            subset = candidates[candidates["aa_codon"] == aa]
            coverage = subset["count"].sum() if not subset.empty else 0

        af = coverage / total if total > 0 else 0
        row["Mutation_AF"] = round(af, 5)
        row["Coverage"] = int(total)

        # Add "INDEL>5%" column
        indel_sum = candidates[candidates["codon"].apply(lambda x: is_insertion_codon(x) or is_deletion_codon(x))]["count"].sum()
        row["INDEL>5%"] = indel_sum / total > 0.05

        updated_rows.append(row)

    final_df = pd.DataFrame(updated_rows)

    return final_df

def parse_resistance_json(sample_name, json_path):
    """Parse Sierra-local JSON and return a pandas DataFrame with drug resistance information."""
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)[0]

    # Dictionary to map abbreviations to full drug names
    drug_fullnames = {
        # PI:
        "ATV/r": "atazanavir/r",
        "DRV/r": "darunavir/r",
        "LPV/r": "lopinavir/r",
        "FPV/r": "fosamprenavir/r",
        "IDV/r": "indinavir/r",
        "NFV": "nelfinavir",
        "SQV/r": "saquinavir/r",
        "TPV/r": "tipranavir/r",
        # NRTI:
        "ABC": "abacavir",
        "AZT": "zidovudine",
        "FTC": "emtricitabine",
        "ISL": "islatravir",
        "3TC": "lamivudine",
        "TDF": "tenofovir",
        "D4T": "stavudine",
        "DDI": "didanosine",
        # NNRTI:
        "DOR": "doravirine",
        "EFV": "efavirenz",
        "ETR": "etravirine",
        "NVP": "nevirapine",
        "RPV": "rilpivirine",
        "DPV": "dapivirine",
        # INSTI:
        "BIC": "bictegravir",
        "CAB": "cabotegravir",
        "DTG": "dolutegravir",
        "EVG": "elvitegravir",
        "RAL": "raltegravir",
    }

    if "drugResistance" not in data:
        logger.warning(f"No 'drugResistance' field found in {json_path}")
        return pd.DataFrame()

    rows = []
    for entry in data["drugResistance"]:
        gene_name = entry.get("gene", {}).get("name", "NA")

        for drugscore in entry.get("drugScores", []):
            drug_class = drugscore.get("drugClass", {}).get("name", "NA")
            drug_abbr = drugscore.get("drug", {}).get("displayAbbr", "NA")
            total_score = drugscore.get("score", "NA")
            res_status = drugscore.get("text", "NA")

            # Map to full drug name (NA if abbreviation not found)
            drug_name = drug_fullnames.get(drug_abbr, "NA")

            row = {
                "Sample_name": sample_name,
                "Gene_name": gene_name,
                "Drug_class": drug_class,
                "Drug_abbr": drug_abbr,
                "Drug_name": drug_name,
                "Total_score": total_score,
                "Res_status": res_status,
            }
            rows.append(row)

    df = pd.DataFrame(rows)
    return df

def filter_mutations(mutation_df, min_allele_frequency):
    """
    Filter mutations with too low allele frequency using iVar consensus-style cumulative frequency logic.

    Sierra-local may report more than two amino acid variants for the same codon
    when more than one ambiguous nucleotide site is present in the consensus sequence.
    This keeps amino acid variants ordered by allele frequency until
    the cumulative allele frequency reaches the minimum allele frequency threshold.
    """
    if mutation_df.empty:
        return mutation_df

    mutation_df = mutation_df.copy()
    mutation_df["_original_order"] = range(len(mutation_df))

    mutation_df = mutation_df[~((mutation_df["Mutations"].str[-1] == "X") & (mutation_df["Mutation_AF"] == 0))]

    mutation_df["_position"] = mutation_df["Mutations"].apply(extract_mutation_position).astype(int)

    mutation_df = mutation_df.sort_values(
        ["Gene_name", "_position", "Mutation_AF"],
        ascending=[True, True, False],
    )

    mutation_df["_cumulative_AF"] = mutation_df.groupby(["Gene_name", "_position"])["Mutation_AF"].cumsum()
    mutation_df["_cumulative_AF_before"] = mutation_df["_cumulative_AF"] - mutation_df["Mutation_AF"]

    mutation_df = mutation_df[mutation_df["_cumulative_AF_before"] < min_allele_frequency]
    mutation_df = mutation_df.sort_values("_original_order")
    return mutation_df.drop(columns=["_original_order", "_position", "_cumulative_AF", "_cumulative_AF_before"])

def main(args=None):
    args = parser_args(args)

    # Load sierra-local JSON files
    sierralocal_df = parse_sierra_json(args.sample_name, args.sierralocal_file)

    # Load codfreq files
    codfreq_df = parse_codfreq(args.codfreq_file)

    # Integrate codfreq values
    mutation_df = integrate_codfreq_info(sierralocal_df, codfreq_df)

    # Filter the DataFrame to remove those rows
    mutation_df = filter_mutations(mutation_df, args.min_allele_frequency)

    if mutation_df.empty:
        logger.warning(f"No mutations found for sample {args.sample_name} or no valid codfreq data.")

    mutation_df.to_csv(args.output_mutation_file, index=False, encoding="utf-8-sig")
    print(f"✅ Resistance table saved to {args.output_mutation_file}")

    # === Generate mutation table with fewer fields===
    filtered_mutation_df = mutation_df[["Sample_name", "Gene_name", "Mutations", "Mutations_type", "Mutations_comments"]]

    filtered_mutation_df.to_csv(args.output_mutation_short, index=False, encoding="utf-8-sig")
    print(f"✅ Resistance table saved to {args.output_mutation_short}")

    # Parse drug resistance info
    resistance_df = parse_resistance_json(args.sample_name, args.sierralocal_file)

    if resistance_df.empty:
        logger.warning(f"No resistance information found for sample {args.sample_name}")
    else:
        resistance_df.to_csv(args.output_resistance_file, index=False, encoding="utf-8-sig")
        print(f"✅ Resistance table saved to {args.output_resistance_file}")

if __name__ == "__main__":
    sys.exit(main())
