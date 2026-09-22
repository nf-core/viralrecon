#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import json
import argparse
import base64
import pandas as pd
from datetime import date
from Bio import SeqIO
from jinja2 import Environment, FileSystemLoader, select_autoescape


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

def parser_args(args=None):
    Description = "Parse Sierra-local JSON reports and corresponding resistance and mutation tables to generate an HTML report per sample."
    Epilog = """Example usage:
    python resistance_report.py --sierralocal_json SAMPLE_resistance.json --mutation_csv SAMPLE_mutation_table.csv
        --resistance_csv SAMPLE_resistance_table.csv --nextclade_csv SAMPLE_nextclade.csv
        --consensus_fasta SAMPLE.fa --ivar_consensus_params "-t 0.8 -q 30 -m 50 -n N"
        --ivar_variant_maf 0.01 --output_html SAMPLE_resistance_report.html
    """
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)

    parser.add_argument(
        "-s",
        "--sierralocal_json",
        type=str,
        required=True,
        help="Sierra-local JSON report.",
    )
    parser.add_argument(
        "-m",
        "--mutation_csv",
        type=str,
        required=True,
        help="Mutation CSV file.",
    )
    parser.add_argument(
        "-r",
        "--resistance_csv",
        type=str,
        required=True,
        help="Resistance CSV file for one sample.",
    )
    parser.add_argument(
        "-n",
        "--nextclade_csv",
        type=str,
        required=True,
        help="Nextclade CSV file.",
    )
    parser.add_argument(
        "-cn",
        "--consensus_fasta",
        type=str,
        required=True,
        help="Consensus FASTA file.",
    )
    parser.add_argument(
        "-ic",
        "--ivar_consensus_params",
        type=str,
        help="Parameters used for ivar consensus calling",
    )
    parser.add_argument(
        "-iv",
        "--ivar_variant_maf",
        type=float,
        default=0.01,
        help="Minor allele frequency threshold used for ivar variant calling",
    )
    parser.add_argument(
        "-d",
        "--deprecated_drugs",
        type=str,
        default="",
        help="Comma-separated list of deprecated drugs that should be removed from the final report (for example: D4T,DDI,DPV,FPV/r,IDV/r,NFV,SQV/r,TPV/r)",
    )
    parser.add_argument(
        "-o",
        "--output_html",
        type=str,
        required=True,
        help="Full path to the sample HTML report file."
    )
    parser.add_argument(
        "--nextclade_dataset_name",
        type=str,
        default="",
        help="Nextclade dataset name, for example neherlab/hiv-1",
    )
    parser.add_argument(
        "--nextclade_dataset_tag",
        type=str,
        default="",
        help="Nextclade dataset tag, for example 2025-09-09--12-13-13Z",
    )
    parser.add_argument(
        "--pipeline_version",
        type=str,
        default="dev",
        help="nf-core/viralrecon pipeline version",
    )

    return parser.parse_args(args)

def read_consensus_records(consensus_fasta):
    records = {}

    for record in SeqIO.parse(consensus_fasta, "fasta"):
        _, separator, region = record.id.rpartition("|")
        if not separator or not region:
            raise ValueError(
                f"Consensus record ID must end with '|REGION': {record.id}"
            )
        records[region] = record

    return records

def estimate_ambiguous_site_proportion(consensus_records):
    """
    Estimate the proportion of ambiguous sites across POL gene.

    Ambiguous sites are IUPAC mixed-base codes (R, Y, S, W, K, M, B, D, H, and V).
    The denominator is the number of non-N sites in both PR_RT and IN.
    """

    ambiguous_codes = set("RYSWKMBDHV")
    regions = ("PR_RT", "IN")

    ambiguous_sites = 0
    non_n_sites = 0

    for region in regions:
        sequence = consensus_records[region].seq.upper()
        ambiguous_sites += sum(base in ambiguous_codes for base in sequence)
        non_n_sites += sum(base != "N" for base in sequence)

    proportion = ambiguous_sites / non_n_sites if non_n_sites else 0.0

    return {
        "ambiguous_sites": ambiguous_sites,
        "non_n_consensus_sites": non_n_sites,
        "proportion": proportion,
    }

def get_nextclade_subtypes(nextclade_file, consensus_records):
    df_next = pd.read_csv(nextclade_file, sep=";")
    clades = df_next.set_index("seqName")["clade"]

    return [
        {
            "region": region,
            "label": region,
            "subtype": clades.loc[record.id],
        }
        for region, record in consensus_records.items()
        if region in ("PR_RT", "IN")
    ]

def parse_sequence_summary(json_path,
                           ivar_consensus_params=None, ivar_variant_maf=None,
                           ambiguous_site_summary=None):
    """
    Extract sequence summary information from Sierra-local JSON.
    - Lists each gene present (PR, RT, IN)
    - Detects missing nucleotide ranges from warnings
    - Detects missing codons at the start if firstAA != 1
    """

    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)[0]

    summary_lines = []
    warnings = {}
    not_sequenced = {}

    # --- Define expected protein lengths
    protein_lengths = {
        "PR": 99,
        "RT": 560,
        "IN": 288
    }

    # --- Parse validation warnings to find missing nucleotide ranges
    for warning in data.get("validationResults", []):
        msg = warning.get("message", "")
        match = re.search(r"\('([A-Z]{2})'.*?,\s*(\d+),\s*(\d+)\)", msg)
        if match:
            gene = match.group(1)
            start_nt = match.group(2)
            end_nt = match.group(3)
            warnings[gene] = (start_nt, end_nt)

    # --- Parse alignedGeneSequences
    for gene_entry in data.get("alignedGeneSequences", []):
        gene = gene_entry.get("gene", {}).get("name", "NA")
        first_aa = gene_entry.get("firstAA", "NA")
        last_aa = gene_entry.get("lastAA", "NA")
        not_sequenced.setdefault(gene, [])

        # --- Detect not sequenced positions from mutations ending with X
        for mut in gene_entry.get("mutations", []):
            if mut.get("text", "").endswith("X") and mut.get("AAs", "") in ["*ACDEFGHIKLMNPQRSTVWY", "ACDFGHILNPRSTVY"]:
                pos = mut.get("position")
                if pos and isinstance(pos, int):
                    not_sequenced[gene].append(pos)

        # --- Build list of contiguous missing intervals
        missing_positions = sorted(not_sequenced.get(gene, []))
        missing_parts = []

        # --- Adjust start
        adjusted_first = first_aa
        for pos in missing_positions:
            if pos == adjusted_first:
                adjusted_first += 1
            else:
                break

        # --- Adjust end
        adjusted_last = last_aa
        true_end = protein_lengths.get(gene, last_aa)
        for pos in reversed(missing_positions):
            if pos == true_end-1:
                missing_positions.append(true_end)
                adjusted_last = pos - 1
            if pos >= adjusted_last:
                adjusted_last = pos - 1
            else:
                break

        # --- Build missing_parts text
        if missing_positions:
            start = end = missing_positions[0]
            blocks = []
            for pos in missing_positions[1:]:
                if pos == end + 1:
                    end = pos
                else:
                    blocks.append((start, end))
                    start = end = pos
            blocks.append((start, end))
            for s, e in blocks:
                missing_parts.append(f"{s}" if s == e else f"{s}-{e}")

        # Add missing at tail if last_aa < expected
        if adjusted_last < true_end:
            tail_range = f"{adjusted_last+1}-{true_end}"
            if tail_range not in missing_parts:
                missing_parts.append(tail_range)

        # Build output line
        line = f"Sequence includes {gene}: codons {adjusted_first} - {adjusted_last}"
        if missing_parts:
            missing_text = ", ".join(sorted(missing_parts, key=lambda x: int(x.split('-')[0])))
            line += f" (missing: {missing_text})"

        summary_lines.append(line)

    # --- Add proportion of ambiguous consensus sites
    summary_lines.append(
        f"% pol polymorphisms (ambiguous sites): "
        f"{ambiguous_site_summary['proportion'] * 100:.2f}%"
    )

    # --- Parse ivar consensus parameters if provided
    # Extract numeric values with regex
    match_t = re.search(r"-t\s*([\d.]+)", ivar_consensus_params)
    match_q = re.search(r"-q\s*(\d+)", ivar_consensus_params)
    match_m = re.search(r"-m\s*(\d+)", ivar_consensus_params)

    t_val = float(match_t.group(1))
    q_val = int(match_q.group(1))
    m_val = int(match_m.group(1))

    summary_lines.append(f"Minimum read depth: ≥{m_val}")
    summary_lines.append(f"Nucleotide mixture threshold (NMT): ≥{ ivar_variant_maf * 100:.0f}%")
    summary_lines.append(f"Mutation detection threshold (MDT): ≥{ (1-t_val) * 100:.0f}%")
    summary_lines.append(f"Minimum quality threshold: {q_val}")

    return summary_lines

def extract_triggered_mutations (mutation):
    """Return the mutation lable using only amino acids that triggered sierra scoring"""
    text = mutation.get("text", "")
    triggered_aas = mutation.get("triggeredAAs", "")

    match = re.match(r"^([A-Z*_-]\d+)", text)
    if not match or not triggered_aas:
        return text

    mutation_prefix = match.group(1)
    return f"{mutation_prefix}{triggered_aas}"

def extract_mutation_scoring(json_path):
    mutation_scores = {}
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)[0]
    for protein_resistance in data.get("drugResistance", []):
        gene = protein_resistance.get("gene", {}).get("name")
        if not gene:
            continue
        if gene not in mutation_scores:
            mutation_scores[gene] = {}
        for drug in protein_resistance.get("drugScores", []):
            drug_name = drug.get("drug", {}).get("displayAbbr", "")
            drug_class = drug.get("drugClass", {}).get("name", "")
            for mutation_block in drug.get("partialScores", []):
                mutations = mutation_block.get("mutations", [])
                score = mutation_block.get("score")

                # Extract all mutation names
                mutation_ids = [extract_triggered_mutations(m) for m in mutations]

                # CASE 1 → single mutation: "M41L"
                if len(mutation_ids) == 1:
                    mutation_key = mutation_ids[0]

                # CASE 2 → combo mutation: "M41L+L210W"
                else:
                    mutation_key = "+".join(mutation_ids)

                # Ensure gene and mutation entry exist
                if mutation_key not in mutation_scores[gene]:
                    mutation_scores[gene][mutation_key] = {}
                mutation_scores[gene][mutation_key][drug_name] = {
                    "score": score,
                    "drug_class": drug_class
                }
    return mutation_scores

def parse_mutation_key(mutation_key):
    """
    Convert mutation key like 'M41L+T215Y' into a sortable tuple of integers.
    Example:
        'M41L+T215Y' → (41, 215)
    """
    parts = mutation_key.split("+")
    positions = []

    for p in parts:
        # Extract the numeric portion from mutations like M41L, T215Y, E40F
        match = re.search(r"(\d+)", p)
        if match:
            positions.append(int(match.group(1)))
        else:
            positions.append(99999)  # fallback if unexpected format

    return tuple(positions)

def sort_mutation_scores(mutation_scores):
    """
    Returns a new mutation_scores dict where each gene's mutations are sorted
    by numeric positions (biological ordering).
    """
    sorted_scores = {}

    for gene, mutations in mutation_scores.items():
        sorted_mutations = dict(
            sorted(
                mutations.items(),
                key=lambda x: parse_mutation_key(x[0])
            )
        )
        sorted_scores[gene] = sorted_mutations

    return sorted_scores

def extract_hivdb_version(json_path):
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)[0]
    try:
        version_text = data["drugResistance"][0]["version"].get("text", "")
        publish_date = data["drugResistance"][0]["version"].get("publishDate", "")
        return {
            "db_version": f"HIVDB {version_text}",
            "publish_date": publish_date
        }
    except Exception:
        return {
            "db_version": "HIVDB version unknown",
            "publish_date": "unknown"
        }

def remove_deprecated_drugs(mutation_scores, deprecated_drugs=None):
    """
    Remove entries for drugs that are no longer in use.
    """
    cleaned_scores = {}

    for gene, mutations in mutation_scores.items():
        cleaned_mutations = {}
        for mutation_key, drugs in mutations.items():
            cleaned_drugs = {drug: info for drug, info in drugs.items() if drug not in deprecated_drugs}
            if cleaned_drugs:
                cleaned_mutations[mutation_key] = cleaned_drugs
        if cleaned_mutations:
            cleaned_scores[gene] = cleaned_mutations

    return cleaned_scores

def parse_resistance_table(resistance_file, deprecated_drugs=None):
    """
    Parse resistance table CSV and remove deprecated drugs.
    """
    df_res = pd.read_csv(resistance_file)

    if deprecated_drugs:
        df_res = df_res[~df_res["Drug_abbr"].isin(deprecated_drugs)]
    return df_res

# ---------------------------------------------------------------------
# Batch processing
# ---------------------------------------------------------------------

def main():
    args = parser_args()

    # Remove from the report drugs that are deprecated or not used anymore
    deprecated_drugs = {d.strip() for d in args.deprecated_drugs.split(",")}

    # Build paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    asset_path = os.path.join(script_dir, "../assets")
    css_path = os.path.join(asset_path, "hiv_template_report.css")

    # --- Load CSS content
    with open(css_path, "r", encoding="utf-8") as css_file:
        css_content = css_file.read()

    # --- Load Jinja2 environment from the template's folder
    template_dir = os.path.abspath(asset_path)

    env = Environment(
        loader=FileSystemLoader(template_dir),
        autoescape=select_autoescape(['html', 'xml'])
    )
    template = env.get_template("hiv_template_report.html")

    logo_path = os.path.join(asset_path, "nf-core-viralrecon_logo_light.png")

    # Load image and encode to base64
    with open(logo_path, "rb") as f:
        logo_bytes = f.read()
        logo_b64 = base64.b64encode(logo_bytes).decode("utf-8")

    hivdb_version_info = extract_hivdb_version(args.sierralocal_json)

    df_mut = pd.read_csv(args.mutation_csv)
    if df_mut.empty:
        raise ValueError(f"No mutation data found in {args.mutation_csv}")

    sample_name = df_mut["Sample_name"].iloc[0]
    res_file = args.resistance_csv
    json_file = args.sierralocal_json
    nextclade_file = args.nextclade_csv
    consensus_file = args.consensus_fasta

    consensus_records = read_consensus_records(consensus_file)
    full_record = consensus_records["FULL"]
    ambiguous_site_summary = estimate_ambiguous_site_proportion(consensus_records)
    subtypes = get_nextclade_subtypes(nextclade_file, consensus_records)

    seq_summary = parse_sequence_summary(
        json_file,
        ivar_consensus_params=args.ivar_consensus_params,
        ivar_variant_maf=args.ivar_variant_maf,
        ambiguous_site_summary=ambiguous_site_summary,
    )

    regional_consensus_sequences = [
        {
            "region": region,
            "genes": region.split("_"),
            "fasta": record.format("fasta").strip(),
        }
        for region, record in consensus_records.items()
        if region != "FULL"
    ]

    df_res = parse_resistance_table(res_file, deprecated_drugs)

    mutation_scores_raw = extract_mutation_scoring(json_file)
    mutation_scores = sort_mutation_scores(mutation_scores_raw)
    mutation_scores = remove_deprecated_drugs(mutation_scores, deprecated_drugs)

    sample_data = {
        "sample_name": sample_name,
        "sequence_summary": seq_summary,
        "subtypes": subtypes,
        "mutation_data": df_mut.to_dict(orient="records"),
        "resistance_data": df_res.to_dict(orient="records"),
        "consensus_genome": full_record.format("fasta").strip(),
        "mutation_scores": mutation_scores,
        "regional_consensus_sequences": regional_consensus_sequences
    }

    # --- Render one HTML report for this sample
    html_content = template.render(
        sample = sample_data,
        hivdb_version = hivdb_version_info,
        nextclade_dataset_name = args.nextclade_dataset_name,
        nextclade_dataset_tag = args.nextclade_dataset_tag,
        pipeline_version = args.pipeline_version,
        date = date.today().strftime("%Y-%m-%d"),
        css_content = css_content,
        logo_b64 = logo_b64
    )

    with open(args.output_html, "w", encoding="utf-8") as f:
        f.write(html_content)

    print(f"✅ Report generated for {sample_name}: {args.output_html}")

if __name__ == "__main__":
    main()
