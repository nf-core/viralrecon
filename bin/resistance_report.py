#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import glob
import json
import argparse
import base64
import functools
import pandas as pd
from datetime import date
from Bio import SeqIO
from jinja2 import Environment, FileSystemLoader, select_autoescape


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

def parser_args(args=None):
    Description = "Parse Sierra-local JSON reports and corresponding resistance and mutation tables to generate an HTML report."
    Epilog = """Example usage:
    python resistance_report.py --sierralocal_folder resistance_jsons --mutation_folder mutation_tables --resistance_folder resistance_tables --nextclade_folder nextclade_folder  --consensus_folder consensus --ivar_consensus_params "-t 0.8 -q 30 -m 50 -n N" --output_html resistance_report.html
    """
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)

    parser.add_argument(
        "-s",
        "--sierralocal_folder",
        type=str,
        default="./sierralocal_json",
        help="Folder containing sierra-local JSON reports (default: ./sierralocal_json)",
    )
    parser.add_argument(
        "-m",
        "--mutation_folder",
        type=str,
        default="./mutation_tables",
        help="Folder containing mutation CSV files (default: ./mutation_tables)",
    )
    parser.add_argument(
        "-r",
        "--resistance_folder",
        type=str,
        default="./resistance_tables",
        help="Folder containing resistance CSV files (default: ./resistance_tables)",
    )
    parser.add_argument(
        "-n",
        "--nextclade_folder",
        type=str,
        default="./nextclade_reports",
        help="Folder containing nextclade CSV files (default: ./nextclade_reports)",
    )
    parser.add_argument(
        "-cn",
        "--consensus_folder",
        type=str,
        default="./consensus",
        help="Folder containing consensus files (default: ./consensus)",
    )
    parser.add_argument(
        "-gf",
        "--gff_folder",
        type=str,
        required=True,
        help="Folder containing GFF files with gene coordinates (required).",
    )
    parser.add_argument(
        "-ig",
        "--interest_genes",
        type=str,
        required=True,
        help="List of genes to extract, organized into groups. Genes separated by commas (',') will be included in the same output FASTA file. Gene groups separated by semicolons (';') will produce separate FASTA files.",
    )
    parser.add_argument(
        "-i",
        "--ivar_consensus_params",
        type=str,
        help="Parameters used for ivar consensus calling",
    )
    parser.add_argument(
        "-d",
        "--deprecated_drugs",
        type=str,
        default="D4T,DDI,DPV,FPV/r,IDV/r,NFV,SQV/r,TPV/r",
        help="Comma-separated list of deprecated drugs that should be removed from the final report (default: D4T,DDI,DPV,FPV/r,IDV/r,NFV,SQV/r,TPV/r)",
    )
    parser.add_argument(
        "-o",
        "--output_html",
        type=str,
        help="Full path to output HTML report file.",
    )
    return parser.parse_args(args)

def get_sample_number(sample_name, all_samples):
    """Return the sample index based on alphabetical order."""
    sorted_samples = sorted(all_samples)
    return sorted_samples.index(sample_name) + 1


def parse_sequence_summary(json_path, subtype_info=None, ivar_params=None):
    """
    Extract sequence summary information from Sierra-local JSON.
    - Lists each gene present (PR, RT, IN)
    - Detects missing nucleotide ranges from warnings
    - Detects missing codons at the start if firstAA != 1
    - Integrates subtype from Nextclade if available
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

    # --- Add subtype info from Nextclade if provided
    summary_lines.append(("subtype", subtype_info))

    # --- Parse ivar consensus parameters if provided
    # Extract numeric values with regex
    match_t = re.search(r"-t\s*([\d.]+)", ivar_params)
    match_q = re.search(r"-q\s*(\d+)", ivar_params)
    match_m = re.search(r"-m\s*(\d+)", ivar_params)

    t_val = float(match_t.group(1))
    q_val = int(match_q.group(1))
    m_val = int(match_m.group(1))

    summary_lines.append(f"Minimum read depth: ≥{m_val}")
    summary_lines.append(f"Mutation detection threshold (MDT): ≥{ (1-t_val) * 100:.0f}%")
    summary_lines.append(f"Minimum quality threshold: {q_val}")

    return summary_lines

def get_nextclade_subtype(nextclade_file, sample_name):
    """
    Extract the subtype (clade) assigned by Nextclade for the given sample.
    """
    if not nextclade_file or not os.path.exists(nextclade_file):
        return None
    try:
        df_next = pd.read_csv(nextclade_file, sep=';')
        # Attempt to find matching sample
        row = df_next[df_next["seqName"].str.contains(sample_name, case=False, na=False)]
        if not row.empty:
            subtype = row.iloc[0]["clade"]
            return subtype
    except Exception as e:
        print(f"⚠️ Could not parse Nextclade file for {sample_name}: {e}")
    return None

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
                mutation_ids = [m.get("text") for m in mutations]

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

def parse_interest_groups(interest_string):
    """
    Parse interest genes string into groups.
    Example:
        "PR,RT;IN" -> [["PR", "RT"], ["IN"]]
        "PR;RT;IN" -> [["PR"], ["RT"], ["IN"]]
        "PR,RT,IN" -> [["PR", "RT", "IN"]]
    """
    groups = []
    for group in interest_string.split(";"):
        genes = [g.strip() for g in group.split(",")]
        groups.append(genes)
    return groups


def read_gff_coordinates(gff_file):
    """
    Read GFF file and return a dict with gene_name -> (start, end, strand)
    Only parse entries with type 'gene'.
    """
    coords = {}

    with open(gff_file) as fh:
        for line in fh:
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            if feature_type != "gene":
                continue

            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attributes = fields[8]

            # Extract gene name from attributes (Name=XXXX)
            gene_name = None
            for attr in attributes.split(";"):
                if attr.startswith("Name="):
                    gene_name = attr.replace("Name=", "")
                    break

            if gene_name:
                coords[gene_name] = (start, end, strand)

    return coords


def extract_sequence(seq_record, start, end, strand):
    """
    Extract subsequence from a SeqRecord considering strand.
    GFF is 1-based inclusive.
    """
    subseq = seq_record.seq[start - 1:end]

    if strand == "-":
        subseq = subseq.reverse_complement()

    return subseq

def extract_protein_sequences(seq_record, coordinates, gene_groups, sample_name):
    """
    Extract sequences for groups of genes defined in gene_groups.

    Parameters
    ----------
    seq_record : Bio.SeqRecord
        Reference genome sequence loaded from FASTA.

    coordinates : dict
        Mapping of gene -> (start, end, strand) extracted from GFF.

    gene_groups : list of list
        List where each element is a group of genes to be extracted together.
        Example: [["PR","RT"], ["IN"]]

    Returns
    -------
    list of dicts
        [
            {
                "group_index": int,
                "genes": [list of genes],
                "fasta": "fasta block as string",
            },
            ...
        ]
    """

    protein_sequences = []

    for group_index, group in enumerate(gene_groups, start=1):

        # --- 1. Obtener coordenadas ordenadas del grupo ---
        valid_genes = [g for g in group if g in coordinates]

        if not valid_genes:
            continue

        coords = [(g, *coordinates[g]) for g in valid_genes]
        # ordenar por start
        coords.sort(key=lambda x: x[1])

        # --- 2. Verificar si son contiguas (end de uno == start del siguiente) ---
        merged = []
        current_genes = [coords[0][0]]
        current_start = coords[0][1]
        current_end = coords[0][2]
        current_strand = coords[0][3]

        for i in range(1, len(coords)):
            gene, start, end, strand = coords[i]

            # misma cadena y contiguos
            if strand == current_strand and start == current_end + 1:
                current_genes.append(gene)
                current_end = end
            else:
                # cerrar intervalo anterior
                merged.append((current_genes, current_start, current_end, current_strand))
                # iniciar uno nuevo
                current_genes = [gene]
                current_start = start
                current_end = end
                current_strand = strand

        # agregar el último intervalo
        merged.append((current_genes, current_start, current_end, current_strand))

        # --- 3. Extraer secuencias fusionadas ---
        fasta_lines = []

        for genes_list, start, end, strand in merged:
            seq = extract_sequence(seq_record, start, end, strand)

            gene_name = "_".join(genes_list)
            fasta_lines.append(f">{sample_name}_{gene_name}\n{seq}")

        fasta_block = "\n".join(fasta_lines)

        protein_sequences.append({
            "group_index": group_index,
            "genes": group,
            "fasta": fasta_block
        })

    return protein_sequences

# ---------------------------------------------------------------------
# Batch processing
# ---------------------------------------------------------------------

def main():
    args = parser_args()

    # Remove from the report drugs that are deprecated or not used anymore
    deprecated_drugs = {d.strip() for d in args.deprecated_drugs.split(",")}

    # Detect all files
    mutation_files = sorted(glob.glob(os.path.join(args.mutation_folder, "*_mutation_table.csv")))
    resistance_files = sorted(glob.glob(os.path.join(args.resistance_folder, "*_resistance_table.csv")))
    json_files = sorted(glob.glob(os.path.join(args.sierralocal_folder, "*_resistance.json")))
    nextclade_files = sorted(glob.glob(os.path.join(args.nextclade_folder, "*.csv")))
    consensus_files = sorted(glob.glob(os.path.join(args.consensus_folder, "*.fa")))
    gff_files = sorted(glob.glob(os.path.join(args.gff_folder, "*.gff")))

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

    hivdb_version_info = extract_hivdb_version(json_files[0]) if json_files else {"db_version": "HIVDB version unknown", "publish_date": "unknown"}

    all_samples_data = []

    for mut_file in mutation_files:
        df_mut = pd.read_csv(mut_file)
        if df_mut.empty:
            continue
        sample_name = df_mut["Sample_name"].iloc[0]

        # search corresponding resistance and JSON files
        res_file = next((r for r in resistance_files if sample_name in r), None)
        json_file = next((j for j in json_files if sample_name in j), None)
        nextclade_file = next((n for n in nextclade_files if sample_name in n), None)
        consensus_file = next((c for c in consensus_files if sample_name in c), None)
        gff_file = next((g for g in gff_files if sample_name in g), None)

        consensus_seq = None
        if consensus_file and os.path.exists(consensus_file):
            with open(consensus_file, "r", encoding="utf-8") as f:
                consensus_seq = f.read().strip()

        if not res_file or not json_file:
            print(f"⚠️ Skipping {sample_name}: missing resistance or JSON file")
            continue

        # --- Extract subtype from Nextclade
        subtype = get_nextclade_subtype(nextclade_file, sample_name)

        # --- Parse sequence summary
        seq_summary = parse_sequence_summary(json_file, subtype_info=subtype, ivar_params=args.ivar_consensus_params)

        # --- Parse resistance table
        df_res = parse_resistance_table(res_file, deprecated_drugs)

        # --- Extract mutation scoring from JSON
        mutation_scores_raw = extract_mutation_scoring(json_file)
        mutation_scores = sort_mutation_scores(mutation_scores_raw)

        # Remove deprecated drugs from mutation_scores
        mutation_scores = remove_deprecated_drugs(mutation_scores, deprecated_drugs)

        # Parse gene groups
        gene_groups = parse_interest_groups(args.interest_genes)

        # Read GFF coordinates
        coordinates = read_gff_coordinates(gff_file)

        # Read FASTA
        seq_record = next(SeqIO.parse(consensus_file, "fasta"))

        protein_sequences = extract_protein_sequences(seq_record, coordinates, gene_groups, sample_name)

        # Guardar toda la info en un dict
        all_samples_data.append({
            "sample_name": sample_name,
            "sequence_summary": seq_summary,
            "mutation_data": df_mut.to_dict(orient="records"),
            "resistance_data": df_res.to_dict(orient="records"),
            "consensus_genome": consensus_seq,
            "mutation_scores": mutation_scores,
            "protein_sequences": protein_sequences
        })

    # Ordenar alfabéticamente por nombre de muestra
    all_samples_data.sort(key=lambda x: x["sample_name"])

    # --- Render full HTML report
    html_content = template.render(
        all_samples=all_samples_data,
        hivdb_version=hivdb_version_info,
        date=date.today().strftime("%Y-%m-%d"),
        css_content=css_content,
        logo_b64=logo_b64
    )

    output_html = args.output_html or "all_samples_report.html"
    with open(output_html, "w", encoding="utf-8") as f:
        f.write(html_content)

    print(f"✅ Report generated: {output_html}")

if __name__ == "__main__":
    main()
