#!/usr/bin/env python3
"""Generate a per-sample HTML HIV drug resistance report.

Combines the Sierra-local JSON report, the mutation / resistance / Nextclade
tables and the consensus FASTA of one sample, and renders them with the
``hiv_template_report.html`` Jinja2 template stored in ``../assets``.
"""

from __future__ import annotations

import argparse
import base64
import json
import re
from datetime import date
from pathlib import Path
from typing import Any

import pandas as pd
from Bio import SeqIO
from jinja2 import Environment, FileSystemLoader, select_autoescape

from Bio.SeqRecord import SeqRecord

GENES = ("PR", "RT", "IN")
POL_REGIONS = ("PR_RT", "IN")
AMBIGUOUS_BASES = frozenset("RYSWKMBDHV")
PROTEIN_LENGTHS = {"PR": 99, "RT": 560, "IN": 288}
# Amino acid calls Sierra reports (as "...X") for codons that were not sequenced.
UNSEQUENCED_AAS = frozenset({"*ACDEFGHIKLMNPQRSTVWY", "ACDFGHILNPRSTVY"})
# Sort position for mutations without a numeric position (they go last).
UNKNOWN_POSITION = 99999

DRUG_CLASSES = (
    ("PI", "Protease Inhibitors (PI)"),
    ("NRTI", "Nucleoside Reverse Transcriptase Inhibitors (NRTI)"),
    ("NNRTI", "Non-nucleoside Reverse Transcriptase Inhibitors (NNRTI)"),
    ("INSTI", "Integrase Strand Transfer Inhibitors (INSTI)"),
)
MUTATION_SECTIONS = {
    "PR": (
        ("Major", "PI Major Mutations"),
        ("Accessory", "PI Accessory Mutations"),
        ("Other", "PR Other Mutations"),
    ),
    "RT": (
        ("NRTI", "NRTI Mutations"),
        ("NNRTI", "NNRTI Mutations"),
        ("Other", "RT Other Mutations"),
    ),
    "IN": (
        ("Major", "INSTI Major Mutations"),
        ("Accessory", "INSTI Accessory Mutations"),
        ("Other", "IN Other Mutations"),
    ),
}

LIST_MUTS_PATTERN = re.compile(r"\$listMutsIn\{(\d+)\(NOT ([A-Z]+)\)\}")
SUBSTITUTION_PATTERN = re.compile(r"([A-Z])(\d+)([A-Z*]+)")
TRIGGERED_PREFIX_PATTERN = re.compile(r"^([A-Z*_-]\d+)")


# ---------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------


def parse_args(args: list[str] | None = None) -> argparse.Namespace:
    """Parse command line arguments."""
    description = (
        "Parse Sierra-local JSON reports and corresponding resistance and mutation tables "
        "to generate an HTML report per sample."
    )
    epilog = """Example usage:
    python resistance_report.py --sierralocal_json SAMPLE_resistance.json --mutation_csv SAMPLE_mutation_table.csv
        --resistance_csv SAMPLE_resistance_table.csv --nextclade_csv SAMPLE_nextclade.csv
        --consensus_fasta SAMPLE.fa --ivar_consensus_params "-t 0.8 -q 30 -m 50 -n N"
        --ivar_variant_maf 0.01 --output_html SAMPLE_resistance_report.html
    """
    parser = argparse.ArgumentParser(
        description=description,
        epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-s", "--sierralocal_json", type=Path, required=True, help="Sierra-local JSON report.")
    parser.add_argument("-m", "--mutation_csv", type=Path, required=True, help="Mutation CSV file.")
    parser.add_argument(
        "-r",
        "--resistance_csv",
        type=Path,
        required=True,
        help="Resistance CSV file for one sample.",
    )
    parser.add_argument("-n", "--nextclade_csv", type=Path, required=True, help="Nextclade CSV file.")
    parser.add_argument("-cn", "--consensus_fasta", type=Path, required=True, help="Consensus FASTA file.")
    parser.add_argument(
        "-ic",
        "--ivar_consensus_params",
        required=True,
        help='Parameters used for ivar consensus calling, for example "-t 0.8 -q 30 -m 50 -n N".',
    )
    parser.add_argument(
        "-iv",
        "--ivar_variant_maf",
        type=float,
        default=0.01,
        help="Minor allele frequency threshold used for ivar variant calling.",
    )
    parser.add_argument(
        "-d",
        "--deprecated_drugs",
        default="",
        help=(
            "Comma-separated list of deprecated drugs that should be removed from the final report "
            "(for example: D4T,DDI,DPV,FPV/r,IDV/r,NFV,SQV/r,TPV/r)."
        ),
    )
    parser.add_argument(
        "-o", "--output_html", type=Path, required=True, help="Full path to the sample HTML report file."
    )
    parser.add_argument(
        "--nextclade_dataset_name",
        default="",
        help="Nextclade dataset name, for example neherlab/hiv-1.",
    )
    parser.add_argument(
        "--nextclade_dataset_tag",
        default="",
        help="Nextclade dataset tag, for example 2025-09-09--12-13-13Z.",
    )
    parser.add_argument("--pipeline_version", default="dev", help="nf-core/viralrecon pipeline version.")
    return parser.parse_args(args)


# ---------------------------------------------------------------------
# Input readers
# ---------------------------------------------------------------------


def load_sierra_json(json_path: Path) -> dict[str, Any]:
    """Return the first (and only) sequence entry of a Sierra-local JSON report."""
    with open(json_path, encoding="utf-8") as handle:
        return json.load(handle)[0]


def read_consensus_records(consensus_fasta: Path) -> dict[str, SeqRecord]:
    """Read the consensus FASTA into a dict keyed by the region after the last '|'."""
    records = {}
    for record in SeqIO.parse(consensus_fasta, "fasta"):
        _, separator, region = record.id.rpartition("|")
        if not separator or not region:
            raise ValueError(f"Consensus record ID must end with '|REGION': {record.id}")
        records[region] = record
    return records


def read_resistance_table(resistance_csv: Path, deprecated_drugs: set[str]) -> pd.DataFrame:
    """Read the resistance table and drop deprecated drugs."""
    df_res = pd.read_csv(resistance_csv)
    return df_res[~df_res["Drug_abbr"].isin(deprecated_drugs)]


# ---------------------------------------------------------------------
# Consensus-derived quality metrics
# ---------------------------------------------------------------------


def estimate_ambiguous_site_percent(consensus_records: dict[str, SeqRecord]) -> float:
    """Estimate the percentage of ambiguous sites across the pol gene.

    Ambiguous sites are IUPAC mixed-base codes (R, Y, S, W, K, M, B, D, H and V).
    The denominator is the number of non-N sites in the available PR_RT and IN regions.
    The result is rounded to two decimals, as displayed in the report.
    """
    ambiguous_sites = 0
    non_n_sites = 0
    for region in POL_REGIONS:
        if region not in consensus_records:
            continue
        sequence = str(consensus_records[region].seq).upper()
        ambiguous_sites += sum(base in AMBIGUOUS_BASES for base in sequence)
        non_n_sites += len(sequence) - sequence.count("N")

    proportion = ambiguous_sites / non_n_sites if non_n_sites else 0.0
    return float(f"{proportion * 100:.2f}")


def get_nextclade_subtypes(nextclade_csv: Path, consensus_records: dict[str, SeqRecord]) -> list[dict[str, str]]:
    """Return the Nextclade clade of each pol region present in the consensus."""
    clades = pd.read_csv(nextclade_csv, sep=";").set_index("seqName")["clade"]
    return [
        {"region": region, "subtype": clades.loc[record.id]}
        for region, record in consensus_records.items()
        if region in POL_REGIONS
    ]


# ---------------------------------------------------------------------
# Sierra-local JSON parsers
# ---------------------------------------------------------------------


def extract_hivdb_version(data: dict[str, Any]) -> dict[str, str]:
    """Return the HIVDB version and publish date used by Sierra."""
    try:
        version = data["drugResistance"][0]["version"]
        return {
            "db_version": f"HIVDB {version.get('text', '')}",
            "publish_date": version.get("publishDate", ""),
        }
    except (KeyError, IndexError, AttributeError):
        return {"db_version": "HIVDB version unknown", "publish_date": "unknown"}


def parse_sequence_warnings(data: dict[str, Any]) -> list[str]:
    """Return validation messages displayed in the report's warning box."""
    return [
        result["message"].strip()
        for result in data.get("validationResults", [])
        if result["level"] in ("WARNING", "NOTE") and result["message"].strip()
    ]


def parse_gene_stop_codons(data: dict[str, Any]) -> dict[str, list[str]]:
    """Return the stop codons (for example ``W26*``) found in each gene."""
    gene_stops: dict[str, list[str]] = {}
    for entry in data.get("alignedGeneSequences", []):
        gene = entry["gene"]["name"]
        if gene not in GENES:
            continue
        stops = [
            f"{mutation['consensus']}{mutation['position']}*"
            for mutation in entry.get("mutations", [])
            if mutation.get("hasStop") and not mutation["text"].endswith("X")
        ]
        if stops:
            gene_stops[gene] = list(dict.fromkeys(gene_stops.get(gene, []) + stops))
    return gene_stops


def _expand_list_muts(text: str, mutation: dict[str, Any]) -> str:
    """Replace Sierra's ``$listMutsIn{pos(NOT XY)}`` placeholders with the called mutations."""

    def replace(match: re.Match) -> str:
        position, excluded = match.groups()
        consensus = mutation["text"][0]
        # Drug-score mutations may omit AAs; their text still contains the call.
        called_aas = mutation.get("AAs") or re.sub(r"^[A-Z*]\d+", "", mutation["text"])
        amino_acids = "".join(aa for aa in called_aas if aa not in excluded and aa != consensus)
        return f"{consensus}{position}{amino_acids}"

    return LIST_MUTS_PATTERN.sub(replace, text)


def parse_mutation_comments(data: dict[str, Any]) -> dict[str, dict[str, list[str]]]:
    """Group drug-score mutation comments by gene and primaryType, avoiding repeats."""
    comments: dict[str, dict[str, list[str]]] = {}
    for entry in data.get("drugResistance", []):
        gene = entry["gene"]["name"]
        if gene not in GENES:
            continue
        mutations = (
            mutation
            for score in entry.get("drugScores", [])
            for partial in score.get("partialScores", [])
            for mutation in partial.get("mutations", [])
        )
        for mutation in mutations:
            if mutation["text"].endswith("X"):
                continue
            for comment in mutation.get("comments", []):
                text = comment.get("text", "").strip()
                if not text:
                    continue
                text = _expand_list_muts(text, mutation)
                texts = comments.setdefault(gene, {}).setdefault(mutation["primaryType"], [])
                if text not in texts:
                    texts.append(text)
    return comments


def _triggered_mutation_label(mutation: dict[str, Any]) -> str:
    """Return the mutation label using only the amino acids that triggered Sierra scoring."""
    text = mutation.get("text", "")
    triggered_aas = mutation.get("triggeredAAs", "")
    match = TRIGGERED_PREFIX_PATTERN.match(text)
    if not match or not triggered_aas:
        return text
    return f"{match.group(1)}{triggered_aas}"


def _mutation_sort_key(mutation_key: str) -> tuple[int, ...]:
    """Convert a key such as ``M41L+T215Y`` into its positions ``(41, 215)`` for sorting."""
    positions = []
    for part in mutation_key.split("+"):
        match = re.search(r"\d+", part)
        positions.append(int(match.group()) if match else UNKNOWN_POSITION)
    return tuple(positions)


def extract_mutation_scores(data: dict[str, Any], deprecated_drugs: set[str]) -> dict[str, dict[str, dict[str, Any]]]:
    """Group Sierra partial scores into one table per gene and drug class.

    Returns ``{gene: {drug_class: {"drugs": [...], "rows": [{"mutation", "scores"}]}}}``
    with mutations in positional order and drugs in alphabetical order.
    """
    scores_by_gene: dict[str, dict[str, dict[str, dict[str, Any]]]] = {}
    for protein_resistance in data.get("drugResistance", []):
        gene = protein_resistance.get("gene", {}).get("name")
        if not gene:
            continue
        gene_scores = scores_by_gene.setdefault(gene, {})
        for drug in protein_resistance.get("drugScores", []):
            drug_name = drug.get("drug", {}).get("displayAbbr", "")
            drug_class = drug.get("drugClass", {}).get("name", "")
            for partial in drug.get("partialScores", []):
                mutation_key = "+".join(_triggered_mutation_label(m) for m in partial.get("mutations", []))
                # Register the key even for deprecated drugs so ties in the positional
                # sort keep Sierra's order.
                drug_scores = gene_scores.setdefault(mutation_key, {})
                if drug_name not in deprecated_drugs:
                    drug_scores[drug_name] = {
                        "score": partial.get("score"),
                        "drug_class": drug_class,
                    }

    tables: dict[str, dict[str, dict[str, Any]]] = {}
    for gene, gene_scores in scores_by_gene.items():
        for mutation_key in sorted(gene_scores, key=_mutation_sort_key):
            for drug_name, info in gene_scores[mutation_key].items():
                table = tables.setdefault(gene, {}).setdefault(info["drug_class"], {"drugs": set(), "rows": {}})
                table["drugs"].add(drug_name)
                table["rows"].setdefault(mutation_key, {})[drug_name] = info["score"]

    for gene_tables in tables.values():
        for table in gene_tables.values():
            table["drugs"] = sorted(table["drugs"], key=str.lower)
            table["rows"] = [{"mutation": mutation, "scores": scores} for mutation, scores in table["rows"].items()]
    return tables


# ---------------------------------------------------------------------
# Sequence summary
# ---------------------------------------------------------------------


def _format_ranges(positions: list[int]) -> list[str]:
    """Collapse sorted positions into ``["5", "10-12"]`` style ranges."""
    ranges = []
    start = end = positions[0]
    for position in positions[1:]:
        if position == end + 1:
            end = position
        else:
            ranges.append(f"{start}" if start == end else f"{start}-{end}")
            start = end = position
    ranges.append(f"{start}" if start == end else f"{start}-{end}")
    return ranges


def summarize_gene_coverage(data: dict[str, Any]) -> list[dict[str, Any]]:
    """Summarize the codons covered by each gene and which codons are missing.

    Codons Sierra reports as unsequenced (``...X``) at the start or end of a gene are
    trimmed from the covered range; every missing stretch is reported in ``missing``.
    """
    coverage = []
    for entry in data.get("alignedGeneSequences", []):
        gene = entry["gene"]["name"]
        if gene not in GENES:
            continue
        first_aa, last_aa = entry["firstAA"], entry["lastAA"]
        true_end = PROTEIN_LENGTHS[gene]

        missing_positions = sorted(
            mutation["position"]
            for mutation in entry.get("mutations", [])
            if mutation.get("text", "").endswith("X")
            and mutation.get("AAs", "") in UNSEQUENCED_AAS
            and mutation.get("position")
            and isinstance(mutation["position"], int)
        )

        missing_parts = []
        if first_aa > 1:
            missing_parts.append("1" if first_aa == 2 else f"1-{first_aa - 1}")

        # Trim unsequenced codons from the start of the covered range.
        adjusted_first = first_aa
        for position in missing_positions:
            if position != adjusted_first:
                break
            adjusted_first += 1

        # Trim unsequenced codons from the end of the covered range.
        adjusted_last = last_aa
        for position in reversed(missing_positions):
            if position == true_end - 1:
                # The last codon is missing too when the one before it is unsequenced.
                missing_positions.append(true_end)
                adjusted_last = position - 1
            if position < adjusted_last:
                break
            adjusted_last = position - 1

        if missing_positions:
            missing_parts.extend(_format_ranges(missing_positions))

        if adjusted_last < true_end:
            tail_range = f"{adjusted_last + 1}-{true_end}"
            if tail_range not in missing_parts:
                missing_parts.append(tail_range)

        missing_parts.sort(key=lambda part: int(part.split("-")[0]))
        coverage.append(
            {
                "gene": gene,
                "first_codon": adjusted_first,
                "last_codon": adjusted_last,
                "missing": ", ".join(missing_parts),
            }
        )
    return coverage


def _ivar_option(params: str, flag: str, value_pattern: str) -> str:
    """Return the value of one ivar consensus option, for example ``-t`` -> ``0.8``."""
    match = re.search(rf"{flag}\s*({value_pattern})", params)
    if match is None:
        raise ValueError(f"Option '{flag}' not found in ivar consensus parameters: {params!r}")
    return match.group(1)


def build_settings(ivar_consensus_params: str, ivar_variant_maf: float) -> list[dict[str, str]]:
    """Return the calling thresholds shown in the report's Settings block."""
    threshold = float(_ivar_option(ivar_consensus_params, "-t", r"[\d.]+"))
    min_quality = int(_ivar_option(ivar_consensus_params, "-q", r"\d+"))
    min_depth = int(_ivar_option(ivar_consensus_params, "-m", r"\d+"))
    return [
        {"label": "Minimum read depth", "value": f"≥{min_depth}"},
        {
            "label": "Nucleotide mixture threshold (NMT)",
            "value": f"≥{ivar_variant_maf * 100:.0f}%",
        },
        {"label": "Minimum quality threshold", "value": str(min_quality)},
        {
            "label": "Mutation detection threshold (MDT)",
            "value": f"≥{(1 - threshold) * 100:.0f}%",
        },
    ]


def group_apobec_drms(mutation_rows: list[dict[str, Any]]) -> dict[str, list[str]]:
    """Combine flagged substitutions by gene and position for the APOBEC warning."""
    grouped: dict[str, dict[Any, set[str] | None]] = {}
    for row in mutation_rows:
        if not row.get("isApobecDRM"):
            continue
        calls = grouped.setdefault(row["Gene_name"], {})
        text = row["Mutations"]
        match = SUBSTITUTION_PATTERN.fullmatch(text)
        if match and not row.get("isInsertion") and not row.get("isDeletion"):
            consensus, position, amino_acids = match.groups()
            calls.setdefault((consensus, position), set()).update(amino_acids)
        else:
            # Preserve non-substitution notation without merging its characters.
            calls.setdefault(text, None)
    return {
        gene: [
            key if amino_acids is None else f"{key[0]}{key[1]}{''.join(sorted(amino_acids))}"
            for key, amino_acids in calls.items()
        ]
        for gene, calls in grouped.items()
    }


# ---------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------


def main(args: list[str] | None = None) -> None:
    args = parse_args(args)

    # Drugs that are deprecated or not used anymore are removed from the report
    deprecated_drugs = {drug.strip() for drug in args.deprecated_drugs.split(",") if drug.strip()}

    assets_dir = Path(__file__).resolve().parent.parent / "assets"
    env = Environment(
        loader=FileSystemLoader(assets_dir),
        autoescape=select_autoescape(["html", "xml"]),
    )
    template = env.get_template("hiv_template_report.html")
    css_content = (assets_dir / "hiv_template_report.css").read_text(encoding="utf-8")
    logo_b64 = base64.b64encode((assets_dir / "nf-core-viralrecon_logo_light.png").read_bytes()).decode("utf-8")

    df_mut = pd.read_csv(args.mutation_csv)
    if df_mut.empty:
        raise ValueError(f"No mutation data found in {args.mutation_csv}")
    mutation_rows = df_mut.to_dict(orient="records")
    for row in mutation_rows:
        match = SUBSTITUTION_PATTERN.fullmatch(row["Mutations"])
        row["unchanged_aa"] = bool(
            match
            and match.group(1) == match.group(3)
            and not row.get("isInsertion")
            and not row.get("isDeletion")
        )
    sample_name = df_mut["Sample_name"].iloc[0]

    sierra = load_sierra_json(args.sierralocal_json)
    consensus_records = read_consensus_records(args.consensus_fasta)
    subtypes = get_nextclade_subtypes(args.nextclade_csv, consensus_records)
    subtype_by_region = {subtype["region"]: subtype["subtype"] for subtype in subtypes}

    sample_data = {
        "sample_name": sample_name,
        "settings": build_settings(args.ivar_consensus_params, args.ivar_variant_maf),
        "subtypes": subtypes,
        "subtypes_differ": (
            all(region in subtype_by_region for region in POL_REGIONS)
            and subtype_by_region["PR_RT"] != subtype_by_region["IN"]
        ),
        "ambiguous_site_percent": estimate_ambiguous_site_percent(consensus_records),
        "gene_coverage": summarize_gene_coverage(sierra),
        "sequence_warnings": parse_sequence_warnings(sierra),
        "gene_stop_codons": parse_gene_stop_codons(sierra),
        "mutation_data": mutation_rows,
        "gene_apobec_drms": group_apobec_drms(mutation_rows),
        "mutation_comments": parse_mutation_comments(sierra),
        "resistance_data": read_resistance_table(args.resistance_csv, deprecated_drugs).to_dict(orient="records"),
        "mutation_scores": extract_mutation_scores(sierra, deprecated_drugs),
        "consensus_genome": consensus_records["FULL"].format("fasta").strip(),
        "regional_consensus_sequences": [
            {"genes": region.split("_"), "fasta": record.format("fasta").strip()}
            for region, record in consensus_records.items()
            if region != "FULL"
        ],
    }

    html_content = template.render(
        sample=sample_data,
        genes=GENES,
        drug_classes=DRUG_CLASSES,
        mutation_sections=MUTATION_SECTIONS,
        hivdb_version=extract_hivdb_version(sierra),
        nextclade_dataset_name=args.nextclade_dataset_name,
        nextclade_dataset_tag=args.nextclade_dataset_tag,
        pipeline_version=args.pipeline_version,
        date=date.today().isoformat(),
        css_content=css_content,
        logo_b64=logo_b64,
    )
    args.output_html.write_text(html_content, encoding="utf-8")

    print(f"✅ Report generated for {sample_name}: {args.output_html}")


if __name__ == "__main__":
    main()
