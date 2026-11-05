#!/usr/bin/env python3
import json
from pathlib import Path
from datetime import datetime
from Bio import SeqIO
import argparse


def parser_args(args=None):
    Description = (
        "Parse a GFF and its reference FASTA to generate a JSON profile file for codfreq"
        "with fragment definitions for HIV genes PR, RT, IN, etc."
    )
    Epilog = """Example usage:
    python gff2json.py --fasta NC_001802.1.fasta --gff NC_001802.1.gff --output NC_001802.1.json
    """

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)

    parser.add_argument(
        "-g",
        "--gff",
        type=str,
        required=True,
        help="Input GFF file with annotated features (required).",
    )
    parser.add_argument(
        "-f",
        "--fasta",
        type=str,
        required=True,
        help="Reference FASTA file corresponding to the GFF (required).",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="profile.json",
        help="Output JSON file path (default: 'profile.json').",
    )

    return parser.parse_args(args)


def parse_gff(gff_path):
    """Extract relevant feature coordinates from a GFF file."""
    fragments = []
    with open(gff_path) as gff:
        for line in gff:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue
            seqid, source, feature, start, end, score, strand, phase, attributes = cols
            start, end = int(start), int(end)

            if feature == "gene":
                attrs = {
                    kv.split("=")[0]: kv.split("=")[1]
                    for kv in attributes.split(";")
                    if "=" in kv
                }
                gene = attrs.get("gene", attrs.get("Name", attrs.get("product", "unknown"))).split()[0]
                fragments.append({
                    "fragmentName": gene,
                    "geneName": gene,
                    "start": start,
                    "end": end
                })
                print(f"Found feature: {gene} ({start}-{end})")
    return fragments


def build_json(reference_fasta, fragments, out_json):
    """Build and save JSON file with fragment configuration."""
    seq_record = SeqIO.read(reference_fasta, "fasta")
    ref_name = seq_record.id
    ref_seq = str(seq_record.seq).lower()
    today = datetime.today().strftime("%Y%m%d")

    fragment_config = []
    fragment_config.append({
        "fragmentName": ref_name,
        "refSequence": ref_seq
    })

    for frag in fragments:
        fragment_config.append({
            "fragmentName": frag["fragmentName"],
            "fromFragment": ref_name,
            "geneName": frag["geneName"],
            "refRanges": [[frag["start"], frag["end"]]]
        })

    data = {
        "version": today,
        "fragmentConfig": fragment_config
    }

    with open(out_json, "w") as f:
        json.dump(data, f, indent=4)

    print(f"✅ JSON file successfully written to: {out_json}")


def main():
    args = parser_args()
    fragments = parse_gff(args.gff)
    build_json(args.fasta, fragments, args.output)


if __name__ == "__main__":
    main()
