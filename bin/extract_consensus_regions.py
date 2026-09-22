#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract gene regions from a consensus FASTA using GFF coordinates."
    )
    parser.add_argument(
        "-cn", "--consensus_fasta", required=True, help="Consensus FASTA file."
    )
    parser.add_argument(
        "-gf", "--gff", required=True, help="GFF file with gene coordinates."
    )
    parser.add_argument(
        "-ig",
        "--interest_genes",
        default="PR,RT;IN",
        help='Gene groups separated by semicolons (default: "PR,RT;IN").',
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output multi-FASTA containing the full and regional consensus sequences.",
    )
    return parser.parse_args()


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

    with open(gff_file, encoding="utf-8") as fh:
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


def extract_regions(seq_record, coordinates, gene_groups):
    regions = []

    for genes in gene_groups:
        missing = [gene for gene in genes if gene not in coordinates]
        if missing:
            raise ValueError(f"Genes not found in GFF: {', '.join(missing)}")

        gene_coordinates = sorted(
            [(gene, *coordinates[gene]) for gene in genes], key=lambda item: item[1]
        )
        start, end, strand = gene_coordinates[0][1:]

        for _, gene_start, gene_end, gene_strand in gene_coordinates[1:]:
            if gene_strand != strand or gene_start != end + 1:
                raise ValueError(f"Genes in group {','.join(genes)} are not contiguous")
            end = gene_end

        regions.append((genes, extract_sequence(seq_record, start, end, strand)))

    return regions


def main():
    args = parse_args()
    seq_record = SeqIO.read(args.consensus_fasta, "fasta")
    coordinates = read_gff_coordinates(args.gff)
    gene_groups = parse_interest_groups(args.interest_genes)
    regions = extract_regions(seq_record, coordinates, gene_groups)

    output_records = [
        SeqRecord(
            seq_record.seq,
            id=f"{seq_record.id}|FULL",
            description=""
        )
    ]

    for genes, sequence in regions:
        region_name = "_".join(genes)
        output_records.append(
            SeqRecord(
                sequence,
                id=f"{seq_record.id}|{region_name}",
                description=""
            )
        )

    SeqIO.write(output_records, args.output, "fasta")


if __name__ == "__main__":
    main()
