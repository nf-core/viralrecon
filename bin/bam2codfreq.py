#!/usr/bin/env python3

# Copyright (c) 2024 Stanford University HIVDB Team
# Modifications (c) 2025 Sarai Varona
#
# Original code licensed under the MIT License (see LICENSE in project root).
# This file includes modifications and is distributed under the MIT License.

"""
bam2codfreq.py

Unified script to convert SAM/BAM files into codon frequency tables (.codfreq)
using a profile JSON defining fragments and genes.

This script combines all helper modules and functions from the HIVdb CodFreq (https://github.com/hivdb/codfreq)
scripts into a single file. It computes codon frequencies and writes the final output to a CSV file.

Usage:
    python bam2codfreq.py --bam input.bam --profile profile.json --output output.codfreq [options]

Inputs:
    --bam          Path to mapping BAM file (required)
    --profile      Path to JSON profile with fragmentConfig (required)
    --output       Path to output codfreq file (required)

Optional arguments:
    --site-quality-cutoff  Minimum base quality to include a base (default: 0)

Example:
    python bam2codfreq.py --bam sample.bam --profile profile.json --output sample.codfreq
"""

import argparse
import json
import csv
from collections import defaultdict
from array import array
from typing import List, Tuple, Optional, Dict
import pysam

# -------------------------------
# Type aliases
# -------------------------------
NAPos = int
NAChar = int
AAPos = int
Header = str
CodonText = bytes

PosNA = Tuple[NAPos, int, NAChar, int]  # (refpos, insertion_index, base, quality)
PosCodon = Tuple[Header, AAPos, CodonText, int]  # fragment, aapos, codon, mean_quality
FragmentInterval = Tuple[List[Tuple[NAPos, NAPos]], Header]  # list of NA ranges + fragment name

# -------------------------------
# Codon translation table
# -------------------------------
CODON_TABLE: Dict[bytes, bytes] = {
    b'TTT': b'F', b'TTC': b'F', b'TTA': b'L', b'TTG': b'L',
    b'CTT': b'L', b'CTC': b'L', b'CTA': b'L', b'CTG': b'L',
    b'ATT': b'I', b'ATC': b'I', b'ATA': b'I', b'ATG': b'M',
    b'GTT': b'V', b'GTC': b'V', b'GTA': b'V', b'GTG': b'V',
    b'TCT': b'S', b'TCC': b'S', b'TCA': b'S', b'TCG': b'S',
    b'CCT': b'P', b'CCC': b'P', b'CCA': b'P', b'CCG': b'P',
    b'ACT': b'T', b'ACC': b'T', b'ACA': b'T', b'ACG': b'T',
    b'GCT': b'A', b'GCC': b'A', b'GCA': b'A', b'GCG': b'A',
    b'TAT': b'Y', b'TAC': b'Y',
    b'CAT': b'H', b'CAC': b'H', b'CAA': b'Q', b'CAG': b'Q',
    b'AAT': b'N', b'AAC': b'N', b'AAA': b'K', b'AAG': b'K',
    b'GAT': b'D', b'GAC': b'D', b'GAA': b'E', b'GAG': b'E',
    b'TGT': b'C', b'TGC': b'C', b'TGG': b'W',
    b'CGT': b'R', b'CGC': b'R', b'CGA': b'R', b'CGG': b'R',
    b'AGT': b'S', b'AGC': b'S', b'AGA': b'R', b'AGG': b'R',
    b'GGT': b'G', b'GGC': b'G', b'GGA': b'G', b'GGG': b'G',
    b'TAA': b'*', b'TGA': b'*', b'TAG': b'*',
}

# -------------------------------
# Helper functions
# -------------------------------

def iter_single_read_posnas(seq: str, qua: Optional[array], aligned_pairs: List[Tuple[Optional[int], Optional[int]]]) -> List[PosNA]:
    ENCODING = "UTF-8"
    GAP = ord(b"-")
    posnas: List[PosNA] = []
    seqchars = bytes(seq, ENCODING)
    prev_refpos = 0
    prev_seqpos0 = 0
    insidx = 0

    for seqpos0, refpos0 in aligned_pairs:
        if refpos0 is None:
            refpos = prev_refpos
            insidx += 1
        else:
            refpos = refpos0 + 1
            insidx = 0
            prev_refpos = refpos

        if seqpos0 is None:
            n = GAP
            q = qua[prev_seqpos0] if qua else 1
        else:
            n = seqchars[seqpos0]
            q = qua[seqpos0] if qua else 1
            prev_seqpos0 = seqpos0

        if refpos == 0:
            continue

        posnas.append((refpos, insidx, n, q))
    return posnas

def group_posnas_by_napos(posnas: List[PosNA]) -> List[Tuple[NAPos, List[PosNA]]]:
    by_napos: List[Tuple[NAPos, List[PosNA]]] = []
    prev_pos = -1
    for posna in posnas:
        if prev_pos == posna[0]:
            by_napos[-1][1].append(posna)
        else:
            prev_pos = posna[0]
            by_napos.append((posna[0], [posna]))
    return by_napos

def group_basepairs(posnas: List[PosNA], fragment_intervals: List[FragmentInterval]) -> List[Tuple[Header, List[Tuple[AAPos, List[PosNA]]]]]:
    posnas_by_napos = group_posnas_by_napos(posnas)
    basepairs: Dict[Header, List[Tuple[AAPos, List[PosNA]]]] = defaultdict(list)

    for napos, na_and_ins in posnas_by_napos:
        for frag_refranges, fragment_name in fragment_intervals:
            rel_napos0 = 0
            for start, end in frag_refranges:
                if start <= napos <= end:
                    aapos = (rel_napos0 + napos - start) // 3 + 1
                    basepairs[fragment_name].append((aapos, na_and_ins))
                rel_napos0 += end - start + 1
    return list(basepairs.items())

def get_comparable_codon(codon_posnas: List[List[PosNA]]) -> Tuple[CodonText, bool]:
    codon_chars: List[NAChar] = []
    for posnas in codon_posnas:
        for pna in posnas:
            codon_chars.append(pna[2])
    is_partial = len(codon_posnas) < 3
    return bytes(codon_chars), is_partial

def group_codons(basepairs: List[Tuple[Header, List[Tuple[AAPos, List[PosNA]]]]]) -> List[Tuple[Header, AAPos, List[List[PosNA]]]]:
    codons: List[Tuple[Header, AAPos, List[List[PosNA]]]] = []
    for fragment_name, fragment_bps in basepairs:
        prev_aapos = -1
        for aapos, na_and_ins in fragment_bps:
            if aapos == prev_aapos:
                codons[-1][2].append(na_and_ins)
            else:
                prev_aapos = aapos
                codons.append((fragment_name, aapos, [na_and_ins]))
    return codons

def posnas2poscodons(posnas: List[PosNA], fragment_intervals: List[FragmentInterval], read_refstart: int, read_refend: int, site_quality_cutoff: int) -> List[PosCodon]:
    fragments = [f for f in fragment_intervals if any(read_refstart <= end and read_refend >= start for start, end in f[0])]
    basepairs = group_basepairs(posnas, fragments)
    poscodons: List[PosCodon] = []

    for fragment_name, aapos, codon_posnas in group_codons(basepairs):
        codon, is_partial = get_comparable_codon(codon_posnas)
        if is_partial:
            continue
        totalq = sum(pna[3] for posnas in codon_posnas for pna in posnas)
        sizeq = sum(len(pna) for pna in codon_posnas)
        meanq = round(totalq / sizeq) if sizeq else 0
        if meanq < site_quality_cutoff:
            continue
        poscodons.append((fragment_name, aapos, codon, meanq))
    return poscodons

def iter_poscodons(samfile: str, fragment_intervals: List[FragmentInterval], site_quality_cutoff: int = 0) -> List[PosCodon]:
    poscodons_all: List[PosCodon] = []
    with pysam.AlignmentFile(samfile, "rb") as samfp:
        for read in samfp:
            if not read.query_sequence:
                continue
            posnas = iter_single_read_posnas(read.query_sequence, read.query_qualities, read.get_aligned_pairs(False))
            poscodons = posnas2poscodons(posnas, fragment_intervals, read.reference_start + 1, read.reference_end, site_quality_cutoff)
            poscodons_all.extend(poscodons)
    return poscodons_all

def expand_ambiguous_na(na: int) -> bytes:
    return AMBIGUOUS_NAS.get(na, bytes([na]))

def translate_codon(nas: bytes) -> str:
    """Translate a nucleotide sequence (A/T/C/G only) to amino acids.
    - Do not translate if contains '-' (deletion) or if the length is not a multiple of 3.
    """
    if b'-' in nas:
        return ""

    if len(nas) % 3 != 0:
        return ""

    aas = []
    for i in range(0, len(nas), 3):
        codon = nas[i:i+3]
        if codon in CODON_TABLE:
            aas.append(CODON_TABLE[codon].decode())
        else:
            aas.append("UNKNOWN")

    return "".join(aas)

# -------------------------------
# bam2codfreq function
# -------------------------------
def bam2codfreq_all(bamfile: str, profile_obj: dict, site_quality_cutoff: int = 0, output_codfreq: str = "output.codfreq"):
    """
    Process a BAM/SAM file and write codon frequencies to a CSV file (.codfreq).
    The output is ordered according to the fragment order in the profile JSON.
    """
    # Create fragment_intervals from JSON
    fragment_intervals: List[FragmentInterval] = []
    for frag in profile_obj["fragmentConfig"]:
        frag_name = frag["fragmentName"]
        ref_ranges = [(r[0], r[1]) for r in frag.get("refRanges", [])]
        fragment_intervals.append((ref_ranges, frag_name))

    # Preserve fragment order as in JSON
    frag_order = {frag_name: idx for idx, (_, frag_name) in enumerate(fragment_intervals)}

    # Generate codon positions from BAM/SAM
    poscodons = iter_poscodons(bamfile, fragment_intervals, site_quality_cutoff)

    # Count codons and accumulate quality
    codon_stats: Dict[Tuple[str, int, bytes], Dict[str, int]] = defaultdict(lambda: {"count": 0, "total_quality_score": 0})
    total_per_pos: Dict[Tuple[str, int], int] = defaultdict(int)

    for frag, aapos, codon, qual in poscodons:
        codon_stats[(frag, aapos, codon)]["count"] += 1
        codon_stats[(frag, aapos, codon)]["total_quality_score"] += qual
        total_per_pos[(frag, aapos)] += 1

    # Write .codfreq CSV file
    with open(output_codfreq, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["gene", "position", "total", "codon", "count", "total_quality_score", "aa_codon"])

        # Group by (fragment, position)
        grouped_positions = defaultdict(list)
        for (frag, aapos, codon), stats in codon_stats.items():
            grouped_positions[(frag, aapos)].append((codon, stats))

        # Iterate keeping fragment order
        for frag, aapos in sorted(total_per_pos.keys(), key=lambda x: (frag_order.get(x[0], 9999), x[1])):
            entries = grouped_positions[(frag, aapos)]

            # Sort by count descending
            entries_sorted = sorted(entries, key=lambda e: e[1]["count"], reverse=True)

            total = total_per_pos[(frag, aapos)]
            for codon, stats in entries_sorted:
                aa_codon = translate_codon(codon) if len(codon) % 3 == 0 else ""
                writer.writerow([frag, aapos, total, codon.decode("utf-8"), stats["count"], stats["total_quality_score"], aa_codon])

    print(f".codfreq file written to {output_codfreq}")

# -------------------------------
# MAIN
# -------------------------------
def main():
    parser = argparse.ArgumentParser(description="Convert SAM/BAM into codfreq CSV")
    parser.add_argument("--bam", required=True, help="Input BAM/SAM file")
    parser.add_argument("--profile", required=True, help="Profile JSON with fragmentConfig and sequenceAssemblyConfig")
    parser.add_argument("--output", required=True, help="Output codfreq CSV file")
    parser.add_argument("--site-quality-cutoff", type=int, default=0, help="Minimum base quality")

    args = parser.parse_args()

    with open(args.profile) as f:
        profile_obj = json.load(f)

    with pysam.AlignmentFile(args.bam, "rb") as bam_fp:
        bam_references = list(bam_fp.references)

    json_fragments = [frag.get("fromFragment") for frag in profile_obj.get("fragmentConfig", []) if frag.get("fromFragment")]

    if not json_fragments:
        raise ValueError("JSON profile file need a 'fromFragment' field to verify references with .bam files.")

    # Verificar si alguna referencia del BAM coincide con los 'fromFragment' del JSON
    unique_json_ref = set(json_fragments)
    if len(unique_json_ref) > 1:
        raise ValueError(f"JSON file contains multiple 'fromFragment': {unique_json_ref}. Only one expected")

    json_ref = list(unique_json_ref)[0]

    if len(set(bam_references)) > 1:
        raise ValueError(f"BAM file contains multiple references: {bam_references}. Only one expected.")

    bam_ref = bam_references[0]

    if bam_ref != json_ref:
        raise ValueError(
            f"Error: BAM reference ({bam_ref}) does not match JSON reference ({json_ref})."
        )
    else:
        print(f"✔ Reference check OK: {bam_ref} matches {json_ref}")

    bam2codfreq_all(args.bam, profile_obj, args.site_quality_cutoff, args.output)

if __name__ == "__main__":
    main()
