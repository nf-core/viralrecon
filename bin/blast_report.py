#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import base64
from io import BytesIO
import datetime
import os
import sys
import argparse
import csv
import textwrap
from jinja2 import Environment, FileSystemLoader, select_autoescape

def wrap_ylabels(ax, width=30):
    """Wrap long y-axis tick labels to avoid figure shrinking."""
    labels = [l.get_text() for l in ax.get_yticklabels()]
    wrapped = ["\n".join(textwrap.wrap(t, width)) for t in labels]
    ax.set_yticklabels(wrapped)

def parser_args(args=None):
    Description = "Generate HTML reports from BLAST results."
    Epilog = """Example usage:
    python blast_report.py --blast_file sample.filter.blastn.txt --fasta_file sample.scaffolds.fa --sample_name sample --id ticket_id --output_html sample_blast_report.html --output_fasta sample_reversed_filtered_contigs.fa --output_genotype sample_genotype.csv
    """
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)


    parser.add_argument(
        "-b",
        "--blast_file",
        required=True,
        type=str,
        help="BLAST results file (required)",
    )
    parser.add_argument(
        "-f",
        "--fasta_file",
        required=True,
        type=str,
        help="FASTA file (required)",
    )
    parser.add_argument(
        "-s",
        "--sample_name",
        required=True,
        type=str,
        help="Sample name (required)",
    )
    parser.add_argument(
        "-i",
        "--id",
        required=False,
        type=str,
        help="Run ID to be in the report (optional)",
    )
    parser.add_argument(
        "-ml",
        "--min_qlen",
        default=200,
        type=int,
        help="Minimum query length to consider a BLAST hit (optional)",
    )
    parser.add_argument(
        "-mc",
        "--min_coverage",
        default=50,
        type=int,
        help="Minimum coverage to consider a BLAST hit (optional)",
    )
    parser.add_argument(
        "-oh",
        "--output_html",
        required=True,
        type=str,
        help="Output HTML report file (required)",
    )
    parser.add_argument(
        "-of",
        "--output_fasta",
        required=True,
        type=str,
        help="Output FASTA file with reversed filtered contigs (required)",
    )
    parser.add_argument(
        "-og",
        "--output_genotype",
        required=True,
        type=str,
        help="CSV file to write predicted genotype summary",
    )
    parser.add_argument(
        "-sr",
        "--suggest_min_rows",
        type=int,
        default=20,
        help='Minimum number of rows (hits) required to consider auto-suggestion (default: 20)'
    )
    parser.add_argument(
        "-si",
        "--suggest_min_identity",
        type=float,
        default=90.0,
        help='Minimum max %% identity required to consider auto-suggestion (default: 90)'
    )
    parser.add_argument(
        "-sb",
        "--suggest_min_bitscore",
        type=float,
        default=400,
        help='Minimum max bitscore required to consider auto-suggestion (default: 400)'
    )
    parser.add_argument(
        "-d",
        "--dpi",
        type=int,
        default=300,
        help='DPI for output plots. Default: 300'
    )

    return parser.parse_args(args)

def extract_contig_headers(fasta_content):
    headers = []
    for line in fasta_content.split('\n'):
        if line.startswith('>'):
            headers.append(line)
    return headers

def encode_plot_to_base64(fig, dpi=300):
    buffer = BytesIO()
    fig.savefig(buffer, format='png', dpi=dpi, bbox_inches='tight', facecolor='white', edgecolor='none')
    buffer.seek(0)
    img_base64 = base64.b64encode(buffer.read()).decode('utf-8')
    buffer.close()
    return f"data:image/png;base64,{img_base64}"

def generate_report_data(df, fasta_file, sample_name, id, unique_contigs, suggest_enabled=True,
                         suggest_min_rows=20, suggest_min_identity=90.0, suggest_min_bitscore=400, plot_dpi=300):

    # Suggestion logic
    suggestion = "Automatic suggestion disabled."
    if suggest_enabled:
        try:
            max_pident = df.loc[df['pident'].idxmax()]
            max_bitscore = df.loc[df['bitscore'].idxmax()]
            grouped_pident_medians = df.groupby('sscinames')['pident'].median()
            highest_median_pident = grouped_pident_medians.idxmax()
            grouped_bitscore_medians = df.groupby('sscinames')['bitscore'].median()
            highest_median_bitscore = grouped_bitscore_medians.idxmax()
            species_counts = df['sscinames'].value_counts()
            top_species_count = species_counts.get(max_pident['sscinames'], 0)
            if (
                (max_pident['sscinames'] == max_bitscore['sscinames']) and
                (top_species_count >= suggest_min_rows) and
                (df['pident'].max() >= suggest_min_identity) and
                (df['bitscore'].max() >= suggest_min_bitscore) and
                (highest_median_pident == max_pident['sscinames']) and
                (highest_median_bitscore == max_bitscore['sscinames'])
            ):
                suggestion = str(max_pident['sscinames'])
            else:
                suggestion = "Please do manual assessment."
        except Exception as e:
            print(f"⚠️ Suggestion logic failed: {e}")
            suggestion = "Please do manual assessment."

    # --------------------------------------------------------------
    # Count number of contigs / genotypes and prepare plotting data
    n_contigs = df['contig'].nunique()
    u_contigs = df[['contig','coverage','sscinames']].drop_duplicates()
    n_sscinames = df['sscinames'].nunique()
    l_contigs = df[['contig','qlen','sscinames']].drop_duplicates()

    if n_contigs == 0 or n_sscinames == 0:
        raise ValueError(f"Invalid contig or genotype count for {sample_name}")

    pident_medians = df.groupby('sscinames')['pident'].median().sort_values(ascending=True)
    bitscore_medians = df.groupby('sscinames')['bitscore'].median().sort_values(ascending=True)
    coverage_medians = df.groupby('contig')['coverage'].median().sort_values(ascending=True)
    length_medians = df.groupby('contig')['qlen'].median().sort_values(ascending=True)

    pident_order = pident_medians.index.tolist()
    bitscore_order = bitscore_medians.index.tolist()
    coverage_order = coverage_medians.index.tolist()
    length_order = length_medians.index.tolist()

    if n_sscinames+n_contigs == 2:
        fig_height = 3
    elif n_sscinames+n_contigs > 8:
        fig_height = (n_sscinames+n_contigs)*0.8
    else:
        fig_height = n_sscinames+n_contigs

    h_ratio = 1 if n_contigs == 1 else n_contigs * 0.7
    legend_status = n_contigs > 1
    fig_width = 14 if legend_status else 10

    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = fig.add_gridspec(2, 2, height_ratios=[n_sscinames, h_ratio])

    # Plot 1: identity score
    ax = fig.add_subplot(gs[0 , 0])
    ax.set_title("Figure 1: Identity per genotype")
    sns.pointplot(
        data=df, x="pident", y="sscinames", hue="contig",
        dodge=.8 - .8 / n_contigs, palette="dark", errorbar=None,
        markers="d", markersize=6, linestyle="none", zorder=10, order=pident_order, legend=False
    )
    sns.stripplot(
        data=df, x="pident", y="sscinames", hue="contig",
        dodge=True, alpha=.4, legend=False, jitter=0.3, order=pident_order
    )
    if n_sscinames > 1:
        for i in range(n_sscinames):
            if i % 2 == 1:
                ax.axhspan(i - 0.5, i + 0.5, facecolor='gray', alpha=0.2, zorder=-1)
    ax.set_ylim(-0.5, n_sscinames-0.3)
    ax.grid(True, axis="x", linestyle="--")
    ax.set(xlabel='BLAST identity (%)', ylabel='Genotype')
    wrap_ylabels(ax, width=30)

    # Plot 2: bit score
    ax = fig.add_subplot(gs[0 , 1])
    ax.set_title("Figure 2: Bit-score per genotype")
    sns.pointplot(
        data=df, x="bitscore", y="sscinames", hue="contig",
        dodge=.8 - .8 / n_contigs, palette="dark", errorbar=None,
        markers="d", markersize=6, linestyle="none", zorder=10, order=bitscore_order, legend=legend_status
    )
    sns.stripplot(
        data=df, x="bitscore", y="sscinames", hue="contig",
        dodge=True, alpha=.4, legend=False, jitter=0.3, order=bitscore_order
    )
    if n_sscinames > 1:
        for i in range(n_sscinames):
            if i % 2 == 1:
                ax.axhspan(i - 0.5, i + 0.5, facecolor='gray', alpha=0.2, zorder=-1)
    ax.set_ylim(-0.5, n_sscinames-0.3)
    ax.grid(True, axis="x", linestyle="--")
    if legend_status:
        ax.legend().set_title("Contig")
        sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    ax.set(xlabel='BLAST Bit-score', ylabel=' ')
    wrap_ylabels(ax, width=30)

    # Plot 3: coverage
    ax = fig.add_subplot(gs[1, 0])
    ax.set_title("Figure 3: Coverage per contig")
    sns.pointplot(
        data=u_contigs, x="coverage", y="contig", hue="contig",
        dodge=.8 - .8 / n_contigs, palette="dark", errorbar=None,
        markers="d", markersize=6, linestyle="none", order=coverage_order
    )
    if n_contigs > 1:
        for i in range(n_contigs):
            if i % 2 == 0:
                ax.axhspan(i - 0.5, i + 0.5, facecolor='gray', alpha=0.2, zorder=-1)
    ax.set_ylim(-0.5, n_contigs-0.3)
    ax.grid(True, axis="x", linestyle="--")
    ax.set(xlabel='Coverage (x)', ylabel='Contig')
    ax.set_xscale('log')
    ax.set_xlim(left=50, right=10*df['coverage'].max())

    # Plot 4: length
    ax = fig.add_subplot(gs[1, 1])
    ax.set_title("Figure 4: Length per contig")
    sns.pointplot(
        data=l_contigs, x="qlen", y="contig", hue="contig",
        dodge=.8 - .8 / n_contigs, palette="dark", errorbar=None,
        markers="d", markersize=6, linestyle="none", order=length_order, legend=False
    )
    if n_contigs > 1:
        for i in range(n_contigs):
            if i % 2 == 0:
                ax.axhspan(i - 0.5, i + 0.5, facecolor='gray', alpha=0.2, zorder=-1)
    ax.set_ylim(-0.5, n_contigs-0.3)
    ax.grid(True, axis="x", linestyle="--")
    ax.set(xlabel='Length (bp)', ylabel=' ')
    ax.set_xlim(left=200, right=1.1*df['qlen'].max())

    fig.subplots_adjust(hspace=20, wspace=20)
    fig.tight_layout()

    img_buffer = BytesIO()
    fig.savefig(img_buffer, format='png', dpi=plot_dpi, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    img_buffer.seek(0)
    img_base64 = base64.b64encode(img_buffer.read()).decode('utf-8')
    img_data_uri = f"data:image/png;base64,{img_base64}"

    # Filter FASTA to only selected contigs
    filtered_fasta = filter_fasta_by_contigs(fasta_file, unique_contigs)

    contigs_content_html = filtered_fasta.replace('\n', '<br>')
    contig_headers = extract_contig_headers(filtered_fasta)
    contigs_summary = '<br>'.join(contig_headers)
    contig_count = len(contig_headers)

    # Warning
    warningtext = ''
    if df['pident'].max() < 90:
        warningtext = "&#9888; Highest identity < 90%"

    return {
        'sample_name': sample_name,
        'id': id,
        'time_stamp': datetime.datetime.now().strftime("%Y-%m-%d"),
        'is_error_report': False,
        'warningtext': warningtext,
        'img_data_uri': img_data_uri,
        'contigs_content': contigs_content_html,
        'contigs_summary': contigs_summary,
        'contig_count': contig_count,
        'suggestion': suggestion
    }


def generate_error_report_data(df, fasta_file, sample_name, id, unique_contigs, plot_dpi=300):

    file_name = os.path.basename(fasta_file)

    img_data_uri = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg=="
    contigs_content = ""
    contigs_summary = ""
    contig_count = 0

    try:
        n_contigs = df['contig'].nunique()
        u_contigs = df[['contig','coverage','sscinames']].drop_duplicates()
        n_sscinames = df['sscinames'].nunique()
        l_contigs = df[['contig','qlen','sscinames']].drop_duplicates()

        if n_contigs == 0 or n_sscinames == 0:
            print(f"Warning: Invalid contig or genotype count for {file_name}")
            return {
                'img_data_uri': img_data_uri,
                'contigs_content': contigs_content,
                'contigs_summary': contigs_summary,
                'contig_count': contig_count
            }

        pident_medians = df.groupby('sscinames')['pident'].median().sort_values(ascending=True)
        bitscore_medians = df.groupby('sscinames')['bitscore'].median().sort_values(ascending=True)
        coverage_medians = df.groupby('contig')['coverage'].median().sort_values(ascending=True)
        length_medians = df.groupby('contig')['qlen'].median().sort_values(ascending=True)

        pident_order = pident_medians.index.tolist()
        bitscore_order = bitscore_medians.index.tolist()
        coverage_order = coverage_medians.index.tolist()
        length_order = length_medians.index.tolist()

        if n_sscinames+n_contigs == 2:
            fig_height = 3
        elif n_sscinames+n_contigs > 8:
            fig_height = (n_sscinames+n_contigs)*0.8
        else:
            fig_height = n_sscinames+n_contigs

        h_ratio = 1 if n_contigs == 1 else n_contigs * 0.7
        legend_status = n_contigs > 1
        fig_width = 14 if legend_status else 10

        fig = plt.figure(figsize=(fig_width, fig_height))
        gs = fig.add_gridspec(2, 2, height_ratios=[n_sscinames, h_ratio])

        # (plots identical to normal path)
        ax = fig.add_subplot(gs[0 , 0])
        ax.set_title("Figure 1: Identity per genotype")
        sns.pointplot(
            data=df, x="pident", y="sscinames", hue="contig",
            dodge=.8 - .8 / n_contigs, palette="dark", errorbar=None,
            markers="d", markersize=6, linestyle="none", zorder=10, order=pident_order, legend=False
        )
        sns.stripplot(
            data=df, x="pident", y="sscinames", hue="contig",
            dodge=True, alpha=.4, legend=False, jitter=0.3, order=pident_order
        )
        if n_sscinames > 1:
            for i in range(n_sscinames):
                if i % 2 == 1:
                    ax.axhspan(i - 0.5, i + 0.5, facecolor='gray', alpha=0.2, zorder=-1)
        ax.set_ylim(-0.5, n_sscinames-0.3)
        ax.grid(True, axis="x", linestyle="--")
        ax.set(xlabel='BLAST identity (%)', ylabel='Genotype')

        ax = fig.add_subplot(gs[0 , 1])
        ax.set_title("Figure 2: Bit-score per genotype")
        sns.pointplot(
            data=df, x="bitscore", y="sscinames", hue="contig",
            dodge=.8 - .8 / n_contigs, palette="dark", errorbar=None,
            markers="d", markersize=6, linestyle="none", zorder=10, order=bitscore_order, legend=legend_status
        )
        sns.stripplot(
            data=df, x="bitscore", y="sscinames", hue="contig",
            dodge=True, alpha=.4, legend=False, jitter=0.3, order=bitscore_order
        )
        if n_sscinames > 1:
            for i in range(n_sscinames):
                if i % 2 == 1:
                    ax.axhspan(i - 0.5, i + 0.5, facecolor='gray', alpha=0.2, zorder=-1)
        ax.set_ylim(-0.5, n_sscinames-0.3)
        ax.grid(True, axis="x", linestyle="--")
        if legend_status:
            ax.legend().set_title("Contig")
            sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
        ax.set(xlabel='BLAST Bit-score', ylabel=' ')

        ax = fig.add_subplot(gs[1, 0])
        ax.set_title("Figure 3: Coverage per contig")
        sns.pointplot(
            data=u_contigs, x="coverage", y="contig", hue="contig",
            dodge=.8 - .8 / n_contigs, palette="dark", errorbar=None,
            markers="d", markersize=6, linestyle="none", order=coverage_order
        )
        if n_contigs > 1:
            for i in range(n_contigs):
                if i % 2 == 0:
                    ax.axhspan(i - 0.5, i + 0.5, facecolor='gray', alpha=0.2, zorder=-1)
        ax.set_ylim(-0.5, n_contigs-0.3)
        ax.grid(True, axis="x", linestyle="--")
        ax.set(xlabel='Coverage (x)', ylabel='Contig')
        ax.set_xscale('log')
        max_coverage = df['coverage'].max()
        if max_coverage > 0:
            ax.set_xlim(left=max(0.1, df['coverage'].min()/2), right=10*max_coverage)
        else:
            ax.set_xlim(left=0.1, right=100)

        ax = fig.add_subplot(gs[1, 1])
        ax.set_title("Figure 4: Length per contig")
        sns.pointplot(
            data=l_contigs, x="qlen", y="contig", hue="contig",
            dodge=.8 - .8 / n_contigs, palette="dark", errorbar=None,
            markers="d", markersize=6, linestyle="none", order=length_order, legend=False
        )
        if n_contigs > 1:
            for i in range(n_contigs):
                if i % 2 == 0:
                    ax.axhspan(i - 0.5, i + 0.5, facecolor='gray', alpha=0.2, zorder=-1)
        ax.set_ylim(-0.5, n_contigs-0.3)
        ax.grid(True, axis="x", linestyle="--")
        ax.set(xlabel='Length (bp)', ylabel=' ')
        max_length = df['qlen'].max()
        if max_length > 0:
            ax.set_xlim(left=max(1, df['qlen'].min()/2), right=1.1*max_length)
        else:
            ax.set_xlim(left=1, right=1000)

        fig.subplots_adjust(hspace = 20, wspace=20)
        fig.tight_layout()

        img_buffer = BytesIO()
        fig.savefig(img_buffer, format='png', dpi=plot_dpi, bbox_inches='tight',
                    facecolor='white', edgecolor='none')
        img_buffer.seek(0)
        img_base64 = base64.b64encode(img_buffer.read()).decode('utf-8')
        img_data_uri = f"data:image/png;base64,{img_base64}"

    except Exception as e:
        print(f"Error generating plots for {sample_name}: {e}")
        img_data_uri = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg=="

    # Read fasta file (unfiltered original)
    try:
        filtered_fasta = filter_fasta_by_contigs(fasta_file, unique_contigs)
        contigs_content_html = filtered_fasta.replace('\n', '<br>')
        if contigs_content_html:
                contig_headers = extract_contig_headers(filtered_fasta)
                contig_count = len(contig_headers)
                contigs_summary = '<br>'.join(contig_headers)
        else:
            print(f"Warning: FASTA file not found: {fasta_file}")
    except Exception as e:
        print(f"Error reading FASTA file for {file_name}: {e}")

    return {
        'sample_name': sample_name,
        'id': id,
        'time_stamp': datetime.datetime.now().strftime("%Y-%m-%d"),
        'is_error_report': True,
        'img_data_uri': img_data_uri,
        'contigs_content': contigs_content_html,
        'contigs_summary': contigs_summary,
        'contig_count': contig_count,
        'suggestion': "No recomendation due to lack of contigs meeting filtering criteria."
    }

def render_report(output_file, template, data, css_content, logo_b64):

    html_content = template.render(
        css_content=css_content,
        data = data,
        logo_b64=logo_b64
    )

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html_content)

def filter_fasta_by_contigs(fasta_path, contigs_to_keep):
    """Read a FASTA file and return only sequences whose headers match the contigs_to_keep list."""
    filtered_seqs = []
    keep = False
    header = None
    seq = []

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                # If we were collecting the previous sequence, save it
                if header and keep:
                    formatted_seq = format_sequence("".join(seq))
                    filtered_seqs.append((header, formatted_seq))

                header = line
                seq = []

                # Check if contig name appears in header
                # Example header: >NODE_1_length_7412_cov_139.22
                contig_name = header.replace('_reverse_complement', '')[1:]

                keep = contig_name in contigs_to_keep
            else:
                seq.append(line)

        # Save last seq
        if header and keep:
            # Format sequence with line breaks every 60 characters
            formatted_seq = format_sequence("".join(seq))
            filtered_seqs.append((header, formatted_seq))

    # Convert into FASTA format
    filtered_fasta = "\n".join(f"{h}\n{s}" for h, s in filtered_seqs)
    return filtered_fasta

def format_sequence(sequence, line_length=60):
    """Format a biological sequence with line breaks every specified number of characters."""
    if not sequence:
        return ""

    # Split sequence into chunks of line_length characters
    formatted_lines = []
    for i in range(0, len(sequence), line_length):
        formatted_lines.append(sequence[i:i+line_length])

    return "\n".join(formatted_lines)

def check_filters(blast_file, min_qlen=200, min_coverage=50):

    df = pd.read_csv(blast_file, sep="\t", header=0, index_col=0)

    unique_contigs = df['qaccver'].unique()
    df[['contig','temp1']] = df['qaccver'].str.split('_length_', expand=True)
    df[['length','coverage']] = df['temp1'].str.split('_cov_', expand=True)
    df['coverage'] = pd.to_numeric(df['coverage'])
    filtered_df = df[df['qlen'] > min_qlen]
    filtered_df = filtered_df[filtered_df['coverage'] > min_coverage]


    if filtered_df.empty or df['contig'].nunique() == 0:
        return df, True, unique_contigs
    else:
        unique_contigs = filtered_df['qaccver'].unique()
        return filtered_df, False, unique_contigs

def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(comp)[::-1]


def get_reverse_contigs(blast_file):
    contigs = set()
    with open(blast_file, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                sstart = int(row["sstart"])
                send = int(row["send"])
            except:
                continue
            if sstart > send:
                contigs.add(row["qaccver"])
    return contigs

def process_fasta(fasta_in, fasta_out, reverse_set):
    """Reverse-complement sequences whose ID is in reverse_set and update header."""
    def write_record(out, header, seq_lines):
        if not header:
            return
        seq = "".join(seq_lines)

        # Extract FASTA ID (everything after ">" up to first space)
        fasta_id = header.split()[0][1:]

        if fasta_id in reverse_set:
            seq = reverse_complement(seq)
            new_header = f">{fasta_id}_reverse_complement"
        else:
            new_header = header

        out.write(new_header + "\n")
        for i in range(0, len(seq), 80):
            out.write(seq[i:i+80] + "\n")

    with open(fasta_in) as fin, open(fasta_out, "w") as fout:
        header = None
        seq_lines = []

        for line in fin:
            line = line.rstrip("\n")
            if line.startswith(">"):
                write_record(fout, header, seq_lines)
                header = line
                seq_lines = []
            else:
                if line:
                    seq_lines.append(line.strip())

        write_record(fout, header, seq_lines)

def main(args=None, is_error=False):
    args = parser_args(args)

    # Build paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    asset_path = os.path.join(script_dir, "../assets")
    css_path = os.path.join(asset_path, "blast_report_template.css")

    # --- Load CSS content
    with open(css_path, "r", encoding="utf-8") as css_file:
        css_content = css_file.read()

    # --- Load Jinja2 environment from the template's folder
    template_dir = os.path.abspath(asset_path)

    env = Environment(
        loader=FileSystemLoader(template_dir),
        autoescape=select_autoescape(['html', 'xml'])
    )

    template = env.get_template("blast_report_template.html")

    logo_path = os.path.join(asset_path, "nf-core-viralrecon_logo_light.png")

    # Load image and encode to base64
    with open(logo_path, "rb") as f:
        logo_bytes = f.read()
        logo_b64 = base64.b64encode(logo_bytes).decode("utf-8")

    blast_results, is_error, unique_contigs = check_filters(args.blast_file, args.min_qlen, args.min_coverage)

    reverse_set = get_reverse_contigs(args.blast_file)

    if reverse_set:
        print(f"Reverse-complementing {len(reverse_set)} contig(s): {', '.join(reverse_set)}")
    else:
        print("No reverse-oriented contigs detected.")

    process_fasta(args.fasta_file, args.output_fasta, reverse_set)

    print(f"Output FASTA written to: {args.output_fasta}")

    if is_error:
        data = generate_error_report_data(blast_results, args.output_fasta, args.sample_name, args.id, unique_contigs, plot_dpi=args.dpi)
    else:
        data = generate_report_data(blast_results, args.output_fasta, args.sample_name, args.id, unique_contigs, plot_dpi=args.dpi)

    render_report(args.output_html, template, data, css_content, logo_b64)

    print(f"📁 Report saved to: {args.output_html}")

    # Optionally write genotype summary CSV (simple f.write style)
    if args.output_genotype:
        with open(args.output_genotype, "w") as f:
            f.write("Identifier,Sample,Genotype\n")
            f.write(f"{data['id']},{data['sample_name']},{data['suggestion']}\n")

        print(f"The predicted genotype saved to: {args.output_genotype}")

if __name__ == "__main__":
    main()
