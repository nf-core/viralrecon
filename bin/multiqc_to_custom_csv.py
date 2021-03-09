#!/usr/bin/env python

import os
import sys
import errno
import argparse
import yaml

def parse_args(args=None):
    Description = 'Create custom spreadsheet for pertinent MultiQC metrics generated by the nf-core/viralrecon pipeline.'
    Epilog = "Example usage: python multiqc_to_custom_tsv.py"
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('-pl', '--platform', type=str, dest="PLATFORM", default='illumina', help="Sequencing platform for input data. Accepted values = 'illumina' or 'nanopore'  (default: 'illumina').")
    parser.add_argument('-md', '--multiqc_data_dir', type=str, dest="MULTIQC_DATA_DIR", default='multiqc_data', help="Full path to directory containing YAML files for each module, as generated by MultiQC. (default: 'multiqc_data').")
    parser.add_argument('-op', '--out_prefix', type=str, dest="OUT_PREFIX", default='summary', help="Full path to output prefix (default: 'summary').")
    return parser.parse_args(args)

def make_dir(path):
    if not len(path) == 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

# Find key in dictionary created from YAML file recursively
# From https://stackoverflow.com/a/37626981
def find_tag(d, tag):
    if tag in d:
        yield d[tag]
    for k,v in d.items():
        if isinstance(v, dict):
            for i in find_tag(v, tag):
                yield i

def yaml_fields_to_dict(yaml_file, append_dict={}, field_mapping_list=[], valid_sample_list=[]):
    integer_fields = [
        'mapped_passed', 'number_of_SNPs', 'number_of_indels', 'MISSENSE',
        '# contigs (>= 0 bp)', '# contigs (>= 5000 bp)', 'Largest contig'
    ]
    with open(yaml_file) as f:
        yaml_dict = yaml.safe_load(f)
        for k in yaml_dict.keys():
            key = k
            if os.path.basename(yaml_file).startswith('multiqc_picard_insertSize'):
                if k[-3:] == '_FR':
                    key = k[:-3]
            if os.path.basename(yaml_file).startswith('multiqc_cutadapt'):
                names = [x for x in valid_sample_list if key[:-2] == x]
                names += [x for x in valid_sample_list if key == x]
                if names != []:
                    key = names[0]
            include_sample = True
            if len(valid_sample_list) != 0 and key not in valid_sample_list:
                include_sample = False
            if include_sample:
                if key not in append_dict:
                    append_dict[key] = {}
                if field_mapping_list != []:
                    for i,j in field_mapping_list:
                        val = list(find_tag(yaml_dict[k], j[0]))
                        ## Fix for Cutadapt reporting reads/pairs as separate values
                        if j[0] == 'r_written' and len(val) == 0:
                            val = [list(find_tag(yaml_dict[k], 'pairs_written'))[0] * 2]
                        if len(val) != 0:
                            val = val[0]
                            if len(j) == 2:
                                val = list(find_tag(val, j[1]))[0]
                            if j[0] in integer_fields:
                                val = int(val)
                            if i not in append_dict[key]:
                                append_dict[key][i] = val
                            else:
                                print('WARNING: {} key already exists in dictionary so will be overwritten. YAML file {}.'.format(i,yaml_file))
                else:
                    append_dict[key] = yaml_dict[k]
    return append_dict


def metrics_dict_to_file(file_field_list, multiqc_data_dir, out_file, valid_sample_list=[]):
    metrics_dict = {}
    field_list = []
    for yaml_file,mapping_list in file_field_list:
        yaml_file = os.path.join(multiqc_data_dir, yaml_file)
        if os.path.exists(yaml_file):
            metrics_dict = yaml_fields_to_dict(yaml_file=yaml_file, append_dict=metrics_dict, field_mapping_list=mapping_list, valid_sample_list=valid_sample_list)
            field_list += [x[0] for x in mapping_list]
        else:
            print('WARNING: File does not exist: {}'.format(yaml_file))

    if metrics_dict != {}:
        make_dir(os.path.dirname(out_file))
        fout = open(out_file,'w')
        header = ['Sample'] + field_list
        fout.write('{}\n'.format(','.join(header)))
        for k in sorted(metrics_dict.keys()):
            row_list = [k]
            for field in field_list:
                if field in metrics_dict[k]:
                    row_list.append(metrics_dict[k][field])
                else:
                    row_list.append('NA')
            fout.write('{}\n'.format(','.join(map(str,row_list))))
        fout.close()
    return metrics_dict


def main(args=None):
    args = parse_args(args)

    ## File names for MultiQC YAML along with fields to fetch from each file
    illumina_variant_files = [
        ('multiqc_fastp.yaml',                                     [('# Input reads', ['before_filtering','total_reads']),
                                                                    ('# Trimmed reads (fastp)', ['after_filtering','total_reads'])]),
        ('multiqc_general_stats.yaml',                             [('% Non-host reads (Kraken 2)', ['PREPROCESS: Kraken 2_mqc-generalstats-preprocess_kraken_2-Unclassified'])]),
        ('multiqc_samtools_flagstat_samtools_bowtie2.yaml',        [('# Mapped reads', ['mapped_passed']),
                                                                    ('% Mapped reads', ['mapped_passed_pct'])]),
        ('multiqc_samtools_flagstat_samtools_ivar.yaml',           [('# Trimmed reads (iVar)', ['flagstat_total'])]),
        ('multiqc_samtools_flagstat_samtools_markduplicates.yaml', [('# Duplicate reads', ['duplicates_passed']),
                                                                    ('# Reads after MarkDuplicates', ['flagstat_total'])]),
        ('multiqc_picard_insertSize.yaml',                         [('Insert size mean', ['MEAN_INSERT_SIZE']),
                                                                    ('Insert size std dev', ['STANDARD_DEVIATION'])]),
        ('multiqc_picard_wgsmetrics.yaml',                         [('Coverage mean (picard)', ['MEAN_COVERAGE']),
                                                                    ('Coverage median (picard)', ['MEDIAN_COVERAGE']),
                                                                    ('Coverage std dev (picard)', ['SD_COVERAGE']),
                                                                    ('% Coverage > 10x (picard)', ['PCT_10X'])]),
        ('multiqc_general_stats.yaml',                             [('Coverage median', ['VARIANTS: mosdepth_mqc-generalstats-variants_mosdepth-median_coverage']),
                                                                    ('% Coverage > 1x', ['VARIANTS: mosdepth_mqc-generalstats-variants_mosdepth-1_x_pc']),
                                                                    ('% Coverage > 10x', ['VARIANTS: mosdepth_mqc-generalstats-variants_mosdepth-10_x_pc'])]),
        ('multiqc_bcftools_stats_bcftools_ivar.yaml',              [('# SNPs (iVar)', ['number_of_SNPs']),
                                                                    ('# INDELs (iVar)', ['number_of_indels'])]),
        ('multiqc_bcftools_stats_bcftools_bcftools.yaml',          [('# SNPs (BCFTools)', ['number_of_SNPs']),
                                                                    ('# INDELs (BCFTools)', ['number_of_indels'])]),
        ('multiqc_snpeff_snpeff_ivar.yaml',                        [('# Missense variants (iVar)', ['MISSENSE'])]),
        ('multiqc_snpeff_snpeff_bcftools.yaml',                    [('# Missense variants (BCFTools)', ['MISSENSE'])]),
        ('multiqc_quast_quast_ivar.yaml',                          [('# Ns per 100kb consensus (iVar)', ["# N's per 100 kbp"])]),
        ('multiqc_quast_quast_bcftools.yaml',                      [('# Ns per 100kb consensus (BCFTools)', ["# N's per 100 kbp"])])
    ]

    illumina_assembly_files = [
        ('multiqc_fastp.yaml',                                     [('# Input reads', ['before_filtering','total_reads'])]),
        ('multiqc_cutadapt.yaml',                                  [('# Trimmed reads (Cutadapt)', ['r_written'])]),
        ('multiqc_general_stats.yaml',                             [('% Non-host reads (Kraken 2)', ['PREPROCESS: Kraken 2_mqc-generalstats-preprocess_kraken_2-Unclassified'])]),
        ('multiqc_quast_quast_spades.yaml',                        [('# Contigs (SPAdes)', ['# contigs (>= 0 bp)']),
                                                                    ('Largest contig (SPAdes)', ['Largest contig']),
                                                                    ('% Genome fraction (SPAdes)', ['Genome fraction (%)']),
                                                                    ('N50 (SPAdes)', ['N50'])]),
        ('multiqc_quast_quast_unicycler.yaml',                     [('# Contigs (Unicycler)', ['# contigs (>= 0 bp)']),
                                                                    ('Largest contig (Unicycler)', ['Largest contig']),
                                                                    ('% Genome fraction (Unicycler)', ['Genome fraction (%)']),
                                                                    ('N50 (Unicycler)', ['N50'])]),
        ('multiqc_quast_quast_minia.yaml',                         [('# Contigs (minia)', ['# contigs (>= 0 bp)']),
                                                                    ('Largest contig (minia)', ['Largest contig']),
                                                                    ('% Genome fraction (minia)', ['Genome fraction (%)']),
                                                                    ('N50 (minia)', ['N50'])])
    ]

    nanopore_variant_files = [
        ('multiqc_samtools_flagstat.yaml',                         [('# Mapped reads', ['mapped_passed'])]),
        ('multiqc_general_stats.yaml',                             [('Coverage median', ['mosdepth_mqc-generalstats-mosdepth-median_coverage']),
                                                                    ('% Coverage > 1x', ['mosdepth_mqc-generalstats-mosdepth-1_x_pc']),
                                                                    ('% Coverage > 10x', ['mosdepth_mqc-generalstats-mosdepth-10_x_pc'])]),
        ('multiqc_bcftools_stats.yaml',                            [('# SNPs', ['number_of_SNPs']),
                                                                    ('# INDELs', ['number_of_indels'])]),
        ('multiqc_snpeff.yaml',                                    [('# Missense variants', ['MISSENSE'])]),
        ('multiqc_quast.yaml',                                     [('# Ns per 100kb consensus', ["# N's per 100 kbp"])])
    ]

    if args.PLATFORM == 'illumina':
        ## Dictionary of samples being single-end/paired-end
        is_pe_dict = {}
        yaml_file = os.path.join(args.MULTIQC_DATA_DIR,'multiqc_fastp.yaml')
        if os.path.exists(yaml_file):
            metrics_dict = yaml_fields_to_dict(yaml_file=yaml_file, append_dict={}, field_mapping_list=[('command', ['command'])], valid_sample_list=[])
            for sample,val in metrics_dict.items():
                if metrics_dict[sample]['command'].find('--out2') != -1:
                    is_pe_dict[sample] = True
                else:
                    is_pe_dict[sample] = False

        ## Write variant calling metrics to file
        metrics_dict_to_file(
            file_field_list   = illumina_variant_files,
            multiqc_data_dir  = args.MULTIQC_DATA_DIR,
            out_file          = args.OUT_PREFIX+'_variants_metrics_mqc.csv',
            valid_sample_list = is_pe_dict.keys()
        )

        ## Write de novo assembly metrics to file
        metrics_dict_to_file(
            file_field_list   = illumina_assembly_files,
            multiqc_data_dir  = args.MULTIQC_DATA_DIR,
            out_file          = args.OUT_PREFIX+'_assembly_metrics_mqc.csv',
            valid_sample_list = is_pe_dict.keys()
        )

    elif args.PLATFORM == 'nanopore':
        ## List of real samples to output in report
        sample_list = []
        yaml_file = os.path.join(args.MULTIQC_DATA_DIR,'multiqc_samtools_flagstat.yaml')
        if os.path.exists(yaml_file):
            metrics_dict = yaml_fields_to_dict(yaml_file=yaml_file, append_dict={}, field_mapping_list=[('# Mapped reads', ['mapped_passed'])], valid_sample_list=[])
            sample_list = metrics_dict.keys()
        
        metrics_dict_to_file(
            file_field_list   = nanopore_variant_files,
            multiqc_data_dir  = args.MULTIQC_DATA_DIR,
            out_file          = args.OUT_PREFIX+'_variants_metrics_mqc.csv',
            valid_sample_list = sample_list
        )

    else:
        print("Unrecognised option passed to --platform: {}. Accepted values = 'illumina' or 'nanopore'".format(args.PLATFORM))
        sys.exit(1)

if __name__ == '__main__':
    sys.exit(main())
