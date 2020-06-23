#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

regexes = {
    'nf-core/viralrecon': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'parallel-fastq-dump': ['v_parallel_fastq_dump.txt', r"parallel-fastq-dump\s:\s(\S+)"],
    'FastQC': ['v_fastqc.txt', r"FastQC\sv(\S+)"],
    'fastp': ['v_fastp.txt', r"fastp\s(\S+)"],
    'Bowtie 2': ['v_bowtie2.txt', r"bowtie2-align-s\sversion\s(\S+)"],
    'Samtools': ['v_samtools.txt', r"samtools\s(\S+)"],
    'BEDTools': ['v_bedtools.txt', r"bedtools\sv(\S+)"],
    'Mosdepth': ['v_mosdepth.txt', r"mosdepth\s(\S+)"],
    'Picard': ['v_picard.txt', r"\n(\S+)"],
    'iVar': ['v_ivar.txt', r"iVar\sversion\s(\S+)"],
    'VarScan 2': ['v_varscan.txt', r"VarScan\sv(\S+)"],
    'BCFTools': ['v_bcftools.txt', r"bcftools\s(\S+)"],
    'SnpEff': ['v_snpeff.txt', r"SnpEff\s(\S+)"],
    'SnpSift': ['v_snpsift.txt', r"SnpSift\sversion\s(\S+)"],
    'QUAST': ['v_quast.txt', r"QUAST\sv(\S+)"],
    'Cutadapt': ['v_cutadapt.txt', r"(\S+)"],
    'Kraken2': ['v_kraken2.txt', r"Kraken\sversion\s(\S+)"],
    'SPAdes': ['v_spades.txt', r"SPAdes\sgenome\sassembler\sv(\S+)"],
    'Unicycler': ['v_unicycler.txt', r"Unicycler\sv(\S+)"],
    'minia': ['v_minia.txt', r"Minia\sversion\s(\S+)"],
    'BLAST': ['v_blast.txt', r"blastn:\s(\S+)"],
    'ABACAS': ['v_abacas.txt', r"ABACAS.(\S+)"],
    'plasmidID': ['v_plasmidid.txt', r"(\S+)"],
    'Bandage': ['v_bandage.txt', r"Version:\s(\S+)"],
    'Minimap2': ['v_minimap2.txt', r"(\S+)"],
    'vg': ['v_vg.txt', r"vg\sversion\sv(\S+)"],
    'R': ['v_R.txt', r"R\sversion\s(\S+)"],
    'MultiQC': ['v_multiqc.txt', r"multiqc,\sversion\s(\S+)"]
}
results = OrderedDict()
results['nf-core/viralrecon'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['parallel-fastq-dump'] = '<span style="color:#999999;\">N/A</span>'
results['FastQC'] = '<span style="color:#999999;\">N/A</span>'
results['fastp'] = '<span style="color:#999999;\">N/A</span>'
results['Bowtie 2'] = '<span style="color:#999999;\">N/A</span>'
results['Samtools'] = '<span style="color:#999999;\">N/A</span>'
results['BEDTools'] = '<span style="color:#999999;\">N/A</span>'
results['Mosdepth'] = '<span style="color:#999999;\">N/A</span>'
results['Picard'] = '<span style="color:#999999;\">N/A</span>'
results['iVar'] = '<span style="color:#999999;\">N/A</span>'
results['VarScan 2'] = '<span style="color:#999999;\">N/A</span>'
results['BCFTools'] = '<span style="color:#999999;\">N/A</span>'
results['SnpEff'] = '<span style="color:#999999;\">N/A</span>'
results['SnpSift'] = '<span style="color:#999999;\">N/A</span>'
results['QUAST'] = '<span style="color:#999999;\">N/A</span>'
results['Cutadapt'] = '<span style="color:#999999;\">N/A</span>'
results['Kraken2'] = '<span style="color:#999999;\">N/A</span>'
results['SPAdes'] = '<span style="color:#999999;\">N/A</span>'
results['Unicycler'] = '<span style="color:#999999;\">N/A</span>'
results['minia'] = '<span style="color:#999999;\">N/A</span>'
results['BLAST'] = '<span style="color:#999999;\">N/A</span>'
results['ABACAS'] = '<span style="color:#999999;\">N/A</span>'
results['plasmidID'] = '<span style="color:#999999;\">N/A</span>'
results['Bandage'] = '<span style="color:#999999;\">N/A</span>'
results['Minimap2'] = '<span style="color:#999999;\">N/A</span>'
results['vg'] = '<span style="color:#999999;\">N/A</span>'
results['R'] = '<span style="color:#999999;\">N/A</span>'
results['MultiQC'] = '<span style="color:#999999;\">N/A</span>'

# Search each file using its regex
for k, v in regexes.items():
    try:
        with open(v[0]) as x:
            versions = x.read()
            match = re.search(v[1], versions)
            if match:
                results[k] = "v{}".format(match.group(1))
    except IOError:
        results[k] = False

# Remove software set to false in results
for k in list(results):
    if not results[k]:
        del(results[k])

# Dump to YAML
print ('''
id: 'software_versions'
section_name: 'nf-core/viralrecon Software Versions'
section_href: 'https://github.com/nf-core/viralrecon'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k,v))
print ("    </dl>")

# Write out regexes as csv file:
with open('software_versions.csv', 'w') as f:
    for k,v in results.items():
        f.write("{}\t{}\n".format(k,v))
