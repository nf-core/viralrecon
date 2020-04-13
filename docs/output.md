# nf-core/viralrecon: Output

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [FastQC](#fastqc) - read quality control
* [fastp](#fastp) - read quality filtering and adapter trimming.
* [Kraken2](#kraken2) - Mapping to host genome.
* Mapping + variant calling + consensus
  * [bowtie2](#bowtie2) - mapping against reference genomes.
  * [SAMtools](#samtools) - Mapping result processing and unmapped reads selection.
  * [Picard](#picard) - Enrichment and alignment metrics.
  * [VarScan](#varscan) - Variant calling.
  * [SnpEff and SnpSift] - Variant calling annotation.
  * [Bcftools](#bcftools) - Variant calling index and consensus genome generation.
  * [Bedtools](#bedtools) - Consensus genome masking.
* De novo assembly
  * [SPADES](#spades) - Viral genome assembly.
  * [MetaSPADES](#metaspades) - Viral genome assembly.
  * [Unicycler](#unicycler) - Viral genome assembly.
  * [QUAST](#quast) - Assembly quality assessment.
  * [Blast](#blast) - Blast alignment.
  * [PlasmidID](#plasmidid) - Visualization of the alignment.
  * [ABACAS](#abacas) - Contig ordering according to reference.
* [MultiQC](#multiqc) - aggregate report, describing results of the whole pipeline

## Preprocessing

### FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality. To see how your reads look after trimming, look at the FastQC reports in the `fastp` directory.

**Output directory: `results/fastqc`**

* `sample_fastqc.html`
  * FastQC report, containing quality metrics for your untrimmed raw fastq files
* `zips/sample_fastqc.zip`
  * zip file containing the FastQC report, tab-delimited data file and plot images

### Quality filtering and adaptor rrimming

[Fastp](https://github.com/OpenGene/fastp) is a tool designed to provide fast all-in-one preprocessing for FastQ files. This tool is developed in C++ with multithreading supported to afford high performance.

**Output directory: `preprocess/fastp`**
* `sample_T1_[12].trim.fastq.gz`
 * Only with `--save_trimmed` param. Paired trimmed reads.
* `sample_T1_[12].trim.fail.gz`
 * Only with `--save_trimmed` param. Unpaired trimmed reads.
* `sample_T1_[12].fastp.html`
 * Fastp report in html format.
* `sample_T1_[12].fastp.json`
 * Fastp report in json format.
* `fastqc/sample_T1_[12].trim.fastqc.html`
 * FastQC report of the trimmed reads.
* `fastqc/zips/sample_T1_[12].trim.fastqc.zip`
 * Zip file with the FastQC report.
 * `logs/sample_T1.fastp.log`
  * Trimming log file.

### Kraken2

[Kraken2](https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual)<a name="kraken"> </a><a href="#Kraken_reference">[3]</a> is a taxonomic sequence classifier that assigns taxonomic labels to DNA sequences. Kraken examines the k-mers within a query sequence and uses the information within those k-mers to query a database. That database maps k-mers to the lowest common ancestor (LCA) of all genomes known to contain a given k-mer.

We mapped the fastq file against the reference host genome.

**Output directory: `preprocess/kraken2`**

* `{sample_id}_T1.kraken2.report.txt`
  * Kraken taxonomic report. The report contains one line per taxonomic classification nd the following columns:
    1. Percentage of fragments covered by the clade rooted at this taxon
    2. Number of fragments covered by the clade rooted at this taxon
    3. Number of fragments assigned directly to this taxon
    4. A rank code, indicating (U)nclassified, (R)oot, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. Taxa that are not at any of these 10 ranks have a rank code that is formed by using the rank code of the closest ancestor rank with a number indicating the distance from that rank.  E.g., "G2" is a rank code indicating a taxon is between genus and species and the grandparent taxon is at the genus rank.
    5. NCBI taxonomic ID number
    6. Indented scientific name

## Mapping + variant calling + consensus

### Bowtie and Samtools
[Bowtie](http://bio-bwa.sourceforge.net/)<a name="bowtie"> </a><a href="#Bowtie_reference">[4]</a> is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences. It is particularly good at aligning reads of about 50 up to 100s of characters to relatively long genomes. Bowtie 2 indexes the genome with an FM Index (based on the Burrows-Wheeler Transform or BWT) to keep its memory footprint small. Bowtie 2 supports gapped, local, and paired-end alignment modes.

The result mapping files are further processed with [SAMtools](http://samtools.sourceforge.net/)<a name="samtools"> </a><a href="#SAMtools_reference">[5]</a>, sam format is converted to bam, sorted and an index .bai is generated.

We mapped the fastq file against the reference viral genome.

**Output directory: `variants/bowtie2`**

* `{sample_id}_T1_sorted.bam`
  * Sorted aligned bam file.
* `{sample_id}_T1_sorted.bam.bai`
  * Index file for soreted aligned bam.
* `samtools_stats/{sample_id}_T1.sorted.bam.flagstat`
  * Samtools flagstats mapping stats summary.

### Picard
[Picard](https://broadinstitute.github.io/picard/index.html)<a name="picard"> </a><a href="#Picard_reference">[6]</a> is a set of command line tools for manipulating high-throughput sequencing (HTS) data. In this case we used it to obtain mapping stats.

**Output directory: `variants/bowtie2/picard_metrics`**

* `{sample_id}_T1.CollectWgsMetrics.coverage_metrics`
  * Picard metrics summary file for evaluating coverage and performance.

Picard documentation: [Picarddocs](https://broadinstitute.github.io/picard/command-line-overview.html)

## Variant calling

In this pipeline to generate the consensus viral genome we use the approach of calling for variants between the mapped reads and the reference viral genome, and adding these variants to the reference viral genome.

### Variant calling

#### VarScan

First of all SAMtools is used to generate the variant calling VCF file. Then [VarScan](http://varscan.sourceforge.net/)<a name="varscan"> </a><a href="#VarScan_reference">[7]</a> is used to call for major and low frequency variants. VarScan is a platform-independent software tool developed at the Genome Institute at Washington University to detect variants in NGS data.

**Output directory: `variants/varscan2`**

* `{sample_id}_T1.highfreq.vcf.gz`
  * High frequency variants VCF file.
* `{sample_id}_T1.lowfreq.vcf.gz`
  * Low frequency variants VCF.
* `logs/{sample_id}_T1.highfreq.varscan2.log`
  * VarScan2 log file.

### Annotation

#### SnpEff and SnpSift

[SnpEff](http://snpeff.sourceforge.net/SnpEff.html)<a name="snpeff"> </a><a href="#SnpEff_reference">[8]</a> is a genetic variant annotation and functional effect prediction toolbox. It annotates and predicts the effects of genetic variants on genes and proteins (such as amino acid changes). [SnpSift](http://snpeff.sourceforge.net/SnpSift.html)</a><a href="#SnpSift_reference">[9]</a> annotates genomic variants using databases, filters, and manipulates genomic annotated variants. Once you annotated your files using SnpEff, you can use SnpSift to help you filter large genomic datasets in order to find the most significant variants for your experiment.

**Output directory: `variants/varscan2/snpeff`**

* `{sample_id}_T1.lowfreq.snpEff.csv`
  * Low frequency variants annotation csv file.
* `{sample_id}_T1.lowfreq.snpSift.table.txt`
  * Low frequency variants SnpSift summary table.
* `{sample_id}_T1.lowfreq.snpEff.vcf.gz`
  * Low frequency variants annotated VCF table.
* `{sample_id}_T1.lowfreq.snpEff.genes.txt`
  * Low frequency variants genes table.
* `{sample_id}_T1.lowfreq.snpEff.summary.html`
  * Low frequency variants summary html file.
* `{sample_id}_T1.highfreq.snpEff.csv`
  * High frequency variants annotation csv file.
* `{sample_id}_T1.highfreq.snpSift.table.txt`
  * High frequency variants SnpSift summary table.
* `{sample_id}_T1.highfreq.snpEff.vcf.gz`
  * High frequency variants annotated VCF table.
* `{sample_id}_T1.highfreq.snpEff.genes.txt`
  * High frequency variants genes table.
* `{sample_id}_T1.highfreq.snpEff.summary.html`
  * High frequency variants summary html file.

## Consensus genome

### BCFtools

[Bcftools](http://samtools.github.io/bcftools/bcftools.html)<a name="bcftools"> ftom the SAMtools project </a><a href="#SAMtools_reference">[4]</a> is a set of utilities that manipulate variant calls in the Variant Call Format (VCF) and its binary counterpart BCF. The resulting variant calling vcf for haploid genomes is indexed and then the consensus genome is created adding the variants to the reference viral genome. This consensus genome was obtained using the predominant variants (majority) of the mapping file.

**Output directory: `variants/bcftools`**

* `{sample_id}_T1.consensus.fa`
  * Consensus viral genome file generated from adding the variants called before to the viral reference genome. These variants are only the majoritarian variants, inlcuding only SNPs and small indels.

### Bedtools

[Bedtools](https://bedtools.readthedocs.io/en/latest/)<a name="bedtools"> </a><a href="#Bedtools_reference">[10]</a> are a swiss-army knife of tools for a wide-range of genomics analysis tasks. In this case we use:
  * bedtools genomecov computes histograms (default), per-base reports (-d) and BEDGRAPH (-bg) summaries of feature coverage (e.g., aligned sequences) for a given genome.
  * bedtools maskfasta masks sequences in a FASTA file based on intervals defined in a feature file. This may be useful fro creating your own masked genome file based on custom annotations or for masking all but your target regions when aligning sequence data from a targeted capture experiment.

  **Output directory: `variants/bcftools`**

* `{sample_id}_{reference_virus_name}.consensus.masked.fasta`
  * Masked consensus fasta file.

## De novo assembly

We selected the reads that didn't cluster using kraken2 with the host genome and assembled them to create a viral genome assembly.

### Cutadapt (*only amplicon*)

[Cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html) is used for clipping primer sequences from reads prior to assembly.

**Output directory:** `preprocess/cutadapt`
  * `fastqc`
    * Fastqc reports for primer clipped reads.
  * `logs`
    * Cutadapt logs.
  * `${sample}_[1,2].ptrim.fastq.gz`
    * only if `--save_trimmed`. Fastq files with primer sequences trimmed.

### SPAdes

[SPAdes](https://kbase.us/applist/apps/kb_SPAdes/run_SPAdes/release?gclid=Cj0KCQiAt_PuBRDcARIsAMNlBdroQS7y2hPFuhagq1QPvQ39FcvGxbhtZwhn8YbxIB4LrGIHKjJ-iPwaAn_lEALw_wcB) is a de Bruijn graph-based assembler. We selected the reads that didn't mapped with the host genome and assembled them using SPAdes to create a viral genome assembly.
 
**Output directory:** `assembly/spades`	
* `{sample_id}.scaffolds.fasta`	
  * Assembled scaffolds.	
  
  
### MetaSPAdes

[MetaSPAdes](https://kbase.us/applist/apps/kb_SPAdes/run_SPAdes/release?gclid=Cj0KCQiAt_PuBRDcARIsAMNlBdroQS7y2hPFuhagq1QPvQ39FcvGxbhtZwhn8YbxIB4LrGIHKjJ-iPwaAn_lEALw_wcB) is a de Bruijn graph-based assembler, with the option `--meta` the assembler works for metagenomics date trying to reconstruct different genomes. 	

**Output directory: `assembly/metaspades`**	
* `{sample_id}.meta.scaffolds.fasta`	
  * Assembled scaffolds.	


### Unicycler	

[Unicycler](https://github.com/rrwick/Unicycler) is an assembly pipeline that works as a spades optimiser.	

**Output directory: `assembly/unicycler`**	
* `{sample_id}.assembly.fasta`	
  * Assembled scaffolds.	

### QUAST	

[QUAST](http://bioinf.spbau.ru/quast) evaluates genome assemblies. We compared the reference genome with the contigs and scaffold assemblies. The html results can be opened with any browser (we recommend using Google Chrome). We have a quast folder for each assembler selected.	

**Output directory: `assembly/${assembler}`**	
* `quast/report.html`	
  * Compressed format of the indexed variants file.	
  * The meaning of the different metrics:	
    * Contigs (≥ x bp): is total number of contigs of length ≥ x bp.	
    * Total length (≥ x bp): is the total number of bases in contigs of length ≥ x bp.	
    * Contigs: is the total number of contigs in the assembly.	
    * Largest contig: is the length of the longest contig in the assembly.	
    * Total length: is the total number of bases in the assembly.	
    * Reference length: is the total number of bases in the reference genome.	
    * GC (%): is the total number of G and C nucleotides in the assembly, divided by the total length of the assembly.	
    * Reference GC (%): is the percentage of G and C nucleotides in the reference genome.	
    * N50: is the length for which the collection of all contigs of that length or longer covers at least half an assembly.	
    * NG50: is the length for which the collection of all contigs of that length or longer covers at least half the reference genome. This metric is computed only if the reference genome is provided.	
    * N75 and NG75: are defined similarly to N50 but with 75 % instead of 50 %.	
    * L50 (L75, LG50, LG75) is the number of contigs equal to or longer than N50 (N75, NG50, NG75). In other words, L50, for example, is the minimal number of contigs that cover half the assembly.	

### Blast alignments	

[NCBI Blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi) is used for aligning the contigs against the reference virus genome.	

**Output directory:** `assembly/${assembler}/blast`	
* {sample_id}_blast_filt_header.txt: blast results against the target virus.	
* {sample_id}_blast_bacteria_filt.txt: blast results for bacteria database.	

### PlasmidID	

[PlasmidID](https://github.com/BU-ISCIII/plasmidID) was used to graphically represent the alignment of the reference genome with the assembly obtained with SPAdes. This helps to visualize the coverage of the reference genome in the assembly. To find more information about the output files go to: https://github.com/BU-ISCIII/plasmidID/wiki/Understanding-the-image:-track-by-track	

**Output directory:** `assembly/${assembler}/plasmidid	
* `${sample_id}/images/{sample_id}_{reference_viral_genome}.png`	
  * PNG file with the visualization of the alignment between the assembled viral genome and the reference viral genome.	
* `${sample_id}/data/`	
  * Files used for drawing the circos images. 	
* `$sample_id/database`	
  * Annotation files used for drawing the circos images.	

### ABACAS

[Abacas](abacas.sourceforge.ne) intended to rapidly contiguate (align, order, orientate), visualize and design primers to close gaps on shotgun assembled contigs based on a reference sequence.	

**Output directory:** `assembly/${assembler}/abacas`	
* `{sample_id}`	
  * {samples_id}_abacas.fasta: Ordered and orientated sequence file.	
  * {sample_id}_abacas.tab: Feature file.	
  * {sample_id}_abacas.bin: Bin file that contains contigs that are not used in ordering.	
  * {sample_id}_abacas.crunch: Comparison file.	
  * {sample_id}_abacas.gaps: Gap information.	
  * unused_contigs.out: Information on contigs that have a mapping information but could not be used in the ordering.	
  * {samples_id}_abacas.MULTIFASTA.fa: A list of ordered and orientated contigs in a multi-fasta format.

## MultiQC

[MultiQC](http://multiqc.info) is a visualisation tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in within the report data directory.

The pipeline has special steps which allow the software versions used to be reported in the MultiQC output for future traceability.

**Output directory: `results/multiqc`**

* `Project_multiqc_report.html`
  * MultiQC report - a standalone HTML file that can be viewed in your web browser
* `Project_multiqc_data/`
  * Directory containing parsed statistics from the different tools used in the pipeline

For more information about how to use MultiQC reports, see [http://multiqc.info](http://multiqc.info)
