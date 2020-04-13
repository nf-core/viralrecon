# nf-core/viralrecon: Output

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [FastQC](#fastqc) - read quality control
* [Kraken2](#kraken2) - Mapping to host genome.
* [bowtie2](#bowtie2) - mapping against reference genomes.
* [SAMtools](#samtools) - Mapping result processing and unmapped reads selection.
* [Picard](#picard) - Enrichment and alignment metrics.
* [VarScan](#varscan) - Variant calling.
* [SnpEff and SnpSift] - Variant calling annotation.
* [Bcftools](#bcftools) - Variant calling index and consensus genome generation.
* [Bedtools](#bedtools) - Consensus genome masking.
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

### Trimming

#### Quality trimming

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
