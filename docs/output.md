# nf-core/viralrecon: Output

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [FastQC](#fastqc) - read quality control
* [Fastp](#fastp) - read quality trimming
* [Kraken2](#kraken2) - Mapping to host genome.
* Mapping + variant calling + consensus
  * [bowtie2](#bowtie2) - mapping against reference genomes.
  * [SAMtools](#samtools) - Mapping result processing and unmapped reads selection.
  * [Picard](#picard) - Enrichment and alignment metrics.
  * [VarScan](#varscan) - Variant calling.
  * [SnpEff and SnpSift](#snpeff-and-snpsift) - Variant calling annotation.
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

* `<SAMPLE>_fastqc.html`
  * FastQC report, containing quality metrics for your untrimmed raw fastq files
* `zips/<SAMPLE>_fastqc.zip`
  * zip file containing the FastQC report, tab-delimited data file and plot images

### Fastp

[Fastp](https://github.com/OpenGene/fastp) is a tool designed to provide fast all-in-one preprocessing for FastQ files. This tool is developed in C++ with multithreading supported to afford high performance. Fastp is used for quality filtering and adapter trimming.

**Output directory: `preprocess/fastp`**

* `<SAMPLE>.trim.fastq.gz`
  * Only with `--save_trimmed` param. Paired trimmed reads.
* `<SAMPLE>.trim.fail.gz`
  * Only with `--save_trimmed` param. Unpaired trimmed reads.
* `<SAMPLE>.fastp.html`
  * Fastp report in html format.
* `<SAMPLE>.fastp.json`
  * Fastp report in json format.
* `fastqc/<SAMPLE>.trim.fastqc.html`
  * FastQC report of the trimmed reads.
* `fastqc/zips/<SAMPLE>.trim.fastqc.zip`
  * Zip file with the FastQC report.
* `logs/<SAMPLE>.fastp.log`
  * Trimming log file.

## Mapping + variant calling + consensus

### kraken2

[Kraken2](https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual) is a taxonomic sequence classifier that assigns taxonomic labels to DNA sequences. Kraken examines the k-mers within a query sequence and uses the information within those k-mers to query a database. That database maps k-mers to the lowest common ancestor (LCA) of all genomes known to contain a given k-mer.

We mapped the fastq file against the reference host genome.

**Output directory: `preprocess/kraken2`**

* `<SAMPLE>.host.fastq.gz`
  * Only with `--save_kraken2_fastq`. Reads that mapped with host taxon.
* `<SAMPLE>.viral.fastq.gz`
  * Only with `--save_kraken2_fastq`. Reads that mapped with viral taxon.
* `<SAMPLE>.kraken2.report.txt`
  * Kraken taxonomic report. The report contains one line per taxonomic classification nd the following columns:
    1. Percentage of fragments covered by the clade rooted at this taxon
    2. Number of fragments covered by the clade rooted at this taxon
    3. Number of fragments assigned directly to this taxon
    4. A rank code, indicating (U)nclassified, (R)oot, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, or (S)pecies. Taxa that are not at any of these 10 ranks have a rank code that is formed by using the rank code of the closest ancestor rank with a number indicating the distance from that rank.  E.g., "G2" is a rank code indicating a taxon is between genus and species and the grandparent taxon is at the genus rank.
    5. NCBI taxonomic ID number
    6. Indented scientific name

### Bowtie2

[Bowtie](http://bio-bwa.sourceforge.net/) is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences. It is particularly good at aligning reads of about 50 up to 100s of characters to relatively long genomes. Bowtie 2 indexes the genome with an FM Index (based on the Burrows-Wheeler Transform or BWT) to keep its memory footprint small. Bowtie 2 supports gapped, local, and paired-end alignment modes.

**Output directory: `variants/bowtie2`**

* `<SAMPLE>.bam`
  * Only if `--save_align_intermeds`. Mapping BAM file
* logs/`<SAMPLE>.log`
  * Bowtie2 mapping log file.
* `<SAMPLE>.sorted.bam`
  * Sorted aligned BAM file.
* `<SAMPLE>.sorted.bam.bai`
  * Index file for soreted aligned BAM file.

### SAMtools

The result mapping files are further processed with [SAMtools](http://samtools.sourceforge.net/), sam format is converted to bam, sorted and an index .bai is generated. Samtools is also used to generate statistics about the mapping process.

**Output directory: `variants/bowtie2/samtools_stats`**

* `<SAMPLE>.sorted.bam.flagstat`
  * Samtools flagstats mapping stats summary.
* `<SAMPLE>.sorted.bam.idxstats`
  * Samtools stats in the mapping index file.
* `<SAMPLE>.sorted.bam.stats`
  * Samtools mapping stats report.

### Picard

[Picard](https://broadinstitute.github.io/picard/index.html) is a set of command line tools for manipulating high-throughput sequencing (HTS) data. In this case we used it to obtain mapping stats. If we are running `--protocol amplicon` Picards is going to be run over iVar files. Else, it will be run over Bowtie+Samtools files.

**Output directory: `variants/[bowtie2/iVar]/picard_metrics`**

* `<SAMPLE>.CollectWgsMetrics.coverage_metrics`
  * Picard metrics summary file for evaluating coverage and performance.
* `<SAMPLE>.CollectMultipleMetrics.alignment_summary_metrics`
  * Summary metrics of the alignment.
* `<SAMPLE>.CollectMultipleMetrics.base_distribution_by_cycle.pdf`
  * PDF file with the Base percentage plotted against the number of cycle.
* `<SAMPLE>.CollectMultipleMetrics.base_distribution_by_cycle_metrics`
  * Metrics file used to plot `<SAMPLE>.CollectMultipleMetrics.base_distribution_by_cycle.pdf`.
* `<SAMPLE>.CollectMultipleMetrics.insert_size_histogram.pdf`
  * PDF file with the counts and the cumulative fraction of reads that is higher than the insert size against the insert size.
* `<SAMPLE>.CollectMultipleMetrics.insert_size_metrics`
  * Metrics about the insert size distribution of a paired-end library used to plot `<SAMPLE>.CollectMultipleMetrics.insert_size_histogram.pdf`.
* `<SAMPLE>.CollectMultipleMetrics.quality_by_cycle.pdf`
  * PDF file with the Mean Quality against the cycle number.
* `<SAMPLE>.CollectMultipleMetrics.quality_by_cycle_metrics`
  * Metrics file used to plot `<SAMPLE>.CollectMultipleMetrics.quality_by_cycle.pdf`.
* `<SAMPLE>.CollectMultipleMetrics.quality_distribution.pdf`
  * PDF file with the distribution of quality scores.
* `<SAMPLE>.CollectMultipleMetrics.quality_distribution_metrics`
  * Metrics file used to plot `<SAMPLE>.CollectMultipleMetrics.quality_distribution.pdf`.

Picard documentation: [Picarddocs](https://broadinstitute.github.io/picard/command-line-overview.html)

### iVar

If we are running the `--protocol amplicon`, [iVar](http://gensoft.pasteur.fr/docs/ivar/1.0/manualpage.html) is used to trim the amplicon primers. iVar uses primer positions supplied in a BED file to soft clip primer sequences from an aligned and sorted BAM file.

**Output directory: `variants/ivar`**

* `<SAMPLE>.trim.sorted.bam`
  * Sorted aligned bam file after trimming.
* `<SAMPLE>.trim.sorted.bam.bai`
  * Index file for sorted aligned trimmed bam.
* `log/<SAMPLE>.trim.ivar.log`
  * iVar log file.
* `samtools_stats/<SAMPLE>.trim.sorted.bam.flagstat`
  * Samtools flagstats summary file.
* `samtools_stats/<SAMPLE>.trim.sorted.bam.idxstats`
  * Samtools stats in the mapping index file.
* `samtools_stats/<SAMPLE>.trim.sorted.bam.stats`
  * Samtools mapping stats report.

iVar can also use the output of the samtools mpileup command to call variants - single nucleotide variants(SNVs) and indels.

**Output directory: `variants/ivar/variants`**

* `<SAMPLE>.tsv`
  * TAB separated file with the variants.
* `<SAMPLE>.vcf.gz`
  * VCF file with the iVar variants in the TSV file.
* `<SAMPLE>.vcf.gz.tbi`
  * VCF file index.
* `<SAMPLE>.bcftools_stats.txt`
  * Bcftools stats on the iVar variant calling.

Finally,iVar generates a consensus genome with the variants:

**Output directory: `variants/ivar/consensus`**

* `<SAMPLE>.consensus.fa`
  * Fasta file with thte consensus gneome.
* `<SAMPLE>.consensus.qual.txt`
  * TXT file with the average quality of each base in the consensus sequence.

### VarScan

First of all SAMtools is used to generate the variant calling VCF file. Then [VarScan](http://varscan.sourceforge.net/) is used to call for major and low frequency variants. VarScan is a platform-independent software tool developed at the Genome Institute at Washington University to detect variants in NGS data.

**Output directory: `variants/varscan2`**

* `<SAMPLE>.pileup`
  * If `--save_pileup`. Samtools pileup file. The pileup file summarizes all data from the reads at each genomic region that is covered by at least one read. Each row of the pileup file gives similar information to a single vertical column of reads in the IGV view.
* `<SAMPLE>.highfreq.vcf.gz`
  * High frequency variants VCF file.
* `<SAMPLE>.highfreq.vcf.gz.tbi`
  * High frequency variants VCF index file.
* `<SAMPLE>.lowfreq.vcf.gz`
  * Low frequency variants VCF file.
* `<SAMPLE>.lowfreq.vcf.gz.tbi`
  * Low frequency variants VCF index file.
* `logs/<SAMPLE>.highfreq.varscan2.log`
  * VarScan2 high frequency variants log file.
* `logs/<SAMPLE>.lowfreq.varscan2.log`
  * VarScan2 low frequency variants log file.

### SnpEff and SnpSift

[SnpEff](http://snpeff.sourceforge.net/SnpEff.html) is a genetic variant annotation and functional effect prediction toolbox. It annotates and predicts the effects of genetic variants on genes and proteins (such as amino acid changes).

[SnpSift](http://snpeff.sourceforge.net/SnpSift.html) annotates genomic variants using databases, filters, and manipulates genomic annotated variants. Once you annotated your files using SnpEff, you can use SnpSift to help you filter large genomic datasets in order to find the most significant variants for your experiment.

SnpEff and SnpSift are run for VarScan and iVar variant callings. The ouput will be in the corresponding variant calling folder. In the case of VarScan `lowfreq` files will correspond to the annotation of low frequency variants and `highfreq` will correspond to the annotation of high frequency variants.

**Output directory: `variants/[varscan2/ivar]/snpeff`**

* `<SAMPLE>.snpEff.csv`
  * Low frequency variants annotation csv file.
* `<SAMPLE>.snpSift.table.txt`
  * Low frequency variants SnpSift summary table.
* `<SAMPLE>.snpEff.vcf.gz`
  * Low frequency variants annotated VCF table.
* `<SAMPLE>.snpEff.genes.txt`
  * Low frequency variants genes table.
* `<SAMPLE>.snpEff.summary.html`
  * Low frequency variants summary html file.

### BCFtools

[Bcftools](http://samtools.github.io/bcftools/bcftools.html) is a set of utilities that manipulate variant calls in the Variant Call Format (VCF) and its binary counterpart BCF. The resulting variant calling vcf for haploid genomes is indexed and then the consensus genome is created adding the variants to the reference viral genome. This consensus genome was obtained using the predominant variants (majority) of the mapping file.

**Output directory: `variants/bcftools`**

* `<SAMPLE>.consensus.fa`
  * Consensus viral genome file generated from adding the variants called before to the viral reference genome. These variants are only the majoritarian variants, including only SNPs and small indels.

### Bedtools

[Bedtools](https://bedtools.readthedocs.io/en/latest/) are a swiss-army knife of tools for a wide-range of genomics analysis tasks. In this case we use:

* bedtools genomecov computes histograms (default), per-base reports (-d) and BEDGRAPH (-bg) summaries of feature coverage (e.g., aligned sequences) for a given genome.
* bedtools maskfasta masks sequences in a FASTA file based on intervals defined in a feature file. This may be useful fro creating your own masked genome file based on custom annotations or for masking all but your target regions when aligning sequence data from a targeted capture experiment.

  **Output directory: `variants/bcftools`**

* `<SAMPLE>.consensus.masked.fa`
  * Masked consensus fasta file.

## De novo assembly

We selected the reads that didn't cluster using kraken2 with the host genome and assembled them to create a viral genome assembly.

### Cutadapt

Only when running `--protocol amplicon`, [Cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html) is used for clipping primer sequences from reads prior to assembly.

**Output directory:** `preprocess/cutadapt`

* `fastqc/<SAMPLE>.ptrim_fastqc.html`
  * Fastqc HTML reports for primer clipped reads.
* `fastqc/zips/<SAMPLE>.ptrim_fastqc.zip`
  * Fastqc HTML reports for primer clipped reads.
* `log/<SAMPLE>.cutadapt.log`
  * Cutadapt logs.
* `<SAMPLE>.ptrim.fastq.gz`
  * Only if `--save_trimmed`. Fastq files with primer sequences trimmed.

### SPAdes

[SPAdes](https://kbase.us/applist/apps/kb_SPAdes/run_SPAdes/release?gclid=Cj0KCQiAt_PuBRDcARIsAMNlBdroQS7y2hPFuhagq1QPvQ39FcvGxbhtZwhn8YbxIB4LrGIHKjJ-iPwaAn_lEALw_wcB) is a de Bruijn graph-based assembler. We selected the reads that didn't mapped with the host genome and assembled them using SPAdes to create a viral genome assembly.

**Output directory:** `assembly/spades`

* `<SAMPLE>.scaffolds.fasta`
  * SPAdes sssembled scaffolds.

### MetaSPAdes

[MetaSPAdes](https://kbase.us/applist/apps/kb_SPAdes/run_SPAdes/release?gclid=Cj0KCQiAt_PuBRDcARIsAMNlBdroQS7y2hPFuhagq1QPvQ39FcvGxbhtZwhn8YbxIB4LrGIHKjJ-iPwaAn_lEALw_wcB) is a de Bruijn graph-based assembler, with the option `--meta` the assembler works for metagenomics date trying to reconstruct different genomes.

**Output directory: `assembly/metaspades`**

* `<SAMPLE>.meta.scaffolds.fasta`
  * MetaSPaded assembled scaffolds.

### Unicycler

[Unicycler](https://github.com/rrwick/Unicycler) is an assembly pipeline that works as a spades optimizer.

**Output directory: `assembly/unicycler`**

* `<SAMPLE>.assembly.fasta`
  * Assembled scaffolds.

### QUAST

[QUAST](http://bioinf.spbau.ru/quast) evaluates genome assemblies. We compared the reference genome with the contigs and scaffold assemblies. The HTML results can be opened with any browser (we recommend using Google Chrome). We have a quast folder for each assembler selected.

Quast is going to be run for all the different assemblers (SPAdes, MetaSPAdes and Unicycler) and for the different genome consensus generation programs (BCFtools and iVar).

**Output directory: `[assembly/variants]/[spades/metaspades/unicycler/ivar/varscan2]/quast`**

* `report.html`
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

**Output directory: `assembly/<ASSEMBLER>/blast`**

* `<SAMPLE>.blast.txt`
  * Blast results against the target virus.
* `<SAMPLE>.blast.filt.header.txt`
  * Filtered Blast results.

### PlasmidID

[PlasmidID](https://github.com/BU-ISCIII/plasmidID) was used to graphically represent the alignment of the reference genome with the assembly obtained with SPAdes. This helps to visualize the coverage of the reference genome in the assembly. To find more information about the output files go to the [plasmidID documentation](https://github.com/BU-ISCIII/plasmidID/wiki/Understanding-the-image:-track-by-track)

**Output directory: `assembly/<ASSEMBLER>/plasmidid**

* `<SAMPLE>/images/<SAMPLE>_<REF_VIR_NAME>.png`
  * PNG file with the visualization of the alignment between the assembled viral genome and the reference viral genome.
* `<SAMPLE>/data/`
  * Files used for drawing the circos images.
* `<SAMPLE>/database`
  * Annotation files used for drawing the circos images.

### ABACAS

[Abacas](abacas.sourceforge.ne) intended to rapidly contiguate (align, order, orientate), visualize and design primers to close gaps on shotgun assembled contigs based on a reference sequence.

**Output directory: `assembly/<ASSEMBLER>/abacas`**

* `<SAMPLE>`
  * `<SAMPLE>_abacas.fasta`: Ordered and orientated sequence file.
  * `<SAMPLE>_abacas.tab`: Feature file.
  * `<SAMPLE>_abacas.bin`: Bin file that contains contigs that are not used in ordering.
  * `<SAMPLE>_abacas.crunch`: Comparison file.
  * `<SAMPLE>_abacas.gaps`: Gap information.
  * `unused_contigs.out`: Information on contigs that have a mapping information but could not be used in the ordering.
  * `<SAMPLE>_abacas.MULTIFASTA.fa`: A list of ordered and orientated contigs in a multi-fasta format.

## MultiQC

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarizing all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in within the report data directory.

The pipeline has special steps which allow the software versions used to be reported in the MultiQC output for future traceability.

**Output directory: `results/multiqc`**

* `Project_multiqc_report.html`
  * MultiQC report - a standalone HTML file that can be viewed in your web browser
* `Project_multiqc_data/`
  * Directory containing parsed statistics from the different tools used in the pipeline

For more information about how to use MultiQC reports, see [http://multiqc.info](http://multiqc.info)
