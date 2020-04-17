# ![nf-core/viralrecon](images/nf-core-viralrecon_logo.png)

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the output directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

* [Preprocessing](#Preprocessing)
  * [parallel-fastq-dump](#parallel-fastq-dump) - Download samples from SRA
  * [FastQC](#fastqc) - Raw read QC
  * [fastp](#fastp) - Adapter and quality trimming
  * [cat](#cat) - Merge re-sequenced FastQ files
* [Variant calling](#variant-calling)
  * [Bowtie 2](#bowtie-2) - Read alignment relative to reference genome
  * [SAMtools](#samtools) - Sort, index and generate metrics for alignments
  * [iVar trim](#ivar=trim) - Primer sequence removal for amplicon data
  * [picard-tools](#picard-tools) - Whole genome coverage and alignment metrics
  * [VarScan 2, BCFTools, BEDTools](#varscan-2-bcftools-bedtools) - Option 1: Variant calling, consensus sequence generation and masking
  * [iVar variants and iVar consensus](#ivar-variants-and-ivar-consensus) - Option 2: Variant calling and consensus sequence generation
  * [SnpEff and SnpSift](#snpeff-and-snpsift) - Genetic variant annotation and functional effect prediction
  * [QUAST](#quast) - Consensus assessment report
* [De novo assembly](#de-novo-assembly)
  * [Cutadapt](#cutadapt) - Primer trimming for amplicon data
  * [Kraken2](#kraken2) - Removal of host reads
  * [SPAdes](#spades) - Option 1: Viral genome assembly
  * [metaSPAdes](#metaspades) - Option 2: Viral genome assembly
  * [Unicycler](#unicycler) - Option 3: Viral genome assembly
  * [BLAST](#blast) - Blast to reference assembly
  * [ABACAS](#abacas) - Order contigs according to reference genome
  * [PlasmidID](#plasmidid) - Assembly report and visualisation
  * [QUAST](#quast) - Assembly quality assessment
* [Workflow reporting and genomes](#workflow-reporting-and-genomes)
  * [MultiQC](#multiqc) - Present QC for raw reads, alignment, assembly and variant calling
  * [Reference genome files](#reference-genome-files) - Saving reference genome indices/files
  * [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

## Preprocessing

### parallel-fastq-dump

TODO: Add some notes here about how we are downloading the data.

### FastQC

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality. To see how your reads look after trimming, look at the FastQC reports in the `fastp` directory.

**Output directory: `preprocess/fastqc/`**

* `preprocess/fastqc/`
  * `<SAMPLE>_fastqc.html`: FastQC report, containing quality metrics for your untrimmed raw fastq files.
* `preprocess/fastqc/zips/`
  * `<SAMPLE>_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

![FastQC per base sequence plot](images/fastqc_per_base_sequence_quality_plot-1.png)

### fastp

[fastp](https://github.com/OpenGene/fastp) is a tool designed to provide fast all-in-one preprocessing for FastQ files. It is developed in C++ with multithreading support to afford high performance. We use fastp for adapter trimming and quality filtering.

**Output directory: `preprocess/fastp/`**


* `<SAMPLE>.trim.fastq.gz`: Paired trimmed reads. Only present if `--save_trimmed` parameter is supplied.
* `<SAMPLE>.trim.fail.gz`: Unpaired trimmed reads. Only present if `--save_trimmed` parameter is supplied.
* `<SAMPLE>.fastp.html`: fastp report in html format.
* `<SAMPLE>.fastp.json`: fastp report in json format.
* `fastqc/<SAMPLE>.trim.fastqc.html`: FastQC report of the trimmed reads.
* `fastqc/zips/<SAMPLE>.trim.fastqc.zip`: Zip archive containing the FastQC report.
* `log/<SAMPLE>.fastp.log`: Trimming log file.

![fastp filtered reads plot](images/fastp_filtered_reads_plot-1.png)

### cat

**Output directory: `preprocess/merged_fastq`**

TODO: Add some notes here about why and how we are merging FastQ files from the same sample.

## Variant calling

### Bowtie 2

[Bowtie 2](http://bio-bwa.sourceforge.net/) is an ultrafast and memory-efficient tool for aligning sequencing reads to long reference sequences. Bowtie 2 supports gapped, local, and paired-end alignment modes.

**Output directory: `variants/bowtie2`**

* `<SAMPLE>.bam`: Original BAM file containing mapped reads. Only present if `--save_align_intermeds` parameter is supplied.
* `log/<SAMPLE>.log`: Bowtie2 mapping log file.
* `<SAMPLE>.sorted.bam`: Sorted aligned BAM file.
* `<SAMPLE>.sorted.bam.bai`: Index file for sorted aligned BAM file.

![Bowtie2 paired end reads quality score plot](images/bowtie2_pe_plot-1.png)

### SAMtools

The resulting BAM files are further processed with [SAMtools](http://samtools.sourceforge.net/) for co-ordinate sorting and indexing of the alignments. SAMtools is also used to generate read mapping statistics.

**Output directory: `variants/bowtie2/samtools_stats`**

* `<SAMPLE>.sorted.bam.flagstat`: Samtools flagstats mapping stats summary.
* `<SAMPLE>.sorted.bam.idxstats`: Samtools stats in the mapping index file.
* `<SAMPLE>.sorted.bam.stats`: Samtools mapping stats report.

![SAMtools alignment quality scores plot](images/samtools_alignment_plot-1.png)

### iVar trim

If `--protocol amplicon` is set then [iVar](http://gensoft.pasteur.fr/docs/ivar/1.0/manualpage.html) is used to trim the amplicon primer sequences from the reads. iVar uses the primer positions supplied in `--amplicon_bed` to soft clip primer sequences from an aligned and sorted BAM file.

**Output directory: `variants/ivar`**

* `<SAMPLE>.trim.sorted.bam`: Sorted aligned BAM file after trimming.
* `<SAMPLE>.trim.sorted.bam.bai`: Index file for sorted aligned trimmed BAM.
* `log/<SAMPLE>.trim.ivar.log`: iVar log file.
* `samtools_stats/<SAMPLE>.trim.sorted.bam.flagstat`: SAMtools flagstats summary file.
* `samtools_stats/<SAMPLE>.trim.sorted.bam.idxstats`: SAMtools stats in the mapping index file.
* `samtools_stats/<SAMPLE>.trim.sorted.bam.stats`: SAMtools mapping stats report.
* `<SAMPLE>.trim.stats`: Picard metrics summary file for evaluating coverage and performance.

### picard-tools

[picard-tools](https://broadinstitute.github.io/picard/index.html) is a set of command-line tools for manipulating high-throughput sequencing data. In this case we use it to obtain mapping and coverage metrics. If `--protocol amplicon` is set then these metrics will be obtained from the IVar trimmed alignments as opposed to the original Bowtie 2 alignments.

**Output directory: `variants/<bowtie2/ivar>/picard_metrics`**

* `<SAMPLE>.CollectWgsMetrics.coverage_metrics`: Picard metrics summary file for evaluating coverage and performance.
* `<SAMPLE>.CollectMultipleMetrics.alignment_summary_metrics`: Summary metrics of the alignment.
* `<SAMPLE>.CollectMultipleMetrics.base_distribution_by_cycle.pdf`: PDF file with the Base percentage plotted against the number of cycle.
* `<SAMPLE>.CollectMultipleMetrics.base_distribution_by_cycle_metrics`: Metrics file used to plot `<SAMPLE>.CollectMultipleMetrics.base_distribution_by_cycle.pdf`.
* `<SAMPLE>.CollectMultipleMetrics.insert_size_histogram.pdf`: PDF file with the counts and the cumulative fraction of reads that is higher than the insert size against the insert size.
* `<SAMPLE>.CollectMultipleMetrics.insert_size_metrics`: Metrics about the insert size distribution of a paired-end library used to plot `<SAMPLE>.CollectMultipleMetrics.insert_size_histogram.pdf`.
* `<SAMPLE>.CollectMultipleMetrics.quality_by_cycle.pdf`: PDF file with the Mean Quality against the cycle number.
* `<SAMPLE>.CollectMultipleMetrics.quality_by_cycle_metrics`: Metrics file used to plot `<SAMPLE>.CollectMultipleMetrics.quality_by_cycle.pdf`.
* `<SAMPLE>.CollectMultipleMetrics.quality_distribution.pdf`: PDF file with the distribution of quality scores.
* `<SAMPLE>.CollectMultipleMetrics.quality_distribution_metrics`: Metrics file used to plot `<SAMPLE>.CollectMultipleMetrics.quality_distribution.pdf`.

Picard documentation: [Picard docs](https://broadinstitute.github.io/picard/command-line-overview.html)

![Picard insert size plot](images/picard_insert_size-1.png)

### VarScan 2, BCFTools, BEDTools

First of all SAMtools is used to generate the variant calling VCF file. Then [VarScan 2](http://varscan.sourceforge.net/) is used to call for major and low frequency variants. VarScan is a platform-independent software tool developed at the Genome Institute at Washington University to detect variants in NGS data.

[Bcftools](http://samtools.github.io/bcftools/bcftools.html) is a set of utilities that manipulate variant calls in the Variant Call Format (VCF) and its binary counterpart BCF. The resulting variant calling vcf for haploid genomes is indexed and then the consensus genome is created adding the variants to the reference viral genome. This consensus genome was obtained using the predominant variants (majority) of the mapping file.

[Bedtools](https://bedtools.readthedocs.io/en/latest/) are a swiss-army knife of tools for a wide-range of genomics analysis tasks. In this case we use:

1. bedtools genomecov computes histograms (default), per-base reports (-d) and BEDGRAPH (-bg) summaries of feature coverage (e.g., aligned sequences) for a given genome.
2. bedtools maskfasta masks sequences in a FASTA file based on intervals defined in a feature file. This may be useful fro creating your own masked genome file based on custom annotations or for masking all but your target regions when aligning sequence data from a targeted capture experiment.

**Output directory: `variants/varscan2/`**

* `<SAMPLE>.pileup`: If `--save_pileup`. Samtools pileup file. The pileup file summarizes all data from the reads at each genomic region that is covered by at least one read. Each row of the pileup file gives similar information to a single vertical column of reads in the IGV view.
* `<SAMPLE>.highfreq.vcf.gz`: High frequency variants VCF file.
* `<SAMPLE>.highfreq.vcf.gz.tbi`: High frequency variants VCF index file.
* `<SAMPLE>.lowfreq.vcf.gz`: Low frequency variants VCF file.
* `<SAMPLE>.lowfreq.vcf.gz.tbi`: Low frequency variants VCF index file.
* `log/<SAMPLE>.highfreq.varscan2.log`: VarScan2 high frequency variants log file.
* `log/<SAMPLE>.lowfreq.varscan2.log`: VarScan2 low frequency variants log file.

**Output directory: `variants/varscan2/consensus/`**

* `<SAMPLE>.consensus.fa`: Consensus viral genome file generated from adding the variants called before to the viral reference genome. These variants are only the majoritarian variants, including only SNPs and small indels.

  **Output directory: `variants/varscan2/consensus/`**

  * `<SAMPLE>.consensus.masked.fa`: Masked consensus fasta file.

### iVar variants and iVar consensus

iVar can also use the output of the samtools mpileup command to call variants - single nucleotide variants and indels.

**Output directory: `variants/ivar/variants/`**

* `<SAMPLE>.tsv`: TAB separated file with the variants.

Finally, iVar generates a consensus genome with the variants:

**Output directory: `variants/ivar/consensus/`**

* `<SAMPLE>.consensus.fa`: Fasta file with thte consensus genome.
* `<SAMPLE>.consensus.qual.txt`: File with the average quality of each base in the consensus sequence.

### SnpEff and SnpSift

[SnpEff](http://snpeff.sourceforge.net/SnpEff.html) is a genetic variant annotation and functional effect prediction toolbox. It annotates and predicts the effects of genetic variants on genes and proteins (such as amino acid changes).

[SnpSift](http://snpeff.sourceforge.net/SnpSift.html) annotates genomic variants using databases, filters, and manipulates genomic annotated variants. Once you have annotated your files using SnpEff, you can use SnpSift to help you filter large genomic datasets in order to find the most significant variants for your experiment.

**Output directory: `variants/varscan2/snpeff`**

* `<SAMPLE>.lowfreq.snpEff.csv`: Low frequency variants annotation csv file.
* `<SAMPLE>.lowfreq.snpSift.table.txt`: Low frequency variants SnpSift summary table.
* `<SAMPLE>.lowfreq.snpEff.vcf.gz`: Low frequency variants annotated VCF table.
* `<SAMPLE>.lowfreq.snpEff.genes.txt`: Low frequency variants genes table.
* `<SAMPLE>.lowfreq.snpEff.summary.html`: Low frequency variants summary html file.
* `<SAMPLE>.highfreq.snpEff.csv`: High frequency variants annotation csv file.
* `<SAMPLE>.highfreq.snpSift.table.txt`: High frequency variants SnpSift summary table.
* `<SAMPLE>.highfreq.snpEff.vcf.gz`: High frequency variants annotated VCF table.
* `<SAMPLE>.highfreq.snpEff.genes.txt`: High frequency variants genes table.
* `<SAMPLE>.highfreq.snpEff.summary.html`: High frequency variants summary html file.

### QUAST

TODO: Add description of what QUAST is doing here

## De novo assembly

### Cutadapt

When `--protocol amplicon` is set [Cutadapt](https://cutadapt.readthedocs.io/en/stable/guide.html) is used to clip primer sequences from reads prior to assembly.

**Output directory:** `assembly/cutadapt/`

* `fastqc/<SAMPLE>.ptrim_fastqc.html`: Fastqc HTML reports for primer clipped reads.
* `fastqc/zips/<SAMPLE>.ptrim_fastqc.zip`: Fastqc HTML reports for primer clipped reads.
* `log/<SAMPLE>.cutadapt.log`: Cutadapt logs.
* `<SAMPLE>.ptrim.fastq.gz`: FastQ files with primer sequences trimmed. Only present if `--save_trimmed` parameter is supplied.

### Kraken2

[Kraken2](https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual) is a sequence classifier that assigns taxonomic labels to DNA sequences. Kraken examines the k-mers within a query sequence and uses the information within those k-mers to query a database. That database maps k-mers to the lowest common ancestor (LCA) of all genomes known to contain a given k-mer.

We use a Kraken2 database in this workflow to filter out reads specific to the host genome. The remainder of the reads are then passed to numerous de novo assembly algorithms in order to reconstruct the viral genome assembly.

**Output directory: `assembly/kraken2`**

* `<SAMPLE>.host.fastq.gz`: Reads that were classified to the host database. Only present if `--save_kraken2_fastq` parameter is supplied.
* `<SAMPLE>.viral.fastq.gz`: Reads that were unclassified to the host database. Only present if `--save_kraken2_fastq` parameter is supplied.
* `<SAMPLE>.kraken2.report.txt`: Kraken taxonomic report. See [here](https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual#sample-report-output-format) for a detailed description of the format.

### SPAdes

[SPAdes](http://cab.spbu.ru/software/spades/) is a de Bruijn graph-based assembler. We selected the reads that didn't mapped with the host genome and assembled them using SPAdes to create a viral genome assembly.

**Output directory:** `assembly/spades`

* `<SAMPLE>.scaffolds.fasta`: SPAdes assembled scaffolds.

### metaSPAdes

[metaSPAdes](http://cab.spbu.ru/software/meta-spades/) is a de Bruijn graph-based assembler that is distributed with SPAdes (run via the `--meta` option). It can be used for the simultaneous reconstruction of multiple genomes as observed in metagenomics data.

**Output directory: `assembly/metaspades`**

* `<SAMPLE>.meta.scaffolds.fasta`: metaSPAdes assembled scaffolds.

### Unicycler

[Unicycler](https://github.com/rrwick/Unicycler) is an assembly pipeline that works as a SPAdes optimizer.

**Output directory: `assembly/unicycler`**

* `<SAMPLE>.assembly.fasta`: Assembled scaffolds.

### BLAST

[blastn](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch) is used to align the assembled contigs against the virus reference genome.

**Output directory: `assembly/<ASSEMBLER>/blast`**

* `<SAMPLE>.blast.txt`: BLAST results against the target virus.
* `<SAMPLE>.blast.filt.header.txt`: Filtered BLAST results.

### ABACAS

[ABACAS](https://www.sanger.ac.uk/science/tools/pagit) was developed to rapidly contiguate (align, order, orientate), visualize and design primers to close gaps on shotgun assembled contigs based on a reference sequence.

**Output directory: `assembly/<ASSEMBLER>/abacas/`**

* `<SAMPLE>/`
  * `<SAMPLE>_abacas.fasta`: Ordered and orientated sequence file.
  * `<SAMPLE>_abacas.tab`: Feature file.
  * `<SAMPLE>_abacas.bin`: Bin file that contains contigs that are not used in ordering.
  * `<SAMPLE>_abacas.crunch`: Comparison file.
  * `<SAMPLE>_abacas.gaps`: Gap information.
  * `unused_contigs.out`: Information on contigs that have a mapping information but could not be used in the ordering.
  * `<SAMPLE>_abacas.MULTIFASTA.fa`: A list of ordered and orientated contigs in a multi-fasta format.

### PlasmidID

[PlasmidID](https://github.com/BU-ISCIII/plasmidID) is used to graphically represent the alignment of the reference genome relative to a given assembly. This helps to visualize the coverage of the reference genome in the assembly. To find more information about the output files refer to the [documentation](https://github.com/BU-ISCIII/plasmidID/wiki/Understanding-the-image:-track-by-track)

**Output directory: `assembly/<ASSEMBLER>/plasmidid/**

* `<SAMPLE>/images/<SAMPLE>_<REF_VIR_NAME>.png`: PNG file with the visualization of the alignment between the assembled viral genome and the reference viral genome.
* `<SAMPLE>/data/`: Files used for drawing the circos images.
* `<SAMPLE>/database/`: Annotation files used for drawing the circos images.

### QUAST

[QUAST](http://bioinf.spbau.ru/quast) is used to evaluate the quality of assemblies across multiple samples. The HTML results can be opened within any browser (we recommend using Google Chrome). A single quast report will be generated to collate the results across all samples for each assembler.

**Output directory: `assembly/<ASSEMBLER>/quast/`**

* `report.html`
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

## Workflow reporting and genomes

### MultiQC

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarizing all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in within the report data directory.

Results generated by MultiQC collate pipeline QC from FastQC, TrimGalore, samtools flagstat, samtools idxstats, samtools stats, picard CollectMultipleMetrics, picard MarkDuplicates, Preseq, deepTools plotProfile, deepTools plotFingerprint and featureCounts. The default [`multiqc config file`](../assets/multiqc_config.yaml) also contains the provision for loading custom-content to report peak counts, FRiP scores, peak-to-gene annnotation proportions, sample-similarity heatmaps and PCA plots.  

The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

**Output directory: `results/multiqc/`**

* `multiqc/`  
  * `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  * `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  * `multiqc_plots/`: directory containing static images from the report in various formats.

### Reference genome files

A number of genome-specific files are generated by the pipeline because they are required for downstream processing of the results. If the `--save_reference` parameter is provided then the Bowtie 2 alignment indices, BLAST and Kraken2 databases generated by the pipeline will be saved in this directory. It is recommended to use the `--save_reference` parameter if you are using the pipeline to build a Kraken2 database for the host genome. This can be quite a time-consuming process and it permits their reuse for future runs of the pipeline or for other purposes.  

**Output directory: `results/genome/`**

* `genome/`  
  * `Bowtie2Index/`: a standalone HTML file that can be viewed in your web browser.
  * `BlastDB/`: directory containing parsed statistics from the different tools used in the pipeline.
  * `kraken2_<kraken2_db_name>/`: directory containing static images from the report in various formats.

### Pipeline information

[Nextflow!](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

*Output directories*:

* `pipeline_info/`  
  * Reports generated by the pipeline - `pipeline_report.html`, `pipeline_report.txt` and `software_versions.csv`.
  * Reports generated by Nextflow - `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.svg`.
  * Reformatted samplesheet files used as input to the pipeline - `samplesheet.pass.csv`.
* `Documentation/`  
  * Documentation for interpretation of results in HTML format - `results_description.html`.
