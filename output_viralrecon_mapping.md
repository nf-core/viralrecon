# Output description for SARS-Cov2 consensus genome pipeline
**viralrecon** is a bioinformatics best-practice analysis pipeline used for SARS-Cov2 viral genome reconstruction. These results are focused in mapping and variants calling.

## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

* <a href="#fastqc">FastQC</a> v0.11.8 - read quality control
* <a href="#trimming">Trimmomatic</a> v.0.38 - adapter and low quality trimming.
* <a href="#kraken">Kraken2</a> v.2.0.8_beta - Mapping to host genome.
* <a href="#bowtie">Bowtie2</a> v2.3.5 - mapping against reference genomes.
* <a href="#samtools">SAMtools</a> v1.2 - Mapping result processing and unmapped reads selection.
* <a href="#picard">Picard</a> v1.140 - Enrichment and alignment metrics.
* <a href="#varscan">VarScan</a> v2.4.4-0 - Variant calling.
* <a href="#snpeff">snpEff and SnpSift</a> 4.3.1t - Variant calling annotation.
* <a href="#bcftools">Bcftools</a> v1.9 - Variant calling index and consensus genome generation.
* <a href="#bedtools">Bedtools</a> v2.26.0 - Consensus genome masking.
* <a href="#muscle">Muscle</a> v3.8.1551 - Multipel sequence alignment.
* <a href="#multiqc">MultiQC</a> v1.7 - aggregate report, describing results of the whole pipeline.

**Reference genome:** NC_045512.2

Depending on the analysis, we will have some ANALYSIS_IDs. This ANALYSIS_IDs are going to be composed of the date of the analysis, and some identification of the type of analysis.

## Preprocessing
### FastQC
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)<a name="fastqc"> </a><a href="#fastqc_reference">[1]</a> gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

**Output directory: `preprocess/fastqc`**

* `{sample_id}_T1_[12]_fastqc.html`
  * html report. This file can be opened in your favorite web browser (Firefox/chrome preferable) and it contains the different graphs that fastqc calculates for QC.
* `zips/{sample_id}_T1_[12]_fastqc.zip`
  * zip file containing the FastQC report, tab-delimited data file and plot images

### Trimming
[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)<a name="trimming"> </a><a href="#trimming_reference">[2]</a> is used for removal of adapter contamination and trimming of low quality regions.
Parameters included for trimming are:
-  Nucleotides with phred quality < 10 in 3'end.
-  Mean phred quality < 20 in a 4 nucleotide window.
-  Read lenght < 50

**Results directory: `preprocess/trimmomatic`**
- Files:
   - `fastqc/{sample_id}_T1_[12].trimmed_fastqc.html`: FastQC report of the trimmed reads.
   - `logs/{sample_id}_T1.trimmomatic.log`: trimming log file.

## Mapping
### kraken2
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

## Stats
### Muscle
[Muscle](http://www.drive5.com/muscle/muscle.html)<a name="muscle"> </a><a href="#Muscle_reference">[11]</a> is a program for creating multiple alignments of amino acid or nucleotide sequences. A range of options is provided that give you the choice of optimizing accuracy, speed, or some compromise between the two. Default parameters are those that give the best average accuracy in our tests. Using versions current at the time of writing, my tests show that MUSCLE can achieve both better average accuracy and better speed than CLUSTALW or T‑Coffee, depending on the chosen options.

**Output directory: `variants/muscle`**
* `all_genomes_masked_reference.fasta`
  * Multifasta with the reference genome and the masked fasta file of the consensus.
* `all_genomes_masked_aligned.fasta`
  * Muscle alignment of the multifasta file.

### MultiQC
[MultiQC](http://multiqc.info)<a name="multiqc"> </a><a href="#MultiQC_reference">[12]</a> is a visualization tool that generates a single HTML report summarizing all samples in your project. Most of the pipeline QC results are visualized in the report and further statistics are available in within the report data directory.

The pipeline has special steps which allow the software versions used to be reported in the MultiQC output for future traceability.

**Output directory: `multiqc`**

* `multiqc_report.html`
  * MultiQC report: standalone HTML file that can be viewed in your web browser
* `multiqc_data/`
  * Directory containing parsed statistics from the different tools used in the pipeline
* `multiqc_plots/`
  * The different plots contained in the report in different image formats.

For more information about how to use MultiQC reports, see [http://multiqc.info](http://multiqc.info)

## References
1. <a name="fastqc_reference"></a> <a href="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/">web FASTQC</a> P. Babraham Bioinformatics - FastQC. A Quality Control tool for High Throughput Sequence Data. 2012
trimming_reference
* <a name="trimming_reference"></a> Bolger AM, Lohse M, Usadel B. <a href="https://www.ncbi.nlm.nih.gov/pubmed/24695404">Trimmomatic: a flexible trimmer for Illumina sequence data.</a> Bioinformatics. 2014 Apr 28;30(15):2114–20.
* <a name="Kraken_reference"></a> Wood D.E., Lu J., Langmead B. <a href="https://www.ncbi.nlm.nih.gov/pubmed/?term=Improved+metagenomic+analysis+with+Kraken+2"> Improved metagenomic analysis with Kraken 2.</a> Genome Biol. 2019 Nov 28;20(1):257.
* <a name="Bowtie_reference"></a> Langmead B., Salzberg S.L.<a href="https://www.ncbi.nlm.nih.gov/pubmed/?term=+Fast+gapped-read+alignment+with+Bowtie+2">Fast gapped-read alignment with Bowtie 2.</a> Nat Methods. 2012 Mar 4;9(4):357-9
* <a name="SAMtools_eference"></a> Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R.; 1000 Genome Project Data Processing Subgroup. <a href="https://www.ncbi.nlm.nih.gov/pubmed/19505943">The Sequence Alignment/Map format and SAMtools.</a> Bioinformatics. 2009 Aug 15;25(16):2078-9.
* <a name="picard_reference"></a> <a href="https://github.com/broadinstitute/picard">web Picard toolkit</a> Broad Institute.
* <a name="VarScan_reference"></a> Koboldt D.C., Chen K., Wylie T., Larson D.E., McLellan M.D., Mardis E.R., Weinstock G.M., Wilson R.K., Ding L.<a href="https://www.ncbi.nlm.nih.gov/pubmed/19542151">VarScan: variant detection in massively parallel sequencing of individual and pooled samples.</a> Bioinformatics. 2009 Sep 1;25(17):2283-5.
* <a name="SnpEff_reference"></a> Cingolani P., Platts A., Wang le L., Coon M., Nguyen T., Wang L., Land S.J., Lu X., Ruden D.M. <a href="https://www.ncbi.nlm.nih.gov/pubmed/?term=A+program+for+annotating+and+predicting+the+effects+of+single+nucleotide+polymorphisms%2C+SnpEff%3A+SNPs+in+the+genome+of+Drosophila+melanogaster+strain+w1118">A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3.</a>Fly (Austin). 2012 Apr-Jun;6(2):80-92
* <a name="SnpSift_reference"></a>Cingolani P., Patel V.M., Coon M., Nguyen T., Land S.J., Ruden D.M., Lu X. <a href="https://www.ncbi.nlm.nih.gov/pubmed/?term=Using+Drosophila+melanogaster+as+a+model+for+genotoxic+chemical+mutational+studies+with+a+new+program%2C+SnpSift">Using Drosophila melanogaster as a Model for Genotoxic Chemical Mutational Studies with a New Program, SnpSift.</a>Front Genet. 2012 Mar 15;3:35.
* <a name="Bedtools_reference"></a> Quinlan A.R., Hall I.M.<a href="https://www.ncbi.nlm.nih.gov/pubmed/?term=BEDTools%3A+a+flexible+suite+of+utilities+for+comparing+genomic+features">BEDTools: a flexible suite of utilities for comparing genomic features.</a>Bioinformatics. 2010 Mar 15;26(6):841-2.
* <a name="Muscle_reference"></a> Edgar R.C.<a href="https://www.ncbi.nlm.nih.gov/pubmed/?term=MUSCLE%3A+a+multiple+sequence+alignment+method+with+reduced+time+and+space+complexity"> MUSCLE: a multiple sequence alignment method with reduced time and space complexity.</a> BMC Bioinformatics. 2004 Aug 19;5:113.
* <a name="MultiQC_reference"></a> Ewels P, Magnusson M, Lundin S, Käller M. <a href="https://www.ncbi.nlm.nih.gov/pubmed/?term=MultiQC%3A+Summarize+analysis+results+for+multiple+tools+and+samples+in+a+single+report">MultiQC: summarize analysis results for multiple tools and samples in a single report.</a> Bioinformatics. 2016 Oct 1;32(19):3047-8.
