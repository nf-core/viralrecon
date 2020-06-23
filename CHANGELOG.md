# nf-core/viralrecon: Changelog

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2020-06-23

### `Added`

* [#112](https://github.com/nf-core/viralrecon/issues/112) - Per-amplicon coverage plot
* [#124](https://github.com/nf-core/viralrecon/issues/124) - Intersect variants across callers
* [nf-core/tools#616](https://github.com/nf-core/tools/pull/616) - Updated GitHub Actions to build Docker image and push to Docker Hub
* Parameters:
    * `--min_mapped_reads` to circumvent failures for samples with low number of mapped reads
    * `--varscan2_strand_filter` to toggle the default Varscan 2 strand filter
    * `--skip_mosdepth` - skip genome-wide and amplicon coverage plot generation from mosdepth output
    * `--amplicon_left_suffix` - to provide left primer suffix used in name field of `--amplicon_bed`
    * `--amplicon_right_suffix` - to provide right primer suffix used in name field of `--amplicon_bed`
    * Unify parameter specification with COG-UK pipeline:
        * `--min_allele_freq` - minimum allele frequency threshold for calling variants
        * `--mpileup_depth` - SAMTools mpileup max per-file depth
        * `--ivar_exclude_reads` renamed to `--ivar_trim_noprimer`
        * `--ivar_trim_min_len` - minimum length of read to retain after primer trimming
        * `--ivar_trim_min_qual` - minimum quality threshold for sliding window to pass
        * `--ivar_trim_window_width` - width of sliding window
* [#118] Updated GitHub Actions AWS workflow for small and full size tests.

### `Removed`

* `--skip_qc` parameter

### `Dependencies`

* Add mosdepth `0.2.6`
* Add bioconductor-complexheatmap `2.2.0`
* Add bioconductor-biostrings `2.54.0`
* Add r-optparse `1.6.6`
* Add r-tidyr `1.1.0`
* Add r-tidyverse `1.3.0`
* Add r-ggplot2 `3.3.1`
* Add r-reshape2 `1.4.4`
* Add r-viridis `0.5.1`
* Update sra-tools `2.10.3` -> `2.10.7`
* Update bowtie2 `2.3.5.1` -> `2.4.1`
* Update picard `2.22.8` -> `2.23.0`
* Update minia `3.2.3` -> `3.2.4`
* Update plasmidid `1.5.2` -> `1.6.3`

## [1.0.0] - 2020-06-01

Initial release of nf-core/viralrecon, created with the [nf-core](http://nf-co.re/) template.

This pipeline is a re-implementation of the [SARS_Cov2_consensus-nf](https://github.com/BU-ISCIII/SARS_Cov2_consensus-nf) and [SARS_Cov2_assembly-nf](https://github.com/BU-ISCIII/SARS_Cov2_assembly-nf) pipelines initially developed by [Sarai Varona](https://github.com/svarona) and [Sara Monzon](https://github.com/saramonzon) from [BU-ISCIII](https://github.com/BU-ISCIII). Porting both of these pipelines to nf-core was an international collaboration between numerous contributors and developers, led by [Harshil Patel](https://github.com/drpatelh) from the [The Bioinformatics & Biostatistics Group](https://www.crick.ac.uk/research/science-technology-platforms/bioinformatics-and-biostatistics/) at [The Francis Crick Institute](https://www.crick.ac.uk/), London. We appreciated the need to have a portable, reproducible and scalable pipeline for the analysis of COVID-19 sequencing samples and so the Avengers Assembled!

### Pipeline summary

1. Download samples via SRA, ENA or GEO ids ([`ENA FTP`](https://ena-docs.readthedocs.io/en/latest/retrieval/file-download.html), [`parallel-fastq-dump`](https://github.com/rvalieris/parallel-fastq-dump); *if required*)
2. Merge re-sequenced FastQ files ([`cat`](http://www.linfo.org/cat.html); *if required*)
3. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
4. Adapter trimming ([`fastp`](https://github.com/OpenGene/fastp))
5. Variant calling
    1. Read alignment ([`Bowtie 2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml))
    2. Sort and index alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    3. Primer sequence removal ([`iVar`](https://github.com/andersen-lab/ivar); *amplicon data only*)
    4. Duplicate read marking ([`picard`](https://broadinstitute.github.io/picard/); *removal optional*)
    5. Alignment-level QC ([`picard`](https://broadinstitute.github.io/picard/), [`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
    6. Choice of multiple variant calling and consensus sequence generation routes ([`VarScan 2`](http://dkoboldt.github.io/varscan/), [`BCFTools`](http://samtools.github.io/bcftools/bcftools.html), [`BEDTools`](https://github.com/arq5x/bedtools2/) *||* [`iVar variants and consensus`](https://github.com/andersen-lab/ivar) *||* [`BCFTools`](http://samtools.github.io/bcftools/bcftools.html), [`BEDTools`](https://github.com/arq5x/bedtools2/))
        * Variant annotation ([`SnpEff`](http://snpeff.sourceforge.net/SnpEff.html), [`SnpSift`](http://snpeff.sourceforge.net/SnpSift.html))
        * Consensus assessment report ([`QUAST`](http://quast.sourceforge.net/quast))
6. _De novo_ assembly
    1. Primer trimming ([`Cutadapt`](https://cutadapt.readthedocs.io/en/stable/guide.html); *amplicon data only*)
    2. Removal of host reads ([`Kraken 2`](http://ccb.jhu.edu/software/kraken2/))
    3. Choice of multiple assembly tools ([`SPAdes`](http://cab.spbu.ru/software/spades/) *||* [`metaSPAdes`](http://cab.spbu.ru/software/meta-spades/) *||* [`Unicycler`](https://github.com/rrwick/Unicycler) *||* [`minia`](https://github.com/GATB/minia))
        * Blast to reference genome ([`blastn`](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch))
        * Contiguate assembly ([`ABACAS`](https://www.sanger.ac.uk/science/tools/pagit))
        * Assembly report ([`PlasmidID`](https://github.com/BU-ISCIII/plasmidID))
        * Assembly assessment report ([`QUAST`](http://quast.sourceforge.net/quast))
        * Call variants relative to reference ([`Minimap2`](https://github.com/lh3/minimap2), [`seqwish`](https://github.com/ekg/seqwish), [`vg`](https://github.com/vgteam/vg), [`Bandage`](https://github.com/rrwick/Bandage))
        * Variant annotation ([`SnpEff`](http://snpeff.sourceforge.net/SnpEff.html), [`SnpSift`](http://snpeff.sourceforge.net/SnpSift.html))
7. Present QC and visualisation for raw read, alignment, assembly and variant calling results ([`MultiQC`](http://multiqc.info/))
