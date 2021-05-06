# nf-core/viralrecon: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [[2.0dev](https://github.com/nf-core/rnaseq/releases/tag/2.0)] - 2021-XX-XX

### :warning: Major enhancements

* Pipeline has been re-implemented in [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html)
* All software containers are now exclusively obtained from [Biocontainers](https://biocontainers.pro/#/registry)
* Updated minimum Nextflow version to `v21.04.0` (see [nextflow#572](https://github.com/nextflow-io/nextflow/issues/1964))
* Default human `--kraken2_db` link has been changed from Zenodo to an AWS S3 bucket for more reliable downloads
* [BCFtools](http://samtools.github.io/bcftools/bcftools.html) and [iVar](https://github.com/andersen-lab/ivar) will be run by default for Illumina metagenomics and amplicon data, respectively. However, this behaviour can be customised with the `--callers` parameter.
* Illumina and Nanopore runs containing the same 48 samples sequenced on both platforms have been uploaded to the nf-core AWS account for full-sized tests on release
* Added [Pangolin](https://github.com/cov-lineages/pangolin) for lineage analysis
* Added [Nextclade](https://github.com/nextstrain/nextclade) for clade assignment, mutation calling and consensus sequence quality checks
* Variant graph processes to call variants relative to the reference genome directly from _de novo_ assemblies have been deprecated and removed
* Variant calling with Varscan 2 has been deprecated and removed due to [licensing restrictions](https://github.com/dkoboldt/varscan/issues/12)

### Other enhancements & fixes

* Updated pipeline template to nf-core/tools `1.13.3`
* Optimise MultiQC configuration and input files for faster run-time on huge sample numbers
* [#122](https://github.com/nf-core/viralrecon/issues/122) - Single SPAdes command to rule them all
* [#138](https://github.com/nf-core/viralrecon/issues/138) - Problem masking the consensus sequence
* [#142](https://github.com/nf-core/viralrecon/issues/142) - Unknown method invocation `toBytes` on String type
* [#169](https://github.com/nf-core/viralrecon/issues/169) - ggplot2 error when generating mosdepth amplicon plot with Swift v2 primers
* [#170](https://github.com/nf-core/viralrecon/issues/170) - ivar trimming of Swift libraries new offset feature
* [#175](https://github.com/nf-core/viralrecon/issues/175) - MultiQC report does not include all the metrics

### Parameters

| Old parameter                 | New parameter                         |
|-------------------------------|---------------------------------------|
| `--amplicon_bed`              | `--primer_bed`                        |
| `--amplicon_fasta`            | `--primer_fasta`                      |
| `--amplicon_left_suffix`      | `--primer_left_suffix`                |
| `--amplicon_right_suffix`     | `--primer_right_suffix`               |
| `--filter_dups`               | `--filter_duplicates`                 |
| `--skip_adapter_trimming`     | `--skip_fastp`                        |
| `--skip_amplicon_trimming`    | `--skip_cutadapt`                     |
|                               | `--artic_minion_aligner`              |
|                               | `--artic_minion_caller`               |
|                               | `--artic_minion_medaka_model`         |
|                               | `--blast_db`                          |
|                               | `--enable_conda`                      |
|                               | `--fast5_dir`                         |
|                               | `--fastq_dir`                         |
|                               | `--ivar_trim_offset`                  |
|                               | `--kraken2_assembly_host_filter`      |
|                               | `--kraken2_variants_host_filter`      |
|                               | `--min_barcode_reads`                 |
|                               | `--min_guppyplex_reads`               |
|                               | `--multiqc_title`                     |
|                               | `--platform`                          |
|                               | `--primer_set`                        |
|                               | `--primer_set_version`                |
|                               | `--public_data_ids`                   |
|                               | `--save_trimmed_fail`                 |
|                               | `--save_unaligned`                    |
|                               | `--sequencing_summary`                |
|                               | `--singularity_pull_docker_container` |
|                               | `--skip_bandage`                      |
|                               | `--skip_consensus`                    |
|                               | `--skip_ivar_trim`                    |
|                               | `--skip_nanoplot`                     |
|                               | `--skip_pangolin`                     |
|                               | `--skip_pycoqc`                       |
|                               | `--skip_nextclade`                    |
|                               | `--skip_sra_fastq_download`           |
|                               | `--spades_hmm`                        |
|                               | `--spades_mode`                       |
| `--cut_mean_quality`          |                                       |
| `--filter_unmapped`           |                                       |
| `--ivar_trim_min_len`         |                                       |
| `--ivar_trim_min_qual`        |                                       |
| `--ivar_trim_window_width`    |                                       |
| `--kraken2_use_ftp`           |                                       |
| `--max_allele_freq`           |                                       |
| `--min_allele_freq`           |                                       |
| `--min_base_qual`             |                                       |
| `--min_coverage`              |                                       |
| `--min_trim_length`           |                                       |
| `--minia_kmer`                |                                       |
| `--mpileup_depth`             |                                       |
| `--name`                      |                                       |
| `--qualified_quality_phred`   |                                       |
| `--save_align_intermeds`      |                                       |
| `--save_kraken2_fastq`        |                                       |
| `--save_sra_fastq`            |                                       |
| `--skip_sra`                  |                                       |
| `--skip_vg`                   |                                       |
| `--unqualified_percent_limit` |                                       |
| `--varscan2_strand_filter`    |                                       |

> **NB:** Parameter has been __updated__ if both old and new parameter information is present.  
> **NB:** Parameter has been __added__ if just the new parameter information is present.  
> **NB:** Parameter has been __removed__ if new parameter information isn't present.

### Software dependencies

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency                    | Old version | New version |
|-------------------------------|-------------|-------------|
| `artic`                       |             | 1.2.1       |
| `bc`                          | 1.07.1      |             |
| `bcftools`                    | 1.9         | 1.11        |
| `bedtools`                    | 2.29.2      | 2.30.0      |
| `bioconductor-biostrings`     | 2.54.0      | 2.58.0      |
| `bioconductor-complexheatmap` | 2.2.0       | 2.6.2       |
| `blast`                       | 2.9.0       | 2.10.1      |
| `bowtie2`                     | 2.4.1       | 2.4.2       |
| `cutadapt`                    | 2.10        | 3.2         |
| `ivar`                        | 1.2.2       | 1.3.1       |
| `kraken2`                     | 2.0.9beta   | 2.1.1       |
| `markdown`                    | 3.2.2       |             |
| `minimap2`                    | 2.17        |             |
| `mosdepth`                    | 0.2.6       | 0.3.1       |
| `multiqc`                     | 1.9         | 1.10.1      |
| `nanoplot`                    |             | 1.36.1      |
| `nextclade_js`                |             | 0.14.2      |
| `pangolin`                    |             | 2.4.2       |
| `parallel-fastq-dump`         | 0.6.6       |             |
| `picard`                      | 2.23.0      | 2.23.9      |
| `pigz`                        | 2.3.4       |             |
| `plasmidid`                   | 1.6.3       | 1.6.4       |
| `pycoqc`                      |             | 2.5.2       |
| `pygments`                    | 2.6.1       |             |
| `pymdown-extensions`          | 7.1         |             |
| `python`                      | 3.6.10      | 3.8.3       |
| `r-base`                      | 3.6.2       | 4.0.3       |
| `r-ggplot2`                   | 3.3.1       | 3.3.3       |
| `r-tidyr`                     | 1.1.0       |             |
| `requests`                    |             | 2.24.0      |
| `samtools`                    | 1.9         | 1.10        |
| `seqwish`                     | 0.4.1       |             |
| `snpeff`                      | 4.5covid19  | 5.0         |
| `spades`                      | 3.14.0      | 3.15.2      |
| `sra-tools`                   | 2.10.7      |             |
| `tabix`                       |             | 0.2.6       |
| `unicycler`                   | 0.4.7       | 0.4.8       |
| `varscan`                     | 2.4.4       |             |
| `vg`                          | 1.24.0      |             |

> **NB:** Dependency has been __updated__ if both old and new version information is present.  
> **NB:** Dependency has been __added__ if just the new version information is present.  
> **NB:** Dependency has been __removed__ if new version information isn't present.  

## [[1.1.0](https://github.com/nf-core/rnaseq/releases/tag/1.1.0)] - 2020-06-23

### Added

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

### Removed

* `--skip_qc` parameter

### Dependencies

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

## [[1.0.0](https://github.com/nf-core/rnaseq/releases/tag/1.0.0)] - 2020-06-01

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
