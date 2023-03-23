# nf-core/viralrecon: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [[2.6.0](https://github.com/nf-core/viralrecon/releases/tag/2.6.0)] - 2023-03-23

### Credits

Special thanks to the following for their code contributions to the release:

- [Friederike Hanssen](https://github.com/FriederikeHanssen)
- [Hugo Tavares](https://github.com/tavareshugo)
- [James Fellows Yates](https://github.com/jfy133)
- [Jessica Wu](https://github.com/wutron)
- [Matthew Wells](https://github.com/mattheww95)
- [Maxime Garcia](https://github.com/maxulysse)
- [Phil Ewels](https://github.com/ewels)
- [Sara Monzón](https://github.com/saramonzon)

Thank you to everyone else that has contributed by reporting bugs, enhancements or in any other way, shape or form.

### Enhancements & fixes

- [[#297](https://github.com/nf-core/viralrecon/issues/297)] - Add tube map for pipeline
- [[#316](https://github.com/nf-core/viralrecon/issues/316)] - Variant calling isn't run when using `--skip_asciigenome` with metagenomic data
- [[#317](https://github.com/nf-core/viralrecon/issues/317)] - `ivar_variants_to_vcf`: Ignore lines without annotation in ivar tsv file
- [[#320](https://github.com/nf-core/viralrecon/issues/320)] - Pipeline fails at email step: Failed to invoke `workflow.onComplete` event handler
- [[#321](https://github.com/nf-core/viralrecon/issues/321)] - `ivar_variants_to_vcf` script: Duplicated positions in tsv file due to overlapping annotations
- [[#334](https://github.com/nf-core/viralrecon/issues/334)] - Longshot thread 'main' panicked at 'assertion failed: p <= 0.0' error
- [[#341](https://github.com/nf-core/viralrecon/issues/341)] - `artic/minion` and `artic/guppyplex`: Update module version 1.2.2 -> 1.2.3
- [[#348](https://github.com/nf-core/viralrecon/issues/348)] - Document full parameters of iVar consensus
- [[#349](https://github.com/nf-core/viralrecon/issues/349)] - ERROR in Script plasmidID
- [[#356](https://github.com/nf-core/viralrecon/issues/356)] - Add NEB SARS-CoV-2 primers
- [[#368](https://github.com/nf-core/viralrecon/issues/368)] - Incorrect depth from ivar variants reported in variants long table
- Updated pipeline template to [nf-core/tools 2.7.2](https://github.com/nf-core/tools/releases/tag/2.7.2)
- Add `tower.yml` for Report rendering in Nextflow Tower
- Use `--skip_plasmidid` by default

### Parameters

| Old parameter | New parameter |
| ------------- | ------------- |
| `--tracedir`  |               |

> **NB:** Parameter has been **updated** if both old and new parameter information is present.
> **NB:** Parameter has been **added** if just the new parameter information is present.
> **NB:** Parameter has been **removed** if new parameter information isn't present.

### Software dependencies

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency  | Old version | New version |
| ----------- | ----------- | ----------- |
| `artic`     | 1.2.2       | 1.2.3       |
| `bcftools`  | 1.51.1      | 1.16        |
| `blast`     | 2.12.0      | 2.13.0      |
| `cutadapt`  | 3.5         | 4.2         |
| `ivar`      | 1.3.1       | 1.4         |
| `multiqc`   | 1.13a       | 1.14        |
| `nanoplot`  | 1.40.0      | 1.41.0      |
| `nextclade` | 2.2.0       | 2.12.0      |
| `pangolin`  | 4.1.1       | 4.2         |
| `picard`    | 2.27.4      | 3.0.0       |
| `samtools`  | 1.15.1      | 1.16.1      |
| `spades`    | 3.15.4      | 3.15.5      |

> **NB:** Dependency has been **updated** if both old and new version information is present.
>
> **NB:** Dependency has been **added** if just the new version information is present.
>
> **NB:** Dependency has been **removed** if new version information isn't present.

## [[2.5](https://github.com/nf-core/viralrecon/releases/tag/2.5)] - 2022-07-13

### Enhancements & fixes

- Default Nextclade dataset shipped with the pipeline has been bumped from `2022-01-18T12:00:00Z` -> `2022-06-14T12:00:00Z`
- [[#234](https://github.com/nf-core/viralrecon/issues/234)] - Remove replacement of dashes in sample name with underscores
- [[#292](https://github.com/nf-core/viralrecon/issues/292)] - Filter empty FastQ files after adapter trimming
- [[#303](https://github.com/nf-core/viralrecon/pull/303)] - New pangolin dbs (4.0.x) not assigning lineages to Sars-CoV-2 samples in MultiQC report correctly
- [[#304](https://github.com/nf-core/viralrecon/pull/304)] - Re-factor code of `ivar_variants_to_vcf` script
- [[#306](https://github.com/nf-core/viralrecon/issues/306)] - Add contig field information in vcf header in ivar_variants_to_vcf and use bcftools sort
- [[#311](https://github.com/nf-core/viralrecon/issues/311)] - Invalid declaration val medaka_model_string
- [[#316](https://github.com/nf-core/viralrecon/issues/316)] - Variant calling isn't run when using --skip_asciigenome with metagenomic data
- [[nf-core/rnaseq#764](https://github.com/nf-core/rnaseq/issues/764)] - Test fails when using GCP due to missing tools in the basic biocontainer
- Updated pipeline template to [nf-core/tools 2.4.1](https://github.com/nf-core/tools/releases/tag/2.4.1)

### Software dependencies

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency  | Old version | New version |
| ----------- | ----------- | ----------- |
| `artic`     | 1.2.1       | 1.2.2       |
| `bcftools`  | 1.14        | 1.15.1      |
| `multiqc`   | 1.11        | 1.13a       |
| `nanoplot`  | 1.39.0      | 1.40.0      |
| `nextclade` | 1.10.2      | 2.2.0       |
| `pangolin`  | 3.1.20      | 4.1.1       |
| `picard`    | 2.26.10     | 2.27.4      |
| `quast`     | 5.0.2       | 5.2.0       |
| `samtools`  | 1.14        | 1.15.1      |
| `spades`    | 3.15.3      | 3.15.4      |
| `vcflib`    | 1.0.2       | 1.0.3       |

> **NB:** Dependency has been **updated** if both old and new version information is present.
>
> **NB:** Dependency has been **added** if just the new version information is present.
>
> **NB:** Dependency has been **removed** if new version information isn't present.

### Parameters

## [[2.4.1](https://github.com/nf-core/viralrecon/releases/tag/2.4.1)] - 2022-03-01

### Enhancements & fixes

- [[#288](https://github.com/nf-core/viralrecon/issues/288)] - `--primer_set_version` only accepts Integers (incompatible with "4.1" Artic primers set)

## [[2.4](https://github.com/nf-core/viralrecon/releases/tag/2.4)] - 2022-02-22

### Enhancements & fixes

- [nf-core/tools#1415](https://github.com/nf-core/tools/issues/1415) - Make `--outdir` a mandatory parameter
- [[#281](https://github.com/nf-core/viralrecon/issues/281)] - Nanopore medaka processing fails with error if model name, not model file, provided
- [[#286](https://github.com/nf-core/viralrecon/issues/286)] - IVAR_VARIANTS silently failing when FAI index is missing

### Parameters

| Old parameter | New parameter        |
| ------------- | -------------------- |
|               | `--publish_dir_mode` |

> **NB:** Parameter has been **updated** if both old and new parameter information is present.
>
> **NB:** Parameter has been **added** if just the new parameter information is present.
>
> **NB:** Parameter has been **removed** if new parameter information isn't present.

## [[2.3.1](https://github.com/nf-core/viralrecon/releases/tag/2.3.1)] - 2022-02-15

### Enhancements & fixes

- [[#277](https://github.com/nf-core/viralrecon/issues/277)] - Misuse of rstrip in make_variants_long_table.py script

### Software dependencies

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| `mosdepth` | 0.3.2       | 0.3.3       |
| `pangolin` | 3.1.19      | 3.1.20      |

## [[2.3](https://github.com/nf-core/viralrecon/releases/tag/2.3)] - 2022-02-04

### :warning: Major enhancements

- Please see [Major updates in v2.3](https://github.com/nf-core/viralrecon/issues/271) for a more detailed list of changes added in this version.
- When using `--protocol amplicon`, in the previous release, iVar was used for both the variant calling and consensus sequence generation. The pipeline will now perform the variant calling and consensus sequence generation with iVar and BCFTools/BEDTools, respectively.
- Bump minimum Nextflow version from `21.04.0` -> `21.10.3`

### Enhancements & fixes

- Port pipeline to the updated Nextflow DSL2 syntax adopted on nf-core/modules
- Updated pipeline template to [nf-core/tools 2.2](https://github.com/nf-core/tools/releases/tag/2.2)
- [[#209](https://github.com/nf-core/viralrecon/issues/209)] - Check that contig in primer BED and genome fasta match
- [[#218](https://github.com/nf-core/viralrecon/issues/218)] - Support for compressed FastQ files for Nanopore data
- [[#232](https://github.com/nf-core/viralrecon/issues/232)] - Remove duplicate variants called by ARTIC ONT pipeline
- [[#235](https://github.com/nf-core/viralrecon/issues/235)] - Nextclade version bump
- [[#244](https://github.com/nf-core/viralrecon/issues/244)] - Fix BCFtools consensus generation and masking
- [[#245](https://github.com/nf-core/viralrecon/issues/245)] - Mpileup file as output
- [[#246](https://github.com/nf-core/viralrecon/issues/246)] - Option to generate consensus with BCFTools / BEDTools using iVar variants
- [[#247](https://github.com/nf-core/viralrecon/issues/247)] - Add strand-bias filtering option and codon fix in consecutive positions in ivar tsv conversion to vcf
- [[#248](https://github.com/nf-core/viralrecon/issues/248)] - New variants reporting table

### Parameters

| Old parameter | New parameter                   |
| ------------- | ------------------------------- |
|               | `--nextclade_dataset`           |
|               | `--nextclade_dataset_name`      |
|               | `--nextclade_dataset_reference` |
|               | `--nextclade_dataset_tag`       |
|               | `--skip_consensus_plots`        |
|               | `--skip_variants_long_table`    |
|               | `--consensus_caller`            |
| `--callers`   | `--variant_caller`              |

> **NB:** Parameter has been **updated** if both old and new parameter information is present.
>
> **NB:** Parameter has been **added** if just the new parameter information is present.
>
> **NB:** Parameter has been **removed** if new parameter information isn't present.

### Software dependencies

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency  | Old version | New version |
| ----------- | ----------- | ----------- |
| `bcftools`  | 1.11        | 1.14        |
| `blast`     | 2.10.1      | 2.12.0      |
| `bowtie2`   | 2.4.2       | 2.4.4       |
| `cutadapt`  | 3.2         | 3.5         |
| `fastp`     | 0.20.1      | 0.23.2      |
| `kraken2`   | 2.1.1       | 2.1.2       |
| `minia`     | 3.2.4       | 3.2.6       |
| `mosdepth`  | 0.3.1       | 0.3.2       |
| `nanoplot`  | 1.36.1      | 1.39.0      |
| `nextclade` |             | 1.10.2      |
| `pangolin`  | 3.1.7       | 3.1.19      |
| `picard`    | 2.23.9      | 2.26.10     |
| `python`    | 3.8.3       | 3.9.5       |
| `samtools`  | 1.10        | 1.14        |
| `spades`    | 3.15.2      | 3.15.3      |
| `tabix`     | 0.2.6       | 1.11        |
| `vcflib`    |             | 1.0.2       |

> **NB:** Dependency has been **updated** if both old and new version information is present.
>
> **NB:** Dependency has been **added** if just the new version information is present.
>
> **NB:** Dependency has been **removed** if new version information isn't present.

## [[2.2](https://github.com/nf-core/viralrecon/releases/tag/2.2)] - 2021-07-29

### Enhancements & fixes

- Updated pipeline template to [nf-core/tools 2.1](https://github.com/nf-core/tools/releases/tag/2.1)
- Remove custom content to render Pangolin report in MultiQC as it was officially added as a module in [v1.11](https://github.com/ewels/MultiQC/pull/1458)
- [[#212](https://github.com/nf-core/viralrecon/issues/212)] - Access to `PYCOQC.out` is undefined
- [[#229](https://github.com/nf-core/viralrecon/issues/229)] - ARTIC Guppyplex settings for 1200bp ARTIC primers with Nanopore data

### Software dependencies

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency | Old version | New version |
| ---------- | ----------- | ----------- |
| `multiqc`  | 1.10.1      | 1.11        |
| `pangolin` | 3.0.5       | 3.1.7       |
| `samtools` | 1.10        | 1.12        |

> **NB:** Dependency has been **updated** if both old and new version information is present.
>
> **NB:** Dependency has been **added** if just the new version information is present.
>
> **NB:** Dependency has been **removed** if new version information isn't present.

## [[2.1](https://github.com/nf-core/viralrecon/releases/tag/2.1)] - 2021-06-15

### Enhancements & fixes

- Removed workflow to download data from public databases in favour of using [nf-core/fetchngs](https://nf-co.re/fetchngs)
- Added Pangolin results to MultiQC report
- Added warning to MultiQC report for samples that have no reads after adapter trimming
- Added docs about structure of data required for running Nanopore data
- Added docs about using other primer sets for Illumina data
- Added docs about overwriting default container definitions to use latest versions e.g. Pangolin
- Dashes and spaces in sample names will be converted to underscores to avoid issues when creating the summary metrics
- [[#196](https://github.com/nf-core/viralrecon/issues/196)] - Add mosdepth heatmap to MultiQC report
- [[#197](https://github.com/nf-core/viralrecon/issues/197)] - Output a .tsv comprising the Nextclade and Pangolin results for all samples processed
- [[#198](https://github.com/nf-core/viralrecon/issues/198)] - ASCIIGenome failing during analysis
- [[#201](https://github.com/nf-core/viralrecon/issues/201)] - Conditional include are not expected to work
- [[#204](https://github.com/nf-core/viralrecon/issues/204)] - Memory errors for SNP_EFF step

### Parameters

| Old parameter               | New parameter |
| --------------------------- | ------------- |
| `--public_data_ids`         |               |
| `--skip_sra_fastq_download` |               |

> **NB:** Parameter has been **updated** if both old and new parameter information is present.
>
> **NB:** Parameter has been **added** if just the new parameter information is present.
>
> **NB:** Parameter has been **removed** if new parameter information isn't present.

### Software dependencies

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency     | Old version | New version |
| -------------- | ----------- | ----------- |
| `nextclade_js` | 0.14.2      | 0.14.4      |
| `pangolin`     | 2.4.2       | 3.0.5       |

> **NB:** Dependency has been **updated** if both old and new version information is present.
>
> **NB:** Dependency has been **added** if just the new version information is present.
>
> **NB:** Dependency has been **removed** if new version information isn't present.

## [[2.0](https://github.com/nf-core/viralrecon/releases/tag/2.0)] - 2021-05-13

### :warning: Major enhancements

- Pipeline has been re-implemented in [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html)
- All software containers are now exclusively obtained from [Biocontainers](https://biocontainers.pro/#/registry)
- Updated minimum Nextflow version to `v21.04.0` (see [nextflow#572](https://github.com/nextflow-io/nextflow/issues/1964))
- [BCFtools](http://samtools.github.io/bcftools/bcftools.html) and [iVar](https://github.com/andersen-lab/ivar) will be run by default for Illumina metagenomics and amplicon data, respectively. However, this behaviour can be customised with the `--callers` parameter.
- Variant graph processes to call variants relative to the reference genome directly from _de novo_ assemblies have been deprecated and removed
- Variant calling with Varscan 2 has been deprecated and removed due to [licensing restrictions](https://github.com/dkoboldt/varscan/issues/12)
- New tools:
  - [Pangolin](https://github.com/cov-lineages/pangolin) for lineage analysis
  - [Nextclade](https://github.com/nextstrain/nextclade) for clade assignment, mutation calling and consensus sequence quality checks
  - [ASCIIGenome](https://asciigenome.readthedocs.io/en/latest/) for individual variant screenshots with annotation tracks

### Other enhancements & fixes

- Illumina and Nanopore runs containing the same 48 samples sequenced on both platforms have been uploaded to the nf-core AWS account for full-sized tests on release
- Initial implementation of a standardised samplesheet JSON schema to use with user interfaces and for validation
- Default human `--kraken2_db` link has been changed from Zenodo to an AWS S3 bucket for more reliable downloads
- Updated pipeline template to nf-core/tools `1.14`
- Optimise MultiQC configuration and input files for faster run-time on huge sample numbers
- [[#122](https://github.com/nf-core/viralrecon/issues/122)] - Single SPAdes command to rule them all
- [[#138](https://github.com/nf-core/viralrecon/issues/138)] - Problem masking the consensus sequence
- [[#142](https://github.com/nf-core/viralrecon/issues/142)] - Unknown method invocation `toBytes` on String type
- [[#169](https://github.com/nf-core/viralrecon/issues/169)] - ggplot2 error when generating mosdepth amplicon plot with Swift v2 primers
- [[#170](https://github.com/nf-core/viralrecon/issues/170)] - ivar trimming of Swift libraries new offset feature
- [[#175](https://github.com/nf-core/viralrecon/issues/175)] - MultiQC report does not include all the metrics
- [[#188](https://github.com/nf-core/viralrecon/pull/188)] - Add and fix EditorConfig linting in entire pipeline

### Parameters

| Old parameter                 | New parameter                         |
| ----------------------------- | ------------------------------------- |
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
|                               | `--asciigenome_read_depth`            |
|                               | `--asciigenome_window_size`           |
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
|                               | `--skip_asciigenome`                  |
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

> **NB:** Parameter has been **updated** if both old and new parameter information is present.
>
> **NB:** Parameter has been **added** if just the new parameter information is present.
>
> **NB:** Parameter has been **removed** if new parameter information isn't present.

### Software dependencies

Note, since the pipeline is now using Nextflow DSL2, each process will be run with its own [Biocontainer](https://biocontainers.pro/#/registry). This means that on occasion it is entirely possible for the pipeline to be using different versions of the same tool. However, the overall software dependency changes compared to the last release have been listed below for reference.

| Dependency                    | Old version | New version |
| ----------------------------- | ----------- | ----------- |
| `artic`                       |             | 1.2.1       |
| `asciigenome`                 |             | 1.16.0      |
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

> **NB:** Dependency has been **updated** if both old and new version information is present.
>
> **NB:** Dependency has been **added** if just the new version information is present.
>
> **NB:** Dependency has been **removed** if new version information isn't present.

## [[1.1.0](https://github.com/nf-core/viralrecon/releases/tag/1.1.0)] - 2020-06-23

### Added

- [#112](https://github.com/nf-core/viralrecon/issues/112) - Per-amplicon coverage plot
- [#124](https://github.com/nf-core/viralrecon/issues/124) - Intersect variants across callers
- [nf-core/tools#616](https://github.com/nf-core/tools/pull/616) - Updated GitHub Actions to build Docker image and push to Docker Hub
- Parameters:
  - `--min_mapped_reads` to circumvent failures for samples with low number of mapped reads
  - `--varscan2_strand_filter` to toggle the default Varscan 2 strand filter
  - `--skip_mosdepth` - skip genome-wide and amplicon coverage plot generation from mosdepth output
  - `--amplicon_left_suffix` - to provide left primer suffix used in name field of `--amplicon_bed`
  - `--amplicon_right_suffix` - to provide right primer suffix used in name field of `--amplicon_bed`
  - Unify parameter specification with COG-UK pipeline:
    - `--min_allele_freq` - minimum allele frequency threshold for calling variants
    - `--mpileup_depth` - SAMTools mpileup max per-file depth
    - `--ivar_exclude_reads` renamed to `--ivar_trim_noprimer`
    - `--ivar_trim_min_len` - minimum length of read to retain after primer trimming
    - `--ivar_trim_min_qual` - minimum quality threshold for sliding window to pass
    - `--ivar_trim_window_width` - width of sliding window
- [#118] Updated GitHub Actions AWS workflow for small and full size tests.

### Removed

- `--skip_qc` parameter

### Dependencies

- Add mosdepth `0.2.6`
- Add bioconductor-complexheatmap `2.2.0`
- Add bioconductor-biostrings `2.54.0`
- Add r-optparse `1.6.6`
- Add r-tidyr `1.1.0`
- Add r-tidyverse `1.3.0`
- Add r-ggplot2 `3.3.1`
- Add r-reshape2 `1.4.4`
- Add r-viridis `0.5.1`
- Update sra-tools `2.10.3` -> `2.10.7`
- Update bowtie2 `2.3.5.1` -> `2.4.1`
- Update picard `2.22.8` -> `2.23.0`
- Update minia `3.2.3` -> `3.2.4`
- Update plasmidid `1.5.2` -> `1.6.3`

## [[1.0.0](https://github.com/nf-core/viralrecon/releases/tag/1.0.0)] - 2020-06-01

Initial release of nf-core/viralrecon, created with the [nf-core](http://nf-co.re/) template.

This pipeline is a re-implementation of the [SARS_Cov2_consensus-nf](https://github.com/BU-ISCIII/SARS_Cov2_consensus-nf) and [SARS_Cov2_assembly-nf](https://github.com/BU-ISCIII/SARS_Cov2_assembly-nf) pipelines initially developed by [Sarai Varona](https://github.com/svarona) and [Sara Monzon](https://github.com/saramonzon) from [BU-ISCIII](https://github.com/BU-ISCIII). Porting both of these pipelines to nf-core was an international collaboration between numerous contributors and developers, led by [Harshil Patel](https://github.com/drpatelh) from the [The Bioinformatics & Biostatistics Group](https://www.crick.ac.uk/research/science-technology-platforms/bioinformatics-and-biostatistics/) at [The Francis Crick Institute](https://www.crick.ac.uk/), London. We appreciated the need to have a portable, reproducible and scalable pipeline for the analysis of COVID-19 sequencing samples and so the Avengers Assembled!

### Pipeline summary

1. Download samples via SRA, ENA or GEO ids ([`ENA FTP`](https://ena-docs.readthedocs.io/en/latest/retrieval/file-download.html), [`parallel-fastq-dump`](https://github.com/rvalieris/parallel-fastq-dump); _if required_)
2. Merge re-sequenced FastQ files ([`cat`](http://www.linfo.org/cat.html); _if required_)
3. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
4. Adapter trimming ([`fastp`](https://github.com/OpenGene/fastp))
5. Variant calling
   1. Read alignment ([`Bowtie 2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml))
   2. Sort and index alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
   3. Primer sequence removal ([`iVar`](https://github.com/andersen-lab/ivar); _amplicon data only_)
   4. Duplicate read marking ([`picard`](https://broadinstitute.github.io/picard/); _removal optional_)
   5. Alignment-level QC ([`picard`](https://broadinstitute.github.io/picard/), [`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
   6. Choice of multiple variant calling and consensus sequence generation routes ([`VarScan 2`](http://dkoboldt.github.io/varscan/), [`BCFTools`](http://samtools.github.io/bcftools/bcftools.html), [`BEDTools`](https://github.com/arq5x/bedtools2/) _||_ [`iVar variants and consensus`](https://github.com/andersen-lab/ivar) _||_ [`BCFTools`](http://samtools.github.io/bcftools/bcftools.html), [`BEDTools`](https://github.com/arq5x/bedtools2/))
      - Variant annotation ([`SnpEff`](http://snpeff.sourceforge.net/SnpEff.html), [`SnpSift`](http://snpeff.sourceforge.net/SnpSift.html))
      - Consensus assessment report ([`QUAST`](http://quast.sourceforge.net/quast))
6. _De novo_ assembly
   1. Primer trimming ([`Cutadapt`](https://cutadapt.readthedocs.io/en/stable/guide.html); _amplicon data only_)
   2. Removal of host reads ([`Kraken 2`](http://ccb.jhu.edu/software/kraken2/))
   3. Choice of multiple assembly tools ([`SPAdes`](http://cab.spbu.ru/software/spades/) _||_ [`metaSPAdes`](http://cab.spbu.ru/software/meta-spades/) _||_ [`Unicycler`](https://github.com/rrwick/Unicycler) _||_ [`minia`](https://github.com/GATB/minia))
      - Blast to reference genome ([`blastn`](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch))
      - Contiguate assembly ([`ABACAS`](https://www.sanger.ac.uk/science/tools/pagit))
      - Assembly report ([`PlasmidID`](https://github.com/BU-ISCIII/plasmidID))
      - Assembly assessment report ([`QUAST`](http://quast.sourceforge.net/quast))
      - Call variants relative to reference ([`Minimap2`](https://github.com/lh3/minimap2), [`seqwish`](https://github.com/ekg/seqwish), [`vg`](https://github.com/vgteam/vg), [`Bandage`](https://github.com/rrwick/Bandage))
      - Variant annotation ([`SnpEff`](http://snpeff.sourceforge.net/SnpEff.html), [`SnpSift`](http://snpeff.sourceforge.net/SnpSift.html))
7. Present QC and visualisation for raw read, alignment, assembly and variant calling results ([`MultiQC`](http://multiqc.info/))
