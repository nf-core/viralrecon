# nf-core/viralrecon: Usage

## Table of contents

* [Table of contents](#table-of-contents)
* [Introduction](#introduction)
* [Running the pipeline](#running-the-pipeline)
  * [Updating the pipeline](#updating-the-pipeline)
  * [Reproducibility](#reproducibility)
* [Main arguments](#main-arguments)
  * [`-profile`](#-profile)
  * [`--input`](#--input)
  * [`--protocol`](#--protocol)
  * [`--amplicon_bed`](#--amplicon_bed)
  * [`--amplicon_fasta`](#--amplicon_fasta)
* [SRA download](#sra-download)
  * [`--ignore_sra_errors`](#--ignore_sra_errors)
  * [`--save_sra_fastq`](#--save_sra_fastq)
  * [`--skip_sra`](#--skip_sra)
* [Reference genomes](#reference-genomes)
  * [`--genome`](#--genome)
  * [`--fasta`](#--fasta)
  * [`--gff`](#--gff)
  * [`--save_reference`](#--save_reference)
* [Kraken](#kraken)
  * [`--kraken2_db`](#--kraken2_db)
  * [`--kraken2_db_name`](#--kraken2_db_name)
  * [`--kraken2_use_ftp`](#--kraken2_use_ftp)
  * [`--save_kraken2_fastq`](#--save_kraken2_fastq)
* [Read trimming](#read-trimming)
  * [`--skip_adapter_trimming`](#--skip_adapter_trimming)
  * [`--skip_amplicon_trimming`](#--skip_amplicon_trimming)
  * [`--save_trimmed`](#--save_trimmed)
* [Variant calling](#variant-calling)
  * [`--callers`](#-callers)
  * [`--ivar_exclude_reads`](#--ivar_exclude_reads)
  * [`--save_align_intermeds`](#--save_align_intermeds)
  * [`--save_pileup`](#--save_pileup)
  * [`--skip_snpeff`](#--skip_snpeff)
  * [`--skip_variants_quast`](#--skip_variants_quast)
  * [`--skip_variants`](#--skip_variants)
* [De novo assembly](#de-novo-assembly)
  * [`--assemblers`](#--assemblers)
  * [`--skip_blast`](#--skip_blast)
  * [`--skip_abacas`](#--skip_abacas)
  * [`--skip_plasmidid`](#--skip_plasmidid)
  * [`--skip_assembly_quast`](#--skip_assembly_quast)
  * [`--skip_assembly`](#--skip_assembly)  
* [Skipping QC steps](#skipping-qc-steps)
  * `--skip_fastqc`
  * `--skip_picard_metrics`
  * `--skip_multiqc`
  * `--skip_qc`
* [Job resources](#job-resources)
  * [Automatic resubmission](#automatic-resubmission)
  * [Custom resource requests](#custom-resource-requests)
* [AWS Batch specific parameters](#aws-batch-specific-parameters)
  * [`--awsqueue`](#--awsqueue)
  * [`--awsregion`](#--awsregion)
  * [`--awscli`](#--awscli)
* [Other command line parameters](#other-command-line-parameters)
  * [`--outdir`](#--outdir)
  * [`--email`](#--email)
  * [`--email_on_fail`](#--email_on_fail)
  * [`--max_multiqc_email_size`](#--max_multiqc_email_size)
  * [`-name`](#-name)
  * [`-resume`](#-resume)
  * [`-c`](#-c)
  * [`--custom_config_version`](#--custom_config_version)
  * [`--custom_config_base`](#--custom_config_base)
  * [`--max_memory`](#--max_memory)
  * [`--max_time`](#--max_time)
  * [`--max_cpus`](#--max_cpus)
  * [`--plaintext_email`](#--plaintext_email)
  * [`--monochrome_logs`](#--monochrome_logs)
  * [`--multiqc_config`](#--multiqc_config)

## Introduction

Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler. Finally, you can use nextflow `-bg` flag to execute nextflow in background.

It is recommended to limit the Nextflow Java virtual machines memory. We recommend adding the following line to your environment (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/viralrecon --input samplesheet.csv --genome 'NC_045512.2' -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/viralrecon
```

### Reproducibility

It's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/viralrecon releases page](https://github.com/nf-core/viralrecon/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

## Main arguments

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Conda) - see below.

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

* `docker`
  * A generic configuration profile to be used with [Docker](http://docker.com/)
  * Pulls software from dockerhub: [`nfcore/viralrecon`](http://hub.docker.com/r/nfcore/viralrecon/)
* `singularity`
  * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
  * Pulls software from DockerHub: [`nfcore/viralrecon`](http://hub.docker.com/r/nfcore/viralrecon/)
* `conda`
  * Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker or Singularity.
  * A generic configuration profile to be used with [Conda](https://conda.io/docs/)
  * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `test`
  * A profile with a complete configuration for automated testing
  * Includes links to test data so needs no other parameters

### `--input`

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

#### Format

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once (e.g. to increase sequencing depth). The pipeline will perform the analysis in parallel, and subsequently merge them when required.

A final design file may look something like the one below. `SRR10903401` was sequenced twice in Illumina PE format, `SRR11241255` was sequenced once in Illumina SE format, and `SRR11092056` and `SRR11177792` need to be downloaded from the SRA before the main pipeline execution.

```bash
sample,fastq_1,fastq_2
SRR10903401,SRR10903401_1.fastq.gz,SRR10903401_2.fastq.gz
SRR10903401,SRR10903402_1.fastq.gz,SRR10903402_2.fastq.gz
SRR11241255,SRR11241255.fastq.gz,
SRR11092056,,
SRR11177792,,
```

| Column    | Description                                                                                                                            |
|-----------|----------------------------------------------------------------------------------------------------------------------------------------|
| `sample`  | Custom sample name or SRA 'SR' or 'PR' identifier. This will be identical for multiple sequencing libraries/runs from the same sample. |
| `fastq_1` | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".             |
| `fastq_2` | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".             |

### `--protocol`

Specifies the type of protocol used for sequencing i.e. 'metagenomic' or 'amplicon' (Default: 'metagenomic').

### `--amplicon_bed`

Viral genome location of primers. Mandatory when `--protocol amplicon` and not `--skip_mapping`.

#### Format

It must be in bed format, with this fields:

```Bash

chr\tstart_primer\tend_primer\tname\tqual\tstrand

```

Example:

```Bash

NC_045512.2 30 54 nCoV-2019_1_LEFT 60 -
NC_045512.2 385 410 nCoV-2019_1_RIGHT 60 +
NC_045512.2 320 342 nCoV-2019_2_LEFT 60 -
NC_045512.2 704 726 nCoV-2019_2_RIGHT 60 +

```

### `--amplicon_fasta`

Primer sequences in fasta format. Mandatory when `--protocol amplicon` and not `--skip_assembly`.
Example:

```Bash

>nCoV-2019_1_LEFT
ACCAACCAACTTTCGATCTCTTGT
>nCoV-2019_1_RIGHT
CATCTTTAAGATGTTGACGTGCCTC
>nCoV-2019_2_LEFT
CTGTTTTACAGGTTCGCGACGT
>nCoV-2019_2_RIGHT
TAAGGATCAGTGCCAAGCTCGT
>nCoV-2019_3_LEFT
CGGTAATAAAGGAGCTGGTGGC
>nCoV-2019_3_RIGHT
AAGGTGTCTGCAATTCATAGCTCT

```

## SRA download

## `--ignore_sra_errors`

Ignore validation errors when checking SRA identifiers that would otherwise cause the pipeline to fail (Default: false).

## `--save_sra_fastq`

Save FastQ files created from SRA identifiers in the results directory (Default: false).

## `--skip_sra`

Skip steps involving the download and validation of FastQ files using SRA identifiers (Default: false).

## Reference genomes

### `--genome`

This parameter allows you to provide a key for the viral genome you would like to use with the pipeline. To run the pipeline, you must specify which to use with the `--genome` flag.

Note that you can use the same configuration setup to save sets of reference files for your own use. See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for instructions on where to save such a file.

The syntax for this reference configuration is as follows:

```nextflow
params {
  // Genome reference file paths
  genomes {
    'NC_045512.2' {
      fasta = "<path to the genome fasta file>"
      gff   = "<path to the genome annotation file>"
    }
    'MN908947.3' {
      fasta = "<path to the genome fasta file>"
      gff   = "<path to the genome annotation file>"
    }
    // Any number of additional genomes, key is used with --genome
  }
}
```

You can find the keys to specify the genomes in the [Genomes config file](https://github.com/nf-core/configs/blob/master/conf/pipeline/viralrecon/genomes.config).

### `--fasta`

Full path to fasta file containing reference genome for the viral species (*mandatory* if `--genome` is not specified). If you don't have a Bowtie2 index available this will be generated for you automatically. Combine with `--save_reference` to save Bowtie2 index for future runs.

```bash
--fasta '[path to FASTA reference]'
```

### `--gff`

Full path to viral gff annotation file (Default: false).

### `--save_reference`

If the Bowtie2 index is generated by the pipeline use this parameter to save it to your results folder. These can then be used for future pipeline runs, reducing processing times (Default: false).

## Kraken2

### `--kraken2_db`

Full path to Kraken2 database built from host genome (Default: 'https://zenodo.org/record/3738199/files/kraken2_human.tar.gz').

### `--kraken2_db_name`

Name for host genome as recognised by Kraken2 when using the `kraken2 build` command (Default: 'human').

### `--kraken2_use_ftp`

Option for kraken using ftp download instead of rsync (Default: false).

### `--save_kraken2_fastq`

Save the host and viral fastq files in the results directory (Default: false).

## Read trimming

### `--skip_adapter_trimming`

Skip the adapter trimming step performed by fastp. Use this if your input FastQ files have already been trimmed outside of the workflow or if you're very confident that there is no adapter contamination in your data (Default: false).

### `--skip_amplicon_trimming`

Skip the amplicon trimming step performed by Cutadapt. Use this if your input FastQ files have already been trimmed outside of the workflow or if you're very confident that there is no primer sequence contamination in your data (Default: false).

### `--save_trimmed`

By default, trimmed FastQ files will not be saved to the results directory. Specify this flag (or set to true in your config file) to copy these files to the results directory when complete (Default: false).

## Alignments

### `--ivar_exclude_reads`

This option unsets the `-e` parameter in `ivar trim` to discard reads without primers (Default: false).

### `--save_align_intermeds`

By default, intermediate BAM files will not be saved. The final BAM files created after the appropriate filtering step are always saved to limit storage usage. Set to true to also save other intermediate BAM files (Default: false).

## Variant calling

### `--callers`

Specify which variant calling algorithms you would like to use. Available options are `varscan2` and `ivar` (Default: 'varscan2,ivar').

### `--save_pileup`

Save Pileup files in the results directory. These tend to be quite large so are not saved by default (Default: false).

### `--skip_snpeff`

Skip SnpEff and SnpSift annotation of variants (Default: false).

### `--skip_variants_quast`

Skip generation of QUAST aggregated report for consensus sequences (Default: false).

### `--skip_variants`

Specify this parameter to skip all of the variant calling and mapping steps in the pipeline (Default: false).

## De novo assembly

### `--assemblers`

Specify which assembly algorithms you would like to use. Available options are `spades`, `metaspades` and `unicycler` (Default: 'spades,metaspades,unicycler').

### `--skip_blast`

Skip blastn of assemblies relative to reference genome (Default: false)

### `--skip_abacas`

Skip ABACUS process for assembly contiguation (Default: false).

### `--skip_plasmidid`

Skip assembly report generation by PlasmidID (Default: false).

### `--skip_assembly_quast`

Skip generation of QUAST aggregated report for assemblies (Default: false).

### `--skip_assembly`

Specify this parameter to skip all of the de novo assembly steps in the pipeline (Default: false).

## Skipping QC steps

The pipeline contains a large number of quality control steps. Sometimes, it may not be desirable to run all of them if time and compute resources are limited.
The following options make this easy:

| Step                      | Description                                              |
|---------------------------|----------------------------------------------------------|
| `--skip_fastqc`           | Skip FastQC                                              |
| `--skip_picard_metrics`   | Skip Picard CollectMultipleMetrics and CollectWgsMetrics |
| `--skip_multiqc`          | Skip MultiQC                                             |
| `--skip_qc`               | Skip all QC steps except for MultiQC                     |

## Job resources

### Automatic resubmission

Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

### Custom resource requests

Wherever process-specific requirements are set in the pipeline, the default value can be changed by creating a custom config file. See the files hosted at [`nf-core/configs`](https://github.com/nf-core/configs/tree/master/conf) for examples.

If you are likely to be running `nf-core` pipelines regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter (see definition below). You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack).

## AWS Batch specific parameters

Running the pipeline on AWS Batch requires a couple of specific parameters to be set according to your AWS Batch configuration. Please use [`-profile awsbatch`](https://github.com/nf-core/configs/blob/master/conf/awsbatch.config) and then specify all of the following parameters.

### `--awsqueue`

The JobQueue that you intend to use on AWS Batch.

### `--awsregion`

The AWS region in which to run your job. Default is set to `eu-west-1` but can be adjusted to your needs.

### `--awscli`

The [AWS CLI](https://www.nextflow.io/docs/latest/awscloud.html#aws-cli-installation) path in your custom AMI (Default: `/home/ec2-user/miniconda/bin/aws`).

Please make sure to also set the `-w/--work-dir` and `--outdir` parameters to a S3 storage bucket of your choice - you'll get an error message notifying you if you didn't.

## Other command line parameters

### `--outdir`

The output directory where the results will be saved.

### `--email`

Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

### `--email_on_fail`

This works exactly as with `--email`, except emails are only sent if the workflow is not successful.

### `--max_multiqc_email_size`

Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB).

### `-name`

Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`

Specify the path to a specific config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### `--custom_config_version`

Provide git commit id for custom Institutional configs hosted at `nf-core/configs`. This was implemented for reproducibility purposes (Default: `master`).

```bash
## Download and use config file with following git commid id
--custom_config_version d52db660777c4bf36546ddb188ec530c3ada1b96
```

### `--custom_config_base`

If you're running offline, nextflow will not be able to fetch the institutional config files
from the internet. If you don't need them, then this is not a problem. If you do need them,
you should download the files from the repo and tell nextflow where to find them with the
`custom_config_base` option. For example:

```bash
## Download and unzip the config files
cd /path/to/my/configs
wget https://github.com/nf-core/configs/archive/master.zip
unzip master.zip

## Run the pipeline
cd /path/to/my/data
nextflow run /path/to/pipeline/ --custom_config_base /path/to/my/configs/configs-master/
```

> Note that the nf-core/tools helper package has a `download` command to download all required pipeline
> files + singularity containers + institutional configs in one go for you, to make this process easier.

### `--max_memory`

Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`

### `--max_time`

Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`

Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

### `--plaintext_email`

Set to receive plain-text e-mails instead of HTML formatted.

### `--monochrome_logs`

Set to disable colourful command line output and live life in monochrome.

### `--multiqc_config`

Specify a path to a custom MultiQC configuration file.
