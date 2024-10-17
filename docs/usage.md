# nf-core/viralrecon: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/viralrecon/usage](https://nf-co.re/viralrecon/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Pipeline parameters

Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration except for parameters; see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

## Samplesheet format

### Illumina

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
SAMPLE_1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
SAMPLE_1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
SAMPLE_2,AEG588A2_S4_L003_R1_001.fastq.gz,
```

| Column    | Description                                                                                                                |
| --------- | -------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample.              |
| `fastq_1` | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz". |
| `fastq_2` | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz". |

> **NB:** Dashes (`-`) and spaces in sample names are automatically converted to underscores (`_`) to avoid downstream issues in the pipeline.

### Nanopore

You have the option to provide a samplesheet to the pipeline that maps sample ids to barcode ids. This allows you to associate barcode ids to clinical/public database identifiers that can be used to QC or pre-process the data with more appropriate sample names.

```console
--input '[path to samplesheet file]'
```

It has to be a comma-separated file with 2 columns. A final samplesheet file may look something like the one below:

```console
sample,barcode
21X983255,1
70H209408,2
49Y807476,3
70N209581,4
```

| Column    | Description                                                                           |
| --------- | ------------------------------------------------------------------------------------- |
| `sample`  | Custom sample name, one per barcode.                                                  |
| `barcode` | Barcode identifier attributed to that sample during multiplexing. Must be an integer. |

> **NB:** Dashes (`-`) and spaces in sample names are automatically converted to underscores (`_`) to avoid downstream issues in the pipeline.

## Nanopore input format

For Nanopore data the pipeline only supports amplicon-based analysis obtained from primer sets created and maintained by the [ARTIC Network](https://artic.network/). The [artic minion](https://artic.readthedocs.io/en/latest/commands/) tool from the [ARTIC field bioinformatics pipeline](https://github.com/artic-network/fieldbioinformatics) is used to align reads, call variants and to generate the consensus sequence.

### Nanopolish

The default variant caller used by artic minion is [Nanopolish](https://github.com/jts/nanopolish) and this requires that you provide `*.fastq`, `*.fast5` and `sequencing_summary.txt` files as input to the pipeline. These files can typically be obtained after demultiplexing and basecalling the sequencing data using [Guppy](https://nanoporetech.com/nanopore-sequencing-data-analysis) (see [ARTIC SOP docs](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html)). This pipeline requires that the files are organised in the format outlined below and gzip compressed files are also accepted:

```console
.
└── fastq_pass
    └── barcode01
        ├── FAP51364_pass_barcode01_97ca62ca_0.fastq
        ├── FAP51364_pass_barcode01_97ca62ca_1.fastq
        ├── FAP51364_pass_barcode01_97ca62ca_2.fastq
        ├── FAP51364_pass_barcode01_97ca62ca_3.fastq
        ├── FAP51364_pass_barcode01_97ca62ca_4.fastq
        ├── FAP51364_pass_barcode01_97ca62ca_5.fastq
    <TRUNCATED>
```

```console
.
└── fast5_pass
    ├── barcode01
        ├── FAP51364_pass_barcode01_97ca62ca_0.fast5
        ├── FAP51364_pass_barcode01_97ca62ca_1.fast5
        ├── FAP51364_pass_barcode01_97ca62ca_2.fast5
        ├── FAP51364_pass_barcode01_97ca62ca_3.fast5
        ├── FAP51364_pass_barcode01_97ca62ca_4.fast5
        ├── FAP51364_pass_barcode01_97ca62ca_5.fast5
    <TRUNCATED>
```

The command to run the pipeline would then be:

```console
nextflow run nf-core/viralrecon \
    --input samplesheet.csv \
    --outdir <OUTDIR> \
    --platform nanopore \
    --genome 'MN908947.3' \
    --primer_set 'artic' \
    --primer_set_version 3 \
    --fastq_dir fastq_pass/ \
    --fast5_dir fast5_pass/ \
    --sequencing_summary sequencing_summary.txt \
    -profile <docker/singularity/podman/conda/institute>
```

### Medaka

You also have the option of using [Medaka](https://github.com/nanoporetech/medaka) as an alternative variant caller to Nanopolish via the `--artic_minion_caller medaka` parameter. Medaka is faster than Nanopolish, performs mostly the same and can be run directly from `fastq` input files as opposed to requiring the `fastq`, `fast5` and `sequencing_summary.txt` files required to run Nanopolish. You must provide the appropriate [Medaka model](https://github.com/nanoporetech/medaka#models) via the `--artic_minion_medaka_model` parameter if using `--artic_minion_caller medaka`. The `fastq` files have to be organised in the same way as for Nanopolish as outlined in the section above.

The command to run the pipeline would then be:

```console
nextflow run nf-core/viralrecon \
    --input samplesheet.csv \
    --outdir <OUTDIR> \
    --platform nanopore \
    --genome 'MN908947.3' \
    --primer_set 'artic' \
    --primer_set_version 3 \
    --fastq_dir fastq_pass/ \
    --artic_minion_caller medaka \
    --artic_minion_medaka_model r941_min_high_g360 \
    -profile <docker/singularity/podman/conda/institute>
```

## Illumina primer sets

The Illumina processing mode of the pipeline has been tested on numerous different primer sets. Where possible we are trying to collate links and settings for standard primer sets to make it easier to run the pipeline with standard parameter keys. If you are able to get permissions from the vendor/supplier to share the primer information then we would be more than happy to support it within the pipeline.

For SARS-CoV-2 data we recommend using the "MN908947.3" genome because it is supported out-of-the-box by the most commonly used primer sets available from the [ARTIC Network](https://artic.network/). For ease of use, we are also maintaining a version of the "MN908947.3" genome along with the appropriate links to the ARTIC primer sets in the [genomes config file](https://github.com/nf-core/configs/blob/master/conf/pipeline/viralrecon/genomes.config) used by the pipeline. The genomes config file can be updated independently from the main pipeline code to make it possible to dynamically extend this file for other viral genomes/primer sets on request.

For further information or help, don't hesitate to get in touch on the [Slack `#viralrecon` channel](https://nfcore.slack.com/channels/viralrecon) (you can join with [this invite](https://nf-co.re/join/slack)).

### ARTIC primer sets

An example command using v3 ARTIC primers with "MN908947.3":

```console
nextflow run nf-core/viralrecon \
    --input samplesheet.csv \
    --outdir <OUTDIR> \
    --platform illumina \
    --protocol amplicon \
    --genome 'MN908947.3' \
    --primer_set artic \
    --primer_set_version 3 \
    --skip_assembly \
    -profile <docker/singularity/podman/conda/institute>
```

### SWIFT primer sets

The [SWIFT amplicon panel](https://swiftbiosci.com/swift-amplicon-sars-cov-2-panel/) is another commonly used method used to prep and sequence SARS-CoV-2 samples. We haven't been able to obtain explicit permission to host standard SWIFT primer sets but you can obtain a masterfile which is freely available from their website that contains the primer sequences as well as genomic co-ordinates. You just need to convert this file to [BED6](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format and provide it to the pipeline with `--primer_bed swift_primers.bed`. Be sure to check the values provided to `--primer_left_suffix` and `--primer_right_suffix` match the primer names defined in the BED file as highlighted in [this issue](https://github.com/nf-core/viralrecon/issues/169). For an explanation behind the usage of the `--ivar_trim_offset 5` for SWIFT primer sets see [this issue](https://github.com/nf-core/viralrecon/issues/170).

An example command using SWIFT primers with "MN908947.3":

```console
nextflow run nf-core/viralrecon \
    --input samplesheet.csv \
    --outdir <OUTDIR> \
    --platform illumina \
    --protocol amplicon \
    --genome 'MN908947.3' \
    --primer_bed swift_primers.bed \
    --primer_left_suffix '_F' \
    --primer_right_suffix '_R' \
    --ivar_trim_offset 5 \
    --skip_assembly \
    -profile <docker/singularity/podman/conda/institute>
```

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/viralrecon --input samplesheet.csv --outdir <OUTDIR> --genome 'MN908947.3' -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/viralrecon -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/viralrecon
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/viralrecon releases page](https://github.com/nf-core/viralrecon/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

## Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.
:::

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

#### Freyja

[Freyja](https://github.com/andersen-lab/Freyja) relies on a dataset of barcodes that use lineage defining mutations (see [UShER](https://usher-wiki.readthedocs.io/en/latest/#)). By default the most recent barcodes will be downloaded and used. However, if analyses need to be compared across multiple datasets, it might be of interest to re-use the same barcodes, or to rerun all Freyja analyses with the most recent dataset. To do this, specify the barcodes and lineages using the `--freyja_barcodes`, `--freyja_lineages` parameters, respectivly. The boostrapping of Freyja can be skipped by specifying `--skip_freyja_boot`.

### Cutadapt

According to [Cutadapt's documentation regarding adapter types](https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types), you can have:

- Regular 3’ adapter: `-a ADAPTER`
  - Set `--skip_noninternal_primers` to `true`
  - Change `modules_illumina.config` > `CUTADAPT` > `ext.args` to use `-a` instead of `-g`
- Regular 5’ adapter: `-g ADAPTER`
  - Set `--skip_noninternal_primers` to `true`
- Non-internal 3’ adapter: `-a ADAPTERX`:
  - Change `modules_illumina.config` > `PREPARE_PRIMER_FASTA` > `ext.args` to use `$` instead of `^` to add the X at the end of the sequence.
  - Change `modules_illumina.config` > `CUTADAPT` > `ext.args` to use `-a` instead of `-g`
- Non-internal 5’ adapter: `-g XADAPTER`: This is the option by default.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
