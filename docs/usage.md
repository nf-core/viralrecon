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

Artic minion requires that you provide `*.fastq` files as input to the pipeline. These files can typically be obtained after demultiplexing and basecalling the sequencing data using [Guppy](https://nanoporetech.com/nanopore-sequencing-data-analysis) (see [ARTIC SOP docs](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html)). This pipeline requires that the files are organised in the format outlined below and gzip compressed files are also accepted:

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
    --sequencing_summary sequencing_summary.txt \
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

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/viralrecon -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/viralrecon
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/viralrecon releases page](https://github.com/nf-core/viralrecon/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

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
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
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

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

#### Freyja

[Freyja](https://github.com/andersen-lab/Freyja) relies on a dataset of barcodes that use lineage defining mutations (see [UShER](https://usher-wiki.readthedocs.io/en/latest/#)). By default the most recent barcodes will be downloaded and used. However, if analyses need to be compared across multiple datasets, it might be of interest to re-use the same barcodes, or to rerun all Freyja analyses with the most recent dataset. To do this, specify the barcodes and lineages using the `--freyja_barcodes`, `--freyja_lineages` parameters, respectivly. The boostrapping of Freyja can be skipped by specifying `--skip_freyja_boot`.

### Cutadapt

According to [Cutadapt's documentation regarding adapter types](https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types), you can have:

- Regular 3’ adapter: `-a ADAPTER`
  - Set `--skip_noninternal_primers` to `true`
  - Set `--threeprime_adapters` to `true`
- Regular 5’ adapter: `-g ADAPTER`
  - Set `--skip_noninternal_primers` to `true`
- Non-internal 3’ adapter: `-a ADAPTERX`:
  - Change `modules_illumina.config` > `PREPARE_PRIMER_FASTA` > `ext.args` to use `$` instead of `^` to add the X at the end of the sequence.
  - Set `--threeprime_adapters` to `true`
- **Non-internal 5’ adapter**: `-g XADAPTER`: **This is the option by default**.

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

## Enterovirus

`nf-core/viralrecon` has been extended to support **Enterovirus typing** through optimized BLASTN searches with taxonomic ID filtering.

The goal of this addition is to enable rapid and accurate **Enterovirus strain/genotype identification** directly from consensus sequences in the de novo assembly track, reducing BLASTN runtime while maintaining high typing accuracy.

### Enterovirus profile and recommended params

The Enterovirus profile introduces specific adjustments compared to a standard `viralrecon` run, primarily to optimize **strain/genotype typing**.

By default, variant calling is deactivated for enterovirus typing because de novo assembly performs better through assembled contigs and BLAST comparison.

#### De novo assembly

De novo assembly is performed using **spades** and/or **unicycler** to generate contigs for BLASTN typing.
Note that **minia** is not a default assembler for enterovirus.

#### BLASTN

Typing through **blastn** is performed through supplying a blast database with taxid mapping and a taxidlist, to improve typing specificity and reduce runtime.
The new optional parameter, `--taxidlist`, allows users to provide a list of **NCBI Taxonomy IDs** corresponding to Enterovirus taxa of interest.

Taxonomy IDs related to enterovirus can be retriewed through // Add command here

Example usage:

```bash
--blastdb path/to/blastdb
--taxidlist path/to/taxidlist_taxids.txt
--profile test_ev
```

#### Reference genome

The default reference used in the enterovirus config of `nf-core/viralrecon` is the REFSEQ genome [NC_002058.3](https://www.ncbi.nlm.nih.gov/nuccore/NC_002058.3/) from NCBI of Enterovirus C. Given that a blast database is supplied, this reference is not significant for the de novo assembly. However, if no blast database is supplied, the reference genome will be used as reference for blast searches.

Users may provide their **own custom reference genomes** using the parameters:

```bash
--fasta <path_to_reference.fasta>
--gff   <path_to_annotation.gff>
```

## HIV

`nf-core/viralrecon` has been optimized to support **HIV genome reconstruction and resistance detection**.

This implementation takes inspiration from the parameterization used in the [**HIVdb** pipeline at Stanford University](https://hivdb.stanford.edu/hivdb/by-reads/).

The goal is to adapt the viralrecon framework to ensure accurate minor variant calling, reliable consensus sequence generation, and compatibility with resistance analysis tools for HIV.

### HIV profile and recomended params

The HIV profile (`-profile hiv`) introduces specific adjustments compared to a standard `viralrecon` run, mainly in the variant calling and consensus generation steps, to align with the requirements of HIV resistance interpretation pipelines.

⚠️ By default, de novo assembly is deactivated for HIV as the resistance detection is only performed for the consensus genome through mapping approach, which has been already benchmarked agains Stanford's [HIVdb](https://hivdb.stanford.edu/hivdb/by-reads/) results.

#### **Reference genome**

The default reference used in the HIV config of `nf-core/viralrecon` is the one **derived from the codfreq JSON profile**. This reference has been specifically generated from the [codfreq .json profile](https://github.com/hivdb/codfreq/blob/main/profiles/HIV1.json). By aligning to this codfreq-based reference, the resulting codon frequencies and codon coverages are directly comparable to those produced by [**HIVdb**](https://hivdb.stanford.edu/hivdb/by-reads/), ensuring accurate interpretation of resistance data.

However, users may provide their **own custom reference genome** using the parameters:

```bash
--fasta <path_to_reference.fasta>
--gff   <path_to_annotation.gff>
```

When a custom reference is supplied, all variant calling, consensus generation, and resistance detection steps will use this user-defined reference.

⚠️ Note that this may produce results that are not directly comparable to those obtained with HIVdb. Users employing custom references do so at their own discretion and risk, as the chosen reference can affect the INDELs calling and their coverage and frequency.

> **Important**: When performing HIV resistance detection, an annotation file (.gff) is mandatory. The GFF file is required to correctly annotate the pol gene regions (PR, RT and IN), which are essential for both sierra-local and codfreq analyses.

#### **Available genome options**

Two HIV genome references are distributed with viralrecon and can be selected using the `--genome` parameter:

- `--genome codfreq`: The codfreq-derived reference (DEFAULT), generated from the original codfreq JSON annotation files. This reference ensures that codon-level frequency tables and coverage results are fully compatible with HIVdb and sierra-local resistance prediction outputs.
- `--genome 'NC_001802.1'`: nextclade refeence, based on NC_001802.1 (HXB2 genome, K03455), the reference used by Nextclade for HIV-1 analysis.

#### **HIV-specific parameters**

The following parameters have been added to handle **HIV resistance detection** and related tools.

- General HIV resistance options:
  - **`perform_hiv_resistance`**: Whether to perform HIV resistance analysis or not.
- Options related with files required by `sierra-local` software for resistance detecton in HIV. If not provided, the files included with the software will be used.
  - **`hivdb_xml`**: Path to the HIVdb ASI2 XML file. Updated files can be found [here](https://github.com/hivdb/hivfacts/tree/main/data/algorithms), where the different algorithm versions are stored. The one included in the software is the one with `HIVDB_9.8.xml` name.
  - **`apobec_drm`**: Path to the JSON HIVdb APOBEC DRM file. Updated files can be found [here](https://github.com/hivdb/hivfacts/tree/main/data/apobecs) with the `apobec_drms.json` name.
  - **`apobec_csv`**: Path to the CSV APOBEC file. Updated files can be found [here](https://github.com/hivdb/hivfacts/tree/main/data/apobecs) with the `apobecs.csv` name.
  - **`unusual_csv`**: Path to the CSV file used to determine unusual mutations. The file included with the software can be found [here](https://github.com/PoonLab/sierra-local/tree/master/sierralocal/data) with the name `rx-all_subtype-all.csv`
  - **`sdrms_csv`**: Path to the CSV file used to define SDRM mutations. Updated files can be found [here](https://github.com/hivdb/hivfacts/tree/main/data/) with the `sdrms_hiv1.csv` name.
  - **`mutation_csv`**: Path to the CSV file defining mutation types. Updated files can be found [here](https://github.com/hivdb/hivfacts/tree/main/data/) with the `mutation-type-pairs_hiv1.csv` name.

#### **Nextclade configuration**

The following default parameters are used in the `--genome codfreq` and `--genome 'NC_001802.1'` for **Nextclade HIV-1 analysis**. Make sure to update `nextclade_dataset_tag` to the latest one so the pipeline can download it.

```bash
nextclade_dataset_tag  = '2025-09-09--12-13-13Z'
nextclade_dataset_name = 'neherlab/hiv-1'
nextclade_dataset      = false
```

#### **HIV-specific configuration**

When using `--genome 'codfreq'` or `--genome 'NC_001802.1'` or `--perform_hiv_resistance`, which are activated in the `-profile hiv` configuration, the following configuration is activated by default.

_Variant calling_

Variant calling is performed using **iVar variants** with parameters optimized for HIV intra-host diversity:

- **Variant frequency threshold:** `0.01`: Variants are called if their allele frequency is at least 1%, allowing detection of minority variants relevant for resistance analysis.
- **Minimum base quality:** `30`: Ensures reliability of detected variants.
- **Minimum depth of coverage:** `50`: Matches the minimum coverage used in HIVdb protocols to support robust variant calling.

_Consensus generation_

Consensus sequences are generated using **iVar consensus**, which supports ambiguous nucleotide codes, a key feature for HIV resistance interpretation tools.  
The parameters are defined as follows:

- **Frequency threshold:** `0.9`: A variant is included in the consensus until a position's allele frequency is representing at least 90% of the allele frequency. This means that variants with frequencies below 50% may still be incorporated if necessary to achieve the 90% threshold in a given position.
- **Minimum base quality:** `30`
- **Minimum depth of coverage:** `50`
- **Low coverage regions:**: Positions with fewer than 10 reads are masked as `N`.

The user can change this configuration by providing a custom config file with `-c` param. This config file should contain this piece of code with the specific configuration for `IVAR_VARIANTS` and `IVAR_CONSENSUS`. `RESISTANCE_REPORT`'s `ext.args2` should be the same as `IVAR_CONSENSUS`'s `ext.args` for the final report to contain the appropriate information.

```bash
process {
    withName: 'IVAR_CONSENSUS' {
        ext.args = '-t 0.9 -q 30 -m 50 -n N'
    }
    withName: 'IVAR_VARIANTS' {
        ext.args = '-t 0.01 -q 30 -m 50'
    }
    withName: 'RESISTANCE_REPORT' {
        ext.args2 = '-t 0.9 -q 30 -m 50 -n N'
    }
}
```

### HIV resistance detection protocol

The detection of HIV drug resistance in `nf-core/viralrecon` is performed using the [**sierra-local**](https://github.com/PoonLab/sierra-local/) software, complemented with codon-level frequency analysis generated by a modified version of [**codfreq**](https://github.com/hivdb/codfreq/).

This integrated approach ensures that both **mutation interpretation** and **codon frequency information** are available for downstream resistance reporting.

#### **1. Sierra-local: HIV resistance prediction**

[**sierra-local**](https://github.com/PoonLab/sierra-local/) is a Python3 implementation of the [Stanford University HIV Drug Resistance Database](https://hivdb.stanford.edu/) (HIVdb) [Sierra web service](https://hivdb.stanford.edu/page/webservice/).

It allows laboratories to generate **HIV-1 drug resistance predictions locally**, without the need to transmit patient data over a network. This provides full control over **data provenance**, **security**, and **regulatory compliance**.

> _Sierra-local computes drug resistance predictions directly from consensus HIV-1 sequences, producing JSON-formatted output that includes key resistance-associated mutations and drug susceptibility scores._

Within the HIV module of `nf-core/viralrecon`, sierra-local is executed using the consensus .FASTA file generated for each sample. The output is a `.json` file containing the main resistance data.

However, sierra-local operates at the **codon level**, and its JSON output does not include **codon-level frequencies** or **codon depth information**.  
These additional metrics are crucial for resistance reporting and neeed consistency with other viralrecon outputs, which are based on nucleotide-level variant calls.

#### **2. Codon-level frequency analysis with codfreq**

To obtain codon frequency information, `nf-core/viralrecon` integrates a modified version of [**codfreq**](https://github.com/hivdb/codfreq/) tool.

The modified version of codfreq produces a codon frequency table in the _CodFreq_ format, which contains seven columns:

| Column                | Description                                                         |
| :-------------------- | :------------------------------------------------------------------ |
| `gene`                | HIV gene (PR, RT, or IN)                                            |
| `position`            | Codon position at protein level                                     |
| `total`               | Total reads covering that position (tiplet)                         |
| `codon`               | Nucleotide triplet or INDELS                                        |
| `count`               | Total reads supporting that codon                                   |
| `total_quality_score` | Qualituy score of that codon assigned by codfreq                    |
| `aa_codon`            | If the codon is multiple of 3, the corresponding aminoacid sequence |

This format enables codon-level analysis of intra-sample diversity.

However, codfreq requires a **custom JSON file** that describes the HIV gene structure (PR, RT, and IN). To generate this JSON dynamically according to the user’s selected reference genome, several preprocessing steps are performed.

#### **3. Generation of codfreq annotation JSON**

The generation of the required JSON file for codfreq involves the following steps:

1. **Conversion of reference annotations**: The original codfreq HIV profile (in JSON format) has been converted into a `.fasta` and `.gff` file containing the necessary annotations for PR, RT, and IN genes. These files are included within `nf-core/viralrecon` assets.

2. **Reference genome alignment with Liftoff**: Only performed when the provided genome is not the one from codfreq. Using the provided reference genome and the `.fasta` and `.gff` from the step 1, annotation files are processed with **[Liftoff](https://github.com/agshumate/Liftoff)** to transfer the HIV gene annotations to the user’s specific reference genome.

3. **Custom JSON generation**: The annotated `.gff` file and the reference `.fasta` are used as input for a **custom Python script** `bin/gff2json.py` within `nf-core/viralrecon`. This script produces the JSON file required by `codfreq`, ensuring accurate gene coordinates and consistency with the reference genome.

4. **Annotation harmonization for variant mapping**: Only performed when the provided genome is not the one from codfreq. Once the new `.gff` file is produced by Liftoff, an additional annotation step is performed. The variant caller outputs (based on **nucleotide-level variants**) are re-annotated using this updated `.gff`, generating a **secondary annotation layer**. This process allows a direct mapping between **codon-based resistance results** (from `sierra-local` and `codfreq`) and **nucleotide-based variant calls**, ensuring that both analyses can be compared and integrated accurately in downstream reviews.

This automated workflow guarantees that the codfreq JSON file always matches the user-defined HIV reference, maintaining compatibility with both sierra-local and the consensus sequence data.

#### **4. Integration of codfreq with mapped reads**

Once the annotation JSON is generated, the **mapped BAM files** and the **codfreq JSON** are processed using a **modified version of codfreq** within the pipeline.  
This version has been adapted to work seamlessly with `nf-core/viralrecon` output formats and file structures.

The result is a set of `.codfreq` tables containing codon-level read frequencies for PR, RT, and IN genes, harmonized with the annotation used by Sierra-local.

#### **5. Integration and final reporting**

Finally, the outputs from **sierra-local** (`.json`) and **codfreq** (`.codfreq`) are merged to produce custom HIV resistance summary tables.

The resulting integrated reports provide a **comprehensive view of HIV drug resistance**, maintaining compatibility with established databases and ensuring consistency with the rest of the `viralrecon` analytical framework.

To know more about the output files generated in the HIV resistance steps, check the [output documentation](output.md#illumina-hiv-resistance-detection).
