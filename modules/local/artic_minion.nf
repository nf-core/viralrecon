// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process ARTIC_MINION {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::artic=1.2.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/artic:1.2.1--py_0"
    } else {
        container "quay.io/biocontainers/artic:1.2.1--py_0"
    }

    input:
    tuple val(meta), path(fastq)
    path  fast5_dir
    path  sequencing_summary
    path  scheme_dir
    val   scheme
    val   scheme_version

    output:
    tuple val(meta), path("${prefix}.*")                          , emit: results
    tuple val(meta), path("${prefix}.primertrimmed.rg.sorted.bam"), emit: bam_primer_trim
    tuple val(meta), path("${prefix}.sorted.bam")                 , emit: bam
    tuple val(meta), path("${prefix}.consensus.fasta")            , emit: fasta
    tuple val(meta), path("${prefix}.pass.vcf.gz")                , emit: vcf
    path  "*.version.txt"                                         , emit: version

    script:
    def software = getSoftwareName(task.process)
    prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def version  = scheme_version.toLowerCase().replaceAll('v','')
    def fast5    = params.fast5_dir          ? "--fast5-directory $fast5_dir"             : ""
    def summary  = params.sequencing_summary ? "--sequencing-summary $sequencing_summary" : ""
    if (options.args.contains('--medaka ')) {
        fast5    = ""
        summary  = ""
    }
    """
    artic \\
        minion \\
        $options.args \\
        --threads $task.cpus \\
        --read-file $fastq \\
        --scheme-directory $scheme_dir \\
        --scheme-version $version \\
        $fast5 \\
        $summary \\
        $scheme \\
        $prefix
        
    echo \$(artic --version 2>&1) | sed 's/^.*artic //; s/ .*\$//' > ${software}.version.txt
    """
}

// usage: artic minion [-h] [-q] [--medaka] [--medaka-model medaka_model]
//                     [--no-longshot] [--minimap2] [--bwa]
//                     [--normalise NORMALISE] [--threads THREADS]
//                     [--scheme-directory scheme_directory]
//                     [--scheme-version scheme_version]
//                     [--max-haplotypes max_haplotypes] [--read-file read_file]
//                     [--fast5-directory FAST5_DIRECTORY]
//                     [--sequencing-summary SEQUENCING_SUMMARY]
//                     [--skip-nanopolish] [--no-indels] [--no-frameshifts]
//                     [--dry-run] [--strict]
//                     scheme sample

// The minion command can now be simplified a bit by dropping --scheme-directory and autodownload schemes, e.g.:
// artic minion --medaka --read-file test-data/*.fast[aq] scov2/V3 MT007544

// positional arguments:
//   scheme                The name of the scheme
//   sample                The name of the sample

// optional arguments:
//   -h, --help            show this help message and exit
//   -q, --quiet           Do not output warnings to stderr
//   --medaka              Use medaka instead of nanopolish for variants
//   --medaka-model medaka_model
//                         The model to use for medaka (required if using
//                         --medaka)
//   --no-longshot         Do not use Longshot for variant filtering after medaka
//   --minimap2            Use minimap2 (default)
//   --bwa                 Use bwa instead of minimap2
//   --normalise NORMALISE
//                         Normalise down to moderate coverage to save runtime
//                         (default: 100, deactivate with `--normalise 0`)
//   --threads THREADS     Number of threads (default: 8)
//   --scheme-directory scheme_directory
//                         Default scheme directory
//   --scheme-version scheme_version
//                         Primer scheme version (default: 1)
//   --max-haplotypes max_haplotypes
//                         max-haplotypes value for nanopolish
//   --read-file read_file
//                         Use alternative FASTA/FASTQ file to <sample>.fasta
//   --fast5-directory FAST5_DIRECTORY
//                         FAST5 Directory
//   --sequencing-summary SEQUENCING_SUMMARY
//                         Path to Guppy sequencing summary
//   --skip-nanopolish
//   --no-indels           Do not report InDels (uses SNP-only mode of
//                         nanopolish/medaka)
//   --no-frameshifts      Remove variants which induce frameshifts (ignored when
//                         --no-indels set)
//   --dry-run
//   --strict              Run with strict filtering of variants against primer
//                         scheme

// process articMinIONMedaka {
//     tag { sampleName }

//     label 'largecpu'

//     publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}*", mode: "copy"

//     input:
//     tuple file(fastq), file(schemeRepo)

//     output:
//     file("${sampleName}*")
//     tuple sampleName, file("${sampleName}.primertrimmed.rg.sorted.bam"), emit: ptrim
//     tuple sampleName, file("${sampleName}.sorted.bam"), emit: mapped
//     tuple sampleName, file("${sampleName}.consensus.fasta"), emit: consensus_fasta
//     tuple sampleName, file("${sampleName}.pass.vcf.gz"), emit: vcf

//     script:
//     // Make an identifier from the fastq filename
//     sampleName = fastq.getBaseName().replaceAll(~/\.fastq.*$/, '')

//     // Configure artic minion pipeline
//     minionRunConfigBuilder = []

//     if ( params.normalise ) {
//     minionRunConfigBuilder.add("--normalise ${params.normalise}")
//     }
    
//     if ( params.bwa ) {
//     minionRunConfigBuilder.add("--bwa")
//     } else {
//     minionRunConfigBuilder.add("--minimap2")
//     }

//     minionFinalConfig = minionRunConfigBuilder.join(" ")

//     """
//     artic minion --medaka \
//     ${minionFinalConfig} \
//     --threads ${task.cpus} \
//     --scheme-directory ${schemeRepo} \
//     --read-file ${fastq} \
//     ${params.scheme}/${params.schemeVersion} \
//     ${sampleName}
//     """
// }

// process articMinIONNanopolish {
//     tag { sampleName }

//     label 'largecpu'

//     publishDir "${params.outdir}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}*", mode: "copy"

//     input:
//     tuple file(fastq), file(schemeRepo), file(fast5Pass), file(seqSummary)

//     output:
//     file("${sampleName}*")
    
//     tuple sampleName, file("${sampleName}.primertrimmed.rg.sorted.bam"), emit: ptrim
//     tuple sampleName, file("${sampleName}.sorted.bam"), emit: mapped
//     tuple sampleName, file("${sampleName}.consensus.fasta"), emit: consensus_fasta
//     tuple sampleName, file("${sampleName}.pass.vcf.gz"), emit: vcf

//     script:
//     // Make an identifier from the fastq filename
//     sampleName = fastq.getBaseName().replaceAll(~/\.fastq.*$/, '')

//     // Configure artic minion pipeline
//     minionRunConfigBuilder = []

//     if ( params.normalise ) {
//     minionRunConfigBuilder.add("--normalise ${params.normalise}")
//     }
    
//     if ( params.bwa ) {
//     minionRunConfigBuilder.add("--bwa")
//     } else {
//     minionRunConfigBuilder.add("--minimap2")
//     }

//     minionFinalConfig = minionRunConfigBuilder.join(" ")


//     """
//     artic minion ${minionFinalConfig} \
//     --threads ${task.cpus} \
//     --scheme-directory ${schemeRepo} \
//     --read-file ${fastq} \
//     --fast5-directory ${fast5Pass} \
//     --sequencing-summary ${seqSummary} \
//     ${params.scheme}/${params.schemeVersion} \
//     ${sampleName}
//     """
// }
