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