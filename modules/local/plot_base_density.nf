// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process PLOT_BASE_DENSITY {
    tag "$fasta"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
                saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }
    // library(optparse)
    // library(ggplot2)
    // library(scales)
    // library(reshape2)
    // library(Biostrings)

    input:
    tuple val(meta), path(fasta)
    
    output:
    tuple val(meta), path('*.pdf'), emit: pdf
    tuple val(meta), path('*.tsv'), emit: tsv

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    plot_base_density.r \\
        --fasta_files $fasta \\
        --prefixes $prefix \\
        --output_dir ./
    """
}



