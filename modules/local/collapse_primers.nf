// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

process COLLAPSE_PRIMERS {
    tag "$bed"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'primers', publish_id:'') }

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    input:
    path bed
    val  left_suffix
    val  right_suffix

    output:
    path  '*.bed', emit: bed
    
    script:
    """
    collapse_primer_bed.py \\
        --left_primer_suffix $left_suffix \\
        --right_primer_suffix $right_suffix \\
        $bed \\
        ${bed.baseName}.collapsed.bed
    """
}