process COLLAPSE_PRIMERS {
    tag "$bed"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.13' :
        'quay.io/biocontainers/python:3.13' }"

    input:
    path bed
    val left_suffix
    val right_suffix

    output:
    path '*.bed'       , emit: bed
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //g"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/viralrecon/bin/
    """
    collapse_primer_bed.py \\
        --left_primer_suffix $left_suffix \\
        --right_primer_suffix $right_suffix \\
        $bed \\
        ${bed.baseName}.collapsed.bed
    """
}
