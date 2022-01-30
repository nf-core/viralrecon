process RENAME_CONSENSUS {
    tag "$meta.id"

    conda (params.enable_conda ? "conda-forge::python=3.9.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fa"), emit: consensus

    script:  // This script is bundled with the pipeline, in nf-core/viralrecon/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"


    """
    sed "s/>/>${meta.id} /g" $fasta > ${prefix}.fa


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
