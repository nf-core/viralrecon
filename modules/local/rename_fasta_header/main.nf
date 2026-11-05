process RENAME_FASTA_HEADER {
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:24.04' :
        'nf-core/ubuntu:24.04' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fa"), emit: fasta
    tuple val("${task.process}"), val('sed'), eval('echo \\$(sed --version 2>&1) | sed "s/^.*GNU sed) //; s/ .*//"'), emit: versions_sed, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    sed "s/>/>${meta.id} /g" $fasta > ${prefix}.fa
    """
}
