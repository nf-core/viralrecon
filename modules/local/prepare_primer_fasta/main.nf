process PREPARE_PRIMER_FASTA {
    tag "$adapters"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:24.04' :
        'nf-core/ubuntu:24.04' }"

    input:
    path adapters

    output:
    path 'adapters.sub.fa', emit: adapters
    tuple val("${task.process}"), val('sed'), eval('echo \\$(sed --version 2>&1) | sed "s/^.*GNU sed) //; s/ .*//"'), emit: versions_sed, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    sed -r '/^[ACTGactg]+\$/ s/$args/X/g' $adapters > adapters.sub.fa
    """
}
