process KRAKEN2_BUILD {
    tag "$library"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:87fc08d11968d081f3e8a37131c1f1f6715b6542-0' :
        'quay.io/biocontainers/mulled-v2-5799ab18b5fc681e75923b2450abaa969907ec98:87fc08d11968d081f3e8a37131c1f1f6715b6542-0' }"

    input:
    val library

    output:
    path 'kraken2_db'  , emit: db
    tuple val("${task.process}"), val('kraken2'), eval('echo \\$(kraken2 --version 2>&1) | sed "s/^.*Kraken version //; s/ .*//"'), emit: versions_kraken2_build, topic: versions
    tuple val("${task.process}"), val('pigz')   , eval('pigz --version 2>&1 | sed "s/pigz //g"')                                  , emit: versions_pigz, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    """
    kraken2-build --db kraken2_db --threads $task.cpus $args  --download-taxonomy
    kraken2-build --db kraken2_db --threads $task.cpus $args2 --download-library $library
    kraken2-build --db kraken2_db --threads $task.cpus $args3 --build
    """
}
