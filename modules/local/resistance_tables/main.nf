process RESISTANCE_TABLES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5e/5ee6e81aff2205d76ad8755d2181f8ea1dd747daa92aaf9dcba943b69aa9f458/data' :
        'community.wave.seqera.io/library/matplotlib_pandas_python_r-sys_pruned:23244d66110fcdf2' }"

    input:
    tuple val(meta), path(json), path(codfreq)

    output:
    tuple val(meta), path("*_mutation_table.csv")       , emit: mutation_csv
    tuple val(meta), path("*_mutation_table_short.csv") , emit: mutation_short_csv
    tuple val(meta), path("*_resistance_table.csv")     , emit: resistance_csv
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //g"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:  // This script is bundled with the pipeline, in nf-core/viralrecon/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    resistance_tables.py \\
        --sierralocal_file $json \\
        --codfreq_file $codfreq \\
        --sample_name ${meta.id} \\
        --output_mutation_file ${prefix}_mutation_table.csv \\
        --output_resistance_file ${prefix}_resistance_table.csv \\
        --output_mutation_short ${prefix}_mutation_table_short.csv \\
        $args
    """
}
