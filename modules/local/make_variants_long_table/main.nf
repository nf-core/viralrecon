process MAKE_VARIANTS_LONG_TABLE {

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5e/5ee6e81aff2205d76ad8755d2181f8ea1dd747daa92aaf9dcba943b69aa9f458/data' :
        'community.wave.seqera.io/library/matplotlib_pandas_python_r-sys_pruned:23244d66110fcdf2' }"

    input:
    path bcftools_query, stageAs: "bcftools_query/*"
    path snpsift, stageAs: "snpsift/*"
    path pangolin, stageAs: "pangolin/*"

    output:
    path "*.csv"       , emit: csv
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //g"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:  // This script is bundled with the pipeline, in nf-core/viralrecon/bin/
    def args = task.ext.args ?: ''
    """
    make_variants_long_table.py \\
        --bcftools_query_dir ./bcftools_query \\
        --snpsift_dir ./snpsift \\
        --pangolin_dir ./pangolin \\
        $args
    """
}
