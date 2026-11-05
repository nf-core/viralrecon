process BAM2CODFREQ {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.23.3--py39hdd5828d_1' :
        'biocontainers/pysam:0.23.3--py39hdd5828d_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(profile)

    output:
    tuple val(meta), path("*.codfreq") , emit: codfreq
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //g"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:  // This script is bundled with the pipeline, in nf-core/viralrecon/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    bam2codfreq.py \\
        --bam $bam \\
        --profile $profile \\
        --output ${prefix}.codfreq \\
        $args
    """
}
