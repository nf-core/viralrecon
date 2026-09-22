process EXTRACT_CONSENSUS_REGIONS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b4/b45d221e51a26945c244afa9bd126a6757289be51b963721ce9f25f7c8662c38/data' :
        'community.wave.seqera.io/library/biopython_jinja2_pandas_python:bf9cf8457c0990de' }"

    input:
    tuple val(meta), path(consensus), path(annotation)

    output:
    tuple val(meta), path("*.all_consensus.fa"), emit: consensus_regions
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //g"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '--interest_genes "PR,RT;IN"'
    def prefix   = task.ext.prefix ?: "${meta.id}"

    """
    extract_consensus_regions.py \
        --consensus_fasta $consensus \
        --gff $annotation \
        $args \
        --output "${prefix}.all_consensus.fa"
    """
}
