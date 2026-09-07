process RESISTANCE_REPORT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b4/b45d221e51a26945c244afa9bd126a6757289be51b963721ce9f25f7c8662c38/data' :
        'community.wave.seqera.io/library/biopython_jinja2_pandas_python:bf9cf8457c0990de' }"

    input:
    tuple val(meta), path(sierralocal_json), path(mutation_csv), path(resistance_csv), path(nextclade_csv), path(consensus), path(annotation)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //g"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:  // This script is bundled with the pipeline, in nf-core/viralrecon/bin/
    def args = task.ext.args ?: ''
    def ivar_consensus_params = task.ext.args2 ?: '-t N/A -q N/A -m N/A -n N'
    def prefix = task.ext.prefix ?: "${meta.id}_resistance_report"

    """
    resistance_report.py \\
        --sierralocal_json $sierralocal_json \\
        --mutation_csv $mutation_csv \\
        --resistance_csv $resistance_csv \\
        --nextclade_csv $nextclade_csv \\
        --consensus_fasta $consensus \\
        --gff $annotation \\
        --ivar_consensus_params "'${ivar_consensus_params}'" \\
        --output_html ${prefix}.html \\
        $args

    """
}
