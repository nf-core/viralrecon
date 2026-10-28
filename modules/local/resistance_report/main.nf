process RESISTANCE_REPORT {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b4/b45d221e51a26945c244afa9bd126a6757289be51b963721ce9f25f7c8662c38/data' :
        'community.wave.seqera.io/library/biopython_jinja2_pandas_python:bf9cf8457c0990de' }"

    input:
    path sierralocal_json, stageAs: "sierralocal_json/*"
    path mutation_csv    , stageAs: "mutation_tables/*"
    path resistance_csv  , stageAs: "resistance_tables/*"
    path nextclade_csv   , stageAs: "nextclade_reports/*"
    path consensus       , stageAs: "consensus/*"
    path annotation      , stageAs: "gff/*"

    output:
    path("*.html")      , emit: mutation_csv
    tuple val("${task.process}"), val('python'), eval('python --version | sed "s/Python //g"'), emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:  // This script is bundled with the pipeline, in nf-core/viralrecon/bin/
    def args = task.ext.args ?: ''
    def ivar_consensus_params = task.ext.args2 ?: '-t N/A -q N/A -m N/A -n N'
    def prefix = task.ext.prefix ?: 'resistance'

    """
    resistance_report.py \\
        --sierralocal_folder ./sierralocal_json \\
        --mutation_folder ./mutation_tables \\
        --resistance_folder ./resistance_tables \\
        --nextclade_folder ./nextclade_reports \\
        --consensus_folder ./consensus \\
        --gff_folder ./gff \\
        --ivar_consensus_params "'${ivar_consensus_params}'" \\
        --output_html ${prefix}.html \\
        $args
    """
}
