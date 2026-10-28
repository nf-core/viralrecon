process MULTIQC {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c8/c8e346f4f6080eadf1253505e6ff09ef004454fc18e8d672006fd7b222cc412e/data'
        : 'community.wave.seqera.io/library/multiqc:1.35--c17fb751507e9dfc'}"

    input:
    tuple val(meta), path(multiqc_files, stageAs: "?/*"), path(multiqc_config, stageAs: "?/*"), path(multiqc_logo), path(replace_names), path(sample_names), path(extra_multiqc_config)

    output:
    path "*multiqc_report.html", emit: report
    tuple val(meta), path("*_data"), emit: data
    tuple val(meta), path("*_plots"), emit: plots, optional: true
    tuple val(meta), path("*.csv"), emit: csv_reports, optional: true
    // MultiQC should not push its versions to the `versions` topic. Its input depends on the versions topic to be resolved thus outputting to the topic will let the pipeline hang forever
    tuple val("${task.process}"), val('multiqc'), eval('multiqc --version | sed "s/.* //g"'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ? "--filename ${task.ext.prefix}.html" : ''
    def config = multiqc_config ? multiqc_config instanceof List ? "--config ${multiqc_config.join(' --config ')}" : "--config ${multiqc_config}" : ""
    def extra_config = extra_multiqc_config ? "--config ${extra_multiqc_config}" : ''
    def logo = multiqc_logo ? "--cl-config 'custom_logo: \"${multiqc_logo}\"'" : ''
    def replace = replace_names ? "--replace-names ${replace_names}" : ''
    def samples = sample_names ? "--sample-names ${sample_names}" : ''
    def platform = params.platform
    """
    multiqc \\
        --force \\
        ${args} \\
        ${config} \\
        ${prefix} \\
        ${logo} \\
        ${replace} \\
        ${samples} \\
        ${extra_config} \\
        .

    ## Parse YAML files dumped by MultiQC to obtain metrics
    multiqc_to_custom_csv.py --platform $platform

    ## Manually remove files that we don't want in the report
    if grep -q ">skip_assembly<" multiqc_report.html; then
        rm -f *assembly_metrics_mqc.csv
    fi

    if grep -q ">skip_variants<" multiqc_report.html; then
        rm -f *variants_metrics_mqc.csv
    fi

    if [ "$platform" = "illumina" ]; then
        rm -f variants/report.tsv
    fi

    if [ "$platform" = "nanopore" ]; then
        rm -rf quast
    fi

    ## Run MultiQC a second time
        multiqc \\
        --force \\
        ${args2} \\
        ${config} \\
        ${prefix} \\
        ${extra_config} \\
        ${logo} \\
        ${replace} \\
        ${samples} \\
        .
    """

    stub:
    """
    mkdir multiqc_data
    touch multiqc_data/.stub
    mkdir multiqc_plots
    touch multiqc_plots/.stub
    touch multiqc_report.html
    """
}
