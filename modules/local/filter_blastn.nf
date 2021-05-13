// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FILTER_BLASTN {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'blastn', meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }

    input:
    tuple val(meta), path(hits)
    path header

    output:
    tuple val(meta), path('*.txt'), emit: txt

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    awk 'BEGIN{OFS=\"\\t\";FS=\"\\t\"}{print \$0,\$5/\$15,\$5/\$14}' $hits | awk 'BEGIN{OFS=\"\\t\";FS=\"\\t\"} \$15 > 200 && \$17 > 0.7 && \$1 !~ /phage/ {print \$0}' > tmp.out
    cat $header tmp.out > ${prefix}.txt
    """
}
