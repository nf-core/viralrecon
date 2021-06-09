// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options = [:]

process MULTIQC_CUSTOM_ONECOL_TXT {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }

    input:
    val data
    val out_prefix

    output:
    path "*.txt"

    script:
    if (data.size() > 0) {
        """
        echo "${data.join('\n')}" >> ${out_prefix}_mqc.txt
        """
    } else {
        """
        touch ${out_prefix}_mqc.txt
        """
    }
}
