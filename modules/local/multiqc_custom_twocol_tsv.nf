// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options = [:]

process MULTIQC_CUSTOM_TWOCOL_TSV {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img"
    } else {
        container "biocontainers/biocontainers:v1.2.0_cv1"
    }
    
    input:
    val tsv_data
    val col1_name
    val col2_name
    val out_prefix

    output:
    path "*.tsv"

    script:
    if (tsv_data.size() > 0) {
        """
        echo "${col1_name}\t${col2_name}" > ${out_prefix}_mqc.tsv
        echo "${tsv_data.join('\n')}" >> ${out_prefix}_mqc.tsv
        """
    } else {
        """
        touch ${out_prefix}_mqc.tsv
        """
    }
}
