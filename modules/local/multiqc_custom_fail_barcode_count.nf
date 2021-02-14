// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options = [:]

process MULTIQC_CUSTOM_FAIL_BARCODE_COUNT {
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
    val fail_barcode_count
    
    output:
    path "*.tsv"

    script:
    if (fail_barcode_count.size() > 0) {
        """
        echo "Sample\tBarcode count" > fail_barcode_count_samples_mqc.tsv
        echo "${fail_barcode_count.join('\n')}" >> fail_barcode_count_samples_mqc.tsv
        """
    } else {
        """
        touch fail_barcode_count_samples_mqc.tsv
        """
    }
}
