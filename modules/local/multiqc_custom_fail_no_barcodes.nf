// Import generic module functions
include { saveFiles; getSoftwareName } from './functions'

params.options = [:]

process MULTIQC_CUSTOM_FAIL_NO_BARCODES {
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
    val fail_no_barcodes
    

    output:
    path "*.tsv"

    script:
    if (fail_no_barcodes.size() > 0) {
        """
        echo "Sample\tMissing barcode" > fail_no_barcode_samples_mqc.tsv
        echo "${fail_no_barcodes.join('\n')}" >> fail_no_barcode_samples_mqc.tsv
        """
    } else {
        """
        touch fail_no_barcode_samples_mqc.tsv
        """
    }
}
