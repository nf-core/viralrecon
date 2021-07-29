// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MULTIQC {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::multiqc=1.11" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/multiqc:1.11--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0"
    }

    input:
    path 'multiqc_config.yaml'
    path multiqc_custom_config
    path software_versions
    path workflow_summary
    path fail_barcodes_no_sample
    path fail_no_barcode_samples
    path fail_barcode_count_samples
    path fail_guppyplex_count_samples
    path 'amplicon_heatmap_mqc.tsv'
    path ('pycoqc/*')
    path ('artic_minion/*')
    path ('samtools_stats/*')
    path ('bcftools_stats/*')
    path ('mosdepth/*')
    path ('quast/*')
    path ('snpeff/*')
    path pangolin_lineage
    path nextclade_clade

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*.csv"               , optional:true, emit: csv
    path "*_plots"             , optional:true, emit: plots

    script:
    def software      = getSoftwareName(task.process)
    def custom_config = params.multiqc_config ? "--config $multiqc_custom_config" : ''
    """
    ## Run MultiQC once to parse tool logs
    multiqc -f $options.args $custom_config .

    ## Parse YAML files dumped by MultiQC to obtain metrics
    multiqc_to_custom_csv.py --platform nanopore

    ## Manually remove files that we don't want in the report
    rm -rf quast

    ## Run MultiQC a second time
    multiqc -f $options.args -e general_stats --ignore *nextclade_clade_mqc.tsv $custom_config .
    """
}
