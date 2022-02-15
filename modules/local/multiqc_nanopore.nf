process MULTIQC {
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::multiqc=1.11" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.11--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0' }"

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
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def custom_config = multiqc_custom_config ? "--config $multiqc_custom_config" : ''
    """
    ## Run MultiQC once to parse tool logs
    multiqc -f $args $custom_config .

    ## Parse YAML files dumped by MultiQC to obtain metrics
    multiqc_to_custom_csv.py --platform nanopore

    ## Manually remove files that we don't want in the report
    rm -rf quast

    ## Run MultiQC a second time
    multiqc -f $args -e general_stats --ignore *nextclade_clade_mqc.tsv $custom_config .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
