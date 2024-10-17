process MULTIQC {
    label 'process_medium'

    conda "bioconda::multiqc=1.19"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.19--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.19--pyhdfd78af_0' }"

    input:
    path  multiqc_files, stageAs: "?/*"
    path(multiqc_config)
    path(extra_multiqc_config)
    path(multiqc_logo)
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
    path ('freyja_demix/*')

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
    def config = multiqc_config ? "--config $multiqc_config" : ''
    def extra_config = extra_multiqc_config ? "--config $extra_multiqc_config" : ''
    def logo = multiqc_logo ? /--cl-config 'custom_logo: "${multiqc_logo}"'/ : ''

    """
    ## Run MultiQC once to parse tool logs
    multiqc -f $args $config $extra_config $logo .

    ## Parse YAML files dumped by MultiQC to obtain metrics
    multiqc_to_custom_csv.py --platform nanopore

    ## Manually remove files that we don't want in the report
    rm -rf quast

    ## Run MultiQC a second time
    multiqc -f $args -e general_stats --ignore *nextclade_clade_mqc.tsv $config $extra_config $logo .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}
