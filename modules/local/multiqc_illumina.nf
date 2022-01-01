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
    path fail_reads_summary
    path fail_mapping_summary
    path 'amplicon_heatmap_mqc.tsv'
    path ('fastqc/*')
    path ('fastp/*')
    path ('kraken2/*')
    path ('bowtie2/*')
    path ('bowtie2/*')
    path ('ivar_trim/*')
    path ('picard_markduplicates/*')
    path ('mosdepth/*')
    path ('variants_ivar/*')
    path ('variants_ivar/*')
    path ('variants_ivar/*')
    path ('variants_ivar/*')
    path ('variants_ivar/*')
    path ('variants_ivar/*')
    path ('variants_bcftools/*')
    path ('variants_bcftools/*')
    path ('variants_bcftools/*')
    path ('variants_bcftools/*')
    path ('variants_bcftools/*')
    path ('cutadapt/*')
    path ('assembly_spades/*')
    path ('assembly_unicycler/*')
    path ('assembly_minia/*')

    output:
    path "*multiqc_report.html"     , emit: report
    path "*_data"                   , emit: data
    path "*variants_metrics_mqc.csv", optional:true, emit: csv_variants
    path "*assembly_metrics_mqc.csv", optional:true, emit: csv_assembly
    path "*_plots"                  , optional:true, emit: plots

    script:
    def software      = getSoftwareName(task.process)
    def custom_config = params.multiqc_config ? "--config $multiqc_custom_config" : ''
    """
    ## Run MultiQC once to parse tool logs
    multiqc -f $options.args $custom_config .

    ## Parse YAML files dumped by MultiQC to obtain metrics
    multiqc_to_custom_csv.py --platform illumina

    ## Manually remove files that we don't want in the report
    if grep -q ">skip_assembly<" workflow_summary_mqc.yaml; then
        rm -f *assembly_metrics_mqc.csv
    fi

    if grep -q ">skip_variants<" workflow_summary_mqc.yaml; then
        rm -f *variants_metrics_mqc.csv
    fi

    rm -f variants_ivar/report.tsv
    rm -f variants_bcftools/report.tsv

    ## Run MultiQC a second time
    multiqc -f $options.args -e general_stats --ignore *nextclade_clade_mqc.tsv $custom_config .
    """
}
