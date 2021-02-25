// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process MULTIQC {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::multiqc=1.9" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/multiqc:1.9--pyh9f0ad1d_0"
    } else {
        container "quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
    }

    input:
    path 'multiqc_config.yaml'
    path multiqc_custom_config
    path software_versions
    path workflow_summary
    path fail_mapping_summary
    path ('fastqc/*')
    path ('fastp/*')
    path ('fastp/*')
    path ('kraken2/*')
    path ('bowtie2/*')
    path ('bowtie2/*')
    path ('bowtie2/*')
    path ('bowtie2/*')
    path ('ivar_trim/*')
    path ('ivar_trim/*')
    path ('ivar_trim/*')
    path ('ivar_trim/*')
    path ('picard_markduplicates/*')
    path ('picard_markduplicates/*')
    path ('picard_markduplicates/*')
    path ('picard_markduplicates/*')
    path ('picard_metrics/*')
    path ('picard_metrics/*')
    path ('mosdepth/*')
    path ('variants_ivar/*')
    path ('variants_ivar/*')
    path ('variants_ivar/*')
    path ('variants_ivar/*')
    path ('variants_bcftools/*')
    path ('variants_bcftools/*')
    path ('variants_bcftools/*')
    path ('cutadapt/*')
    path ('cutadapt/*')
    path ('assembly_spades/*')
    path ('assembly_unicycler/*')
    path ('assembly_minia/*')
        
    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots

    script:
    def software      = getSoftwareName(task.process)
    def custom_config = params.multiqc_config ? "--config $multiqc_custom_config" : ''
    """
    multiqc -f $options.args $custom_config .
    multiqc_to_custom_tsv_illumina.py
    multiqc -f $options.args $custom_config .
    """
}
