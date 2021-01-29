// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process PLOT_MOSDEPTH_REGIONS {
    tag "$bed"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'mosdepth', publish_id:'') }

    conda (params.enable_conda ? "conda-forge::r-base=4.0.3 conda-forge::r-reshape2=1.4.4 conda-forge::r-optparse=1.6.6 conda-forge::r-ggplot2=3.3.3 conda-forge::r-scales=1.1.1 conda-forge::r-viridis=0.5.1 conda-forge::r-tidyverse=1.3.0 bioconda::bioconductor-biostrings=2.58.0 bioconda::bioconductor-complexheatmap=2.6.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/"
    } else {
        container "quay.io/biocontainers/"
    }
    
    input:
    path bed
    
    output:
    path '*.pdf', emit: pdf
    path '*.tsv', emit: tsv
    
    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    plot_mosdepth_regions.r \\
        --input_files ${bed.join(',')} \\
        --input_suffix ${suffix}.regions.bed.gz \\
        --output_dir ./ \\
        --output_suffix ${suffix}.regions
    """
}