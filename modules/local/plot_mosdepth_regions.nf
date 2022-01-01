// Import generic module functions
include { initOptions; saveFiles } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PLOT_MOSDEPTH_REGIONS {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'mosdepth', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::r-base=4.0.3 conda-forge::r-reshape2=1.4.4 conda-forge::r-optparse=1.6.6 conda-forge::r-ggplot2=3.3.3 conda-forge::r-scales=1.1.1 conda-forge::r-viridis=0.5.1 conda-forge::r-tidyverse=1.3.0 bioconda::bioconductor-biostrings=2.58.0 bioconda::bioconductor-complexheatmap=2.6.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-ad9dd5f398966bf899ae05f8e7c54d0fb10cdfa7:05678da05b8e5a7a5130e90a9f9a6c585b965afa-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-ad9dd5f398966bf899ae05f8e7c54d0fb10cdfa7:05678da05b8e5a7a5130e90a9f9a6c585b965afa-0"
    }

    input:
    path beds

    output:
    path '*coverage.pdf', emit: coverage_pdf
    path '*coverage.tsv', emit: coverage_tsv
    path '*heatmap.pdf' , optional:true, emit: heatmap_pdf
    path '*heatmap.tsv' , optional:true, emit: heatmap_tsv

    script:
    def prefix = options.suffix ?: "mosdepth"
    """
    plot_mosdepth_regions.r \\
        --input_files ${beds.join(',')} \\
        --output_dir ./ \\
        --output_suffix $prefix \\
        $options.args
    """
}
