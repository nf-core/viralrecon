// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

process PLOT_MOSDEPTH_REGIONS {
    tag "$bed"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'mosdepth', publish_id:'') }

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }
    // library(optparse)
    // library(ggplot2)
    // library(scales)
    // library(ComplexHeatmap)
    // library(viridis)
    // library(tidyverse)

    input:
    path bed
    
    output:
    path  '*.pdf', emit: pdf
    path  '*.tsv', emit: tsv
    
    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    plot_mosdepth_regions.r \\
        --input_files ${bed.join(',')} \\
        --input_suffix ${suffix}.regions.bed.gz \\
        --output_dir ./ \\
        --output_suffix ${suffix}.regions
    """
}

//     plot_mosdepth_regions.r \\
//         --input_files ${prefix}.regions.bed.gz \\
//         --input_suffix ${plot_suffix}.regions.bed.gz \\
//         --output_dir ./ \\
//         --output_suffix ${plot_suffix}.regions

//         plot_mosdepth_regions.r \\
//             --input_files ${bed.join(',')} \\
//             --input_suffix ${suffix}.regions.bed.gz \\
//             --output_dir ./ \\
//             --output_suffix ${suffix}.regions
