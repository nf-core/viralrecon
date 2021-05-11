// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process IVAR_VARIANTS_TO_VCF {
    tag "$meta.id"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    input:
    tuple val(meta), path(tsv)
    path  header

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val(meta), path("*.log"), emit: log
    tuple val(meta), path("*.tsv"), emit: tsv

    script:  // This script is bundled with the pipeline, in nf-core/viralrecon/bin/
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    ivar_variants_to_vcf.py \\
        $tsv \\
        ${prefix}.vcf \\
        $options.args \\
        > ${prefix}.variant_counts.log

    cat $header ${prefix}.variant_counts.log > ${prefix}.variant_counts_mqc.tsv
    """
}
