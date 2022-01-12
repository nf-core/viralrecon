process IVAR_VARIANTS_TO_VCF {
    tag "$meta.id"

    conda (params.enable_conda ? "conda-forge::python=3.9.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(tsv)
    path  header

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val(meta), path("*.log"), emit: log
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    script:  // This script is bundled with the pipeline, in nf-core/viralrecon/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ivar_variants_to_vcf.py \\
        $tsv \\
        ${prefix}.vcf \\
        $args \\
        > ${prefix}.variant_counts.log

    cat $header ${prefix}.variant_counts.log > ${prefix}.variant_counts_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
