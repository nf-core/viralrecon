process IVAR_VARIANTS_TO_VCF {
    tag "$meta.id"

    conda "conda-forge::python=3.9.5 conda-forge::matplotlib=3.5.1 conda-forge::pandas=1.3.5 conda-forge::r-sys=3.4 conda-forge::regex=2021.11.10 conda-forge::scipy=1.7.3 conda-forge::biopython=1.79"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-ff46c3f421ca930fcc54e67ab61c8e1bcbddfe22:1ad3da14f705eb0cdff6b5a44fea4909307524b4-0' :
        'quay.io/biocontainers/mulled-v2-ff46c3f421ca930fcc54e67ab61c8e1bcbddfe22:1ad3da14f705eb0cdff6b5a44fea4909307524b4-0' }"

    input:
    tuple val(meta), path(tsv)
    path fasta
    path header

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val(meta), path("*.log"), emit: log
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:  // This script is bundled with the pipeline, in nf-core/viralrecon/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ivar_variants_to_vcf.py \\
        $tsv \\
        ${prefix}.vcf \\
        --fasta $fasta \\
        $args \\
        > ${prefix}.variant_counts.log

    cat $header ${prefix}.variant_counts.log > ${prefix}.variant_counts_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
