process PLOT_BASE_DENSITY {
    tag "$fasta"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5d/5d4923bc2b9c2d8c83468eddf1e32c0bfecffab8192e0e2ee69340fe267560cd/data' :
        'community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-complexheatmap_r-base_r-ggplot2_pruned:6c26995a32d99713' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('*.pdf'), emit: pdf
    tuple val(meta), path('*.tsv'), emit: tsv
    tuple val("${task.process}"), val('r-base'), eval('echo \\$(R --version 2>&1) | sed "s/^.*R version //; s/ .*//"'), emit: versions_r_base, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/viralrecon/bin/
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plot_base_density.r \\
        --fasta_files $fasta \\
        --prefixes $prefix \\
        --output_dir ./
    """
}
