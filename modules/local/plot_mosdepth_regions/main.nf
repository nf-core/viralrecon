process PLOT_MOSDEPTH_REGIONS {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/5d/5d4923bc2b9c2d8c83468eddf1e32c0bfecffab8192e0e2ee69340fe267560cd/data' :
        'community.wave.seqera.io/library/bioconductor-biostrings_bioconductor-complexheatmap_r-base_r-ggplot2_pruned:6c26995a32d99713' }"

    input:
    path beds

    output:
    path '*coverage.pdf', emit: coverage_pdf
    path '*coverage.tsv', emit: coverage_tsv
    path '*heatmap.pdf' , optional:true, emit: heatmap_pdf
    path '*heatmap.tsv' , optional:true, emit: heatmap_tsv
    tuple val("${task.process}"), val('r-base'), eval('echo \\$(R --version 2>&1) | sed "s/^.*R version //; s/ .*//"'), emit: versions_r_base, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in nf-core/viralrecon/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "mosdepth"
    """
    plot_mosdepth_regions.r \\
        --input_files ${beds.join(',')} \\
        --output_dir ./ \\
        --output_suffix $prefix \\
        $args
    """
}
