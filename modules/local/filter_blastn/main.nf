process FILTER_BLASTN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:24.04' :
        'nf-core/ubuntu:24.04' }"

    input:
    tuple val(meta), path(hits)
    path header
    path filtered_header

    output:
    tuple val(meta), path('*filter.blastn.txt')  , emit: txt
    tuple val(meta), path('*.results.blastn.txt'), emit: blast
    tuple val("${task.process}"), val('sed'), eval('echo \\$(sed --version 2>&1) | sed "s/^.*GNU sed) //; s/ .*//"'), emit: versions_sed, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_contig_length = params.min_contig_length
    def min_perc_contig_aligned = params.min_perc_contig_aligned

    """
    cat $header $hits > ${prefix}.results.blastn.txt
    awk 'BEGIN{OFS=\"\\t\";FS=\"\\t\"}{print \$0,\$7/\$17,\$7/\$16}' $hits | awk 'BEGIN{OFS=\"\\t\";FS=\"\\t\"} \$17 > ${min_contig_length} && \$19 > ${min_perc_contig_aligned} && \$1 !~ /phage/ {print \$0}' > tmp.out
    cat $filtered_header tmp.out > ${prefix}.filter.blastn.txt
    """
}
