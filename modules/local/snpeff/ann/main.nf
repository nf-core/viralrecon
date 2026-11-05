process SNPEFF_ANN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snpeff:5.0--hdfd78af_1' :
        'quay.io/biocontainers/snpeff:5.0--hdfd78af_1' }"

    input:
    tuple val(meta), path(vcf)
    path  db
    path  config
    path  fasta

    output:
    tuple val(meta), path("*.vcf")      , emit: vcf
    tuple val(meta), path("*.csv")      , emit: csv
    tuple val(meta), path("*.genes.txt"), emit: txt
    tuple val(meta), path("*.html")     , emit: html
    tuple val("${task.process}"), val('snpeff'), eval('echo \\$(snpEff -version 2>&1) | cut -f 2 -d " "'), emit: versions_snpeff, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def avail_mem = 4
    if (!task.memory) {
        log.info '[snpEff] Available memory not known - defaulting to 4GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    snpEff \\
        -Xmx${avail_mem}g \\
        ${fasta.baseName} \\
        -config $config \\
        -dataDir $db \\
        $args \\
        $vcf \\
        -csvStats ${prefix}.snpeff.csv \\
        > ${prefix}.snpeff.vcf
    mv snpEff_summary.html ${prefix}.snpeff.summary.html
    """
}
