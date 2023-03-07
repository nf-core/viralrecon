process CUTADAPT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::cutadapt=4.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cutadapt:4.2--py39hbf8eff0_0' :
        'quay.io/biocontainers/cutadapt:4.2--py39hbf8eff0_0' }"

    input:
    tuple val(meta), path(reads)
    path adapters

    output:
    tuple val(meta), path('*.fastq.gz'), emit: reads
    tuple val(meta), path('*.log')     , emit: log
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def paired = meta.single_end ? "-a file:adapters.sub.fa" : "-a file:adapters.sub.fa -A file:adapters.sub.fa"
    def trimmed = meta.single_end ? "-o ${prefix}.fastq.gz" : "-o ${prefix}_1.fastq.gz -p ${prefix}_2.fastq.gz"
    """
    sed -r '/^[ACTGactg]+\$/ s/\$/X/g' $adapters > adapters.sub.fa

    cutadapt \\
        --cores $task.cpus \\
        $args \\
        $paired \\
        $trimmed \\
        $reads \\
        > ${prefix}.cutadapt.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """
}
