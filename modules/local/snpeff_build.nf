process SNPEFF_BUILD {
    tag "$fasta"
    label 'process_low'

    conda "bioconda::snpeff=5.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snpeff:5.0--hdfd78af_1' :
        'quay.io/biocontainers/snpeff:5.0--hdfd78af_1' }"

    input:
    path fasta
    path gff

    output:
    path 'snpeff_db'   , emit: db
    path '*.config'    , emit: config
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def basename = fasta.baseName

    def avail_mem = 4
    if (!task.memory) {
        log.info '[snpEff] Available memory not known - defaulting to 4GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    mkdir -p snpeff_db/genomes/
    cd snpeff_db/genomes/
    ln -s ../../$fasta ${basename}.fa

    cd ../../
    mkdir -p snpeff_db/${basename}/
    cd snpeff_db/${basename}/
    ln -s ../../$gff genes.gff

    cd ../../
    echo "${basename}.genome : ${basename}" > snpeff.config

    snpEff \\
        -Xmx${avail_mem}g \\
        build \\
        -config snpeff.config \\
        -dataDir ./snpeff_db \\
        -gff3 \\
        -v \\
        ${basename}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(echo \$(snpEff -version 2>&1) | cut -f 2 -d ' ')
    END_VERSIONS
    """
}
