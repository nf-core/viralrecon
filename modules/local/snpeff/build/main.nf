process SNPEFF_BUILD {
    tag "$fasta"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snpeff:5.0--hdfd78af_1' :
        'quay.io/biocontainers/snpeff:5.0--hdfd78af_1' }"

    input:
    path fasta
    path gff

    output:
    path 'snpeff_db'   , emit: db
    path '*.config'    , emit: config
    tuple val("${task.process}"), val('snpeff'), eval('echo \\$(snpEff -version 2>&1) | cut -f 2 -d " "'), emit: versions_snpeff, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def basename = fasta.baseName
    def extension = gff.getExtension()
    if (extension == "gtf") {
        format = "gtf22"
    } else {
        format = "gff3"
    }

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
    ln -s ../../$gff genes.$extension

    cd ../../
    echo "${basename}.genome : ${basename}" > snpeff.config

    snpEff \\
        -Xmx${avail_mem}g \\
        build \\
        -config snpeff.config \\
        -dataDir ./snpeff_db \\
        -${format} \\
        $args \\
        -v \\
        ${basename}
    """
}
