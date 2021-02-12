// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SNPEFF_BUILD {
    tag "$fasta"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? 'bioconda::snpeff=5.0' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/snpeff:5.0--0'
    } else {
        container 'quay.io/biocontainers/snpeff:5.0--0'
    }

    input:
    path fasta
    path gff

    output:
    path 'snpeff_db'    , emit: db
    path '*.config'     , emit: config
    path '*.version.txt', emit: version

    script:
    def software = getSoftwareName(task.process)
    def basename = fasta.baseName
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

    snpEff build -config snpeff.config -dataDir ./snpeff_db -gff3 -v ${basename}

    echo \$(snpEff -version 2>&1) | sed 's/^.*SnpEff //; s/ .*\$//' > ${software}.version.txt
    """
}