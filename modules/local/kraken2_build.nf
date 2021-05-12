// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process KRAKEN2_BUILD {
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? 'bioconda::kraken2=2.1.1' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/kraken2:2.1.1--pl526hc9558a2_0'
    } else {
        container 'quay.io/biocontainers/kraken2:2.1.1--pl526hc9558a2_0'
    }

    input:
    val library

    output:
    path 'kraken2_db'   , emit: db
    path '*.version.txt', emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    kraken2-build --db kraken2_db --threads $task.cpus $options.args  --download-taxonomy
    kraken2-build --db kraken2_db --threads $task.cpus $options.args2 --download-library $library
    kraken2-build --db kraken2_db --threads $task.cpus $options.args3 --build

    echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//' > ${software}.version.txt
    """
}
