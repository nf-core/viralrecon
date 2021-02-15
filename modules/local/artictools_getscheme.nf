// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process ARTICTOOLS_GETSCHEME {
    tag "$scheme/$scheme_version"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::artic=1.2.1" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/artic:1.2.1--py_0"
    } else {
        container "quay.io/biocontainers/artic:1.2.1--py_0"
    }

    input:
    val scheme
    val scheme_version
    
    output:
    path "*.fasta"       , emit: fasta
    path "*.bed"         , emit: bed
    path "*.version.txt" , emit: version

    script:
    def software = getSoftwareName(task.process)
    """
    artic-tools \\
        get_scheme \\
        $scheme \\
        --schemeVersion $scheme_version \\
        --outDir ./
    
    echo \$(artic --version 2>&1) | sed 's/^.*artic //; s/ .*\$//' > ${software}.version.txt
    """
}
