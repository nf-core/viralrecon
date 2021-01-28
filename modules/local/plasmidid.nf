// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process PLASMIDID {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }
        
    conda (params.enable_conda ? 'bioconda::plasmidid=1.6.3' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/plasmidid:1.6.3--1'
    } else {
        container 'quay.io/biocontainers/plasmidid:1.6.3--1'
    }

    input:
    tuple val(meta), path(scaffold)
    path  fasta
    
    output:
    tuple val(meta), path("${prefix}"), emit: results
    path '*.version.txt'              , emit: version

    script:
    def software = getSoftwareName(task.process)
    prefix       = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    plasmidID \\
        -d $fasta \\
        -s $prefix \\
        -c $scaffold \\
        $options.args \\
        -o .
    
    mv NO_GROUP/$prefix ./$prefix

    echo \$(plasmidID --version 2>&1) > ${software}.version.txt
    """
}