// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process IVAR_VARIANTS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? 'bioconda::ivar=1.3' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/ivar:1.3--h089eab3_0"
    } else {
        container "quay.io/biocontainers/ivar:1.3--h089eab3_0"
    }

    input:
    tuple val(meta), path(mpileup)
    path  fasta
    path  gff

    output:
    tuple val(meta), path('*.tsv'), emit: tsv
    path  '*.version.txt'         , emit: version
    
    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def features = params.gff ? "-g $gff" : ""
    """
    cat $mpileup \\
    | ivar \\
        variants \\
        $options.args \\
        -r $fasta \\
        $features \\
        -p $prefix

    echo \$(ivar version 2>&1) | sed 's/^.*iVar version //; s/ .*\$//' > ${software}.version.txt
    """
}