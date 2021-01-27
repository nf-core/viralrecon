// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process VARSCAN_MPILEUP2CNS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    conda (params.enable_conda ? "bioconda::varscan:2.3.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/varscan:2.3.7--3"
    } else {
        container "quay.io/biocontainers/varscan:2.3.7--3"
    }

    input:
    tuple val(meta), path(mpileup)
    path  fasta

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    tuple val(meta), path("*.log"), emit: log
    path  "*.version.txt"         , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    echo "$meta.id" > sample_name.list

    varscan mpileup2cns \\
        $mpileup \\
        --vcf-sample-list sample_name.list \\
        $options.args \\
        > ${prefix}.vcf \\
        2> ${prefix}.varscan2.log
    
    sed -i.bak '/LC_ALL/d' ${prefix}.varscan2.log

    echo \$(varscan 2>&1) | sed 's/^.*VarScan v//; s/ .*\$//' > ${software}.version.txt
    """
}
