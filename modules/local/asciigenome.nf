// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ASCIIGENOME {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }
        
    conda (params.enable_conda ? "bioconda::asciigenome=1.16.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/asciigenome:1.16.0--0"
    } else {
        container "quay.io/biocontainers/asciigenome:1.16.0--0"
    }
        
    input:
    tuple val(meta), path(bam)
    path  fasta
    path  vcf
    path  bed
    val   use_bed
    val   window
    val   track_height

    // output:
    // tuple val(meta), path("*_metrics"), emit: metrics
    // path  "*.version.txt"             , emit: version

    script:
    def software   = getSoftwareName(task.process)
    def prefix     = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def bed_track  = use_bed ? "$bed" : ""
    def paired_end = meta.single_end ? "" : "&& readsAsPairs -on"
    """
    zcat $vcf \\
        | grep -v '#' \\
        | awk -v FS='\t' -v OFS='\t' '{print \$1, (\$2-$windown), (\$2+$window)}' \\
        > variants.bed

    ASCIIGenome \\
        -ni \\
        -x "trackHeight $track_height $paired_end && save ${meta.id}.%r.pdf" \\
        --batchFile variants.bed \\
        --fasta $fasta \\
        $bam \\
        $vcf \\
        $bed_track \\
        > /dev/null

    """
}

    // echo \$(picard CollectWgsMetrics --version 2>&1) | grep -o 'Version.*' | cut -f2- -d: > ${software}.version.txt
