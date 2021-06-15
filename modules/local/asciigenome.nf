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

    conda (params.enable_conda ? "bioconda::asciigenome=1.16.0 bioconda::bedtools=2.30.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-093691b47d719890dc19ac0c13c4528e9776897f:27211b8c38006480d69eb1be3ef09a7bf0a49d76-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-093691b47d719890dc19ac0c13c4528e9776897f:27211b8c38006480d69eb1be3ef09a7bf0a49d76-0"
    }

    input:
    tuple val(meta), path(bam), path(vcf)
    path  fasta
    path  sizes
    path  gff
    path  bed
    val   window
    val   track_height

    output:
    tuple val(meta), path("*pdf"), emit: pdf
    path  "*.version.txt"        , emit: version

    script:
    def software   = getSoftwareName(task.process)
    def prefix     = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def gff_track  = gff ? "$gff" : ''
    def bed_track  = bed ? "$bed" : ''
    def paired_end = meta.single_end ? '' : '&& readsAsPairs -on'
    """
    zcat $vcf \\
        | grep -v '#' \\
        | awk -v FS='\t' -v OFS='\t' '{print \$1, (\$2-1), (\$2)}' \\
        > variants.bed

    bedtools \\
        slop \\
        -i variants.bed \\
        -g $sizes \\
        -b $window \\
        > variants.slop.bed

    ASCIIGenome \\
        -ni \\
        -x "trackHeight 0 bam#1 && trackHeight $track_height bam@2 $paired_end && filterVariantReads && save ${prefix}.%r.pdf" \\
        --batchFile variants.slop.bed \\
        --fasta $fasta \\
        $bam \\
        $vcf \\
        $bed_track \\
        $gff_track \\
        > /dev/null

    echo \$(ASCIIGenome -ni --version 2>&1) | sed -e "s/ASCIIGenome //g" > ${software}.version.txt
    """
}
