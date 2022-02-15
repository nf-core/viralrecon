process ASCIIGENOME {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::asciigenome=1.16.0 bioconda::bedtools=2.30.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-093691b47d719890dc19ac0c13c4528e9776897f:27211b8c38006480d69eb1be3ef09a7bf0a49d76-0' :
        'quay.io/biocontainers/mulled-v2-093691b47d719890dc19ac0c13c4528e9776897f:27211b8c38006480d69eb1be3ef09a7bf0a49d76-0' }"

    input:
    tuple val(meta), path(bam), path(vcf)
    path fasta
    path sizes
    path gff
    path bed
    val window
    val track_height

    output:
    tuple val(meta), path("*pdf"), emit: pdf
    path "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gff_track = gff ? "$gff" : ''
    def bed_track = bed ? "$bed" : ''
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        asciigenome: \$(echo \$(ASCIIGenome -ni --version 2>&1) | sed -e "s/ASCIIGenome //g")
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
