process MAKE_BED_MASK {
    tag "$meta.id"

    conda (params.enable_conda ? "conda-forge::python=3.9.5 conda-forge::gawk=5.1.0 bioconda::samtools=1.12" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plasmidid:1.6.5--hdfd78af_0' :
        'quay.io/biocontainers/plasmidid:1.6.5--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(vcf)
    path fasta
    val save_mpileup

    output:
    tuple val(meta), path("*.bed")    , emit: bed
    tuple val(meta), path("*.mpileup"), optional:true, emit: mpileup
    path "versions.yml"               , emit: versions

    script:  // This script is bundled with the pipeline, in nf-core/viralrecon/bin/
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: 10
    def prefix = task.ext.prefix ?: "${meta.id}"
    def mpileup = save_mpileup ? "| tee ${prefix}.mpileup" : ""
    """
    samtools \\
        mpileup \\
        $args \\
        --reference $fasta \\
        $bam \\
        $mpileup  \\
        | awk -v OFS='\\t' '{print \$1, \$2-1, \$2, \$4}' | awk '\$4 < $args2' > lowcov_positions.txt

    make_bed_mask.py \\
        $vcf \\
        lowcov_positions.txt \\
        ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
