// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ABACAS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }
        
    conda (params.enable_conda ? 'hcc::abacas=1.3.1' : null)
    container 'biocontainers/abacas:v1.3.1-5-deb_cv1'
    
    input:
    tuple val(meta), path(scaffold)
    path  fasta
    
    output:
    tuple val(meta), path('*.abacas*'), emit: results
    path '*.version.txt'              , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    abacas.pl \\
        -r $fasta \\
        -q $scaffold \\
        $options.args \\
        -o ${prefix}.abacas
    
    mv nucmer.delta ${prefix}.abacas.nucmer.delta
    mv nucmer.filtered.delta ${prefix}.abacas.nucmer.filtered.delta
    mv nucmer.tiling ${prefix}.abacas.nucmer.tiling
    mv unused_contigs.out ${prefix}.abacas.unused.contigs.out

    echo \$(abacas.pl -v 2>&1) | sed 's/^.*blastn: //; s/ .*\$//' > ${software}.version.txt
    """
}