// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
def options    = initOptions(params.options)

process SPADES {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:meta.id) }

    container "staphb/spades:3.15.0"
    
    input:
    tuple val(meta), path(reads)
    path  hmm

    output:
    tuple val(meta), path('*.scaffolds.fa'), emit: scaffolds
    tuple val(meta), path('*.assembly.gfa'), emit: graph
    tuple val(meta), path('*.log')         , emit: log
    path  '*.version.txt'                  , emit: version
    
    script:
    def software    = getSoftwareName(task.process)
    def prefix      = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def input_reads = meta.single_end ? "-s $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
    def custom_hmms = params.spades_hmm ? "--custom-hmms $hmm" : ""
    """
    spades.py \\
        $options.args \\
        --threads $task.cpus \\
        $custom_hmms \\
        $input_reads \\
        -o ./

    mv scaffolds.fasta ${prefix}.scaffolds.fa
    mv assembly_graph_with_scaffolds.gfa ${prefix}.assembly.gfa
    
    echo \$(spades.py --version 2>&1) | sed 's/^.*SPAdes genome assembler v//; s/ .*\$//' > ${software}.version.txt
    """
}

