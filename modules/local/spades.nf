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
    val   coronaspades

    output:
    tuple val(meta), path('*.scaffolds.fa')    , emit: scaffolds
    tuple val(meta), path('*.contigs.fa')      , emit: contigs
    tuple val(meta), path('*.assembly.gfa')    , emit: gfa
    tuple val(meta), path('*.assembly.fastg')  , emit: fastg
    tuple val(meta), path('*.log')             , emit: log
    tuple val(meta), path('*.gene_clusters.fa'), optional:true, emit: gene_clusters
    path  '*.version.txt'                      , emit: version
    
    script:
    def software    = getSoftwareName(task.process)
    def prefix      = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def input_reads = meta.single_end ? "-s $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
    def custom_hmms = params.spades_hmm ? "--custom-hmms $hmm" : ""
    def command     = coronaspades ? "coronaspades.py" : "spades.py"
    """
    $command \\
        $options.args \\
        --threads $task.cpus \\
        $custom_hmms \\
        $input_reads \\
        -o ./

    mv scaffolds.fasta ${prefix}.scaffolds.fa
    mv assembly_graph_with_scaffolds.gfa ${prefix}.assembly.gfa
    mv contigs.fasta ${prefix}.contigs.fa
    mv assembly_graph.fastg ${prefix}.assembly.fastg
    
    if [ -f gene_clusters.fasta ]; then
        mv gene_clusters.fasta ${prefix}.gene_clusters.fa
    fi
    
    echo \$(spades.py --version 2>&1) | sed 's/^.*SPAdes genome assembler v//; s/ .*\$//' > ${software}.version.txt
    """
}
