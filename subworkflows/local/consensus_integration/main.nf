//
// This subworkflow merges consensus sequences from mapping 
// and de novo (ABACAS) strategies using prior evolutionary information
// to generate a single, more complete consensus sequence
// with improved coverage via PriorCons. 
//

include { CAT_CAT     } from '../../../modules/nf-core/cat/cat/main'
include { MAFFT_ALIGN } from '../../../modules/nf-core/mafft/align/main' 
// TODO uncomment the following once priorcons is a nf-core module
// include { PRIORCONS } from '../../../modules/nf-core/'

workflow CONSENSUS_INTEGRATION {

    take:
    ch_fasta_reference // channel: /path/to/genome.fasta
    ch_fasta_mapping   // channel: [ val(meta), path(fasta) ]
    ch_fasta_assembly  // channel: [ val(meta), path(fasta) ]

    main:

    ch_versions = Channel.empty()
    
    //
    // Cat together all sample rerlated fasta
    //
    ch_fasta_combined = ch_fasta_mapping
        .join(ch_fasta_assembly)
        .combine(ch_fasta_reference)
        .map { meta, mapping_fasta, assembly_fasta, reference_fasta ->
            [meta, [reference_fasta, mapping_fasta, assembly_fasta]]
        }

    // TODO add CAT_CAT to modules_illumina.config
    /*
    CAT_CAT(ch_fasta_combined)
    ch_aignment = CAT_CAT.out.file_out
    ch_versions = ch_versions.mix(CAT_CAT.out.versions.first())
    */
    // TODO add the remaining processes 
    /*
    SAMTOOLS_SORT ( ch_bam )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    // TODO nf-core: edit emitted channels
    bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
    */
}
