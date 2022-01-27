//
// Consensus calling with BCFTools and downstream processing QC
//

include { BCFTOOLS_FILTER    } from '../../modules/nf-core/modules/bcftools/filter/main'
include { BEDTOOLS_MERGE     } from '../../modules/nf-core/modules/bedtools/merge/main'
include { BEDTOOLS_MASKFASTA } from '../../modules/nf-core/modules/bedtools/maskfasta/main'
include { BCFTOOLS_CONSENSUS } from '../../modules/nf-core/modules/bcftools/consensus/main'
include { MAKE_BED_MASK      } from '../../modules/local/make_bed_mask'
include { CONSENSUS_QC       } from './consensus_qc'

workflow CONSENSUS_BCFTOOLS {
    take:
    bam          // channel: [ val(meta), [ bam ] ]
    vcf          // channel: [ val(meta), [ vcf ] ]
    tbi          // channel: [ val(meta), [ tbi ] ]
    fasta        // channel: /path/to/genome.fasta
    gff          // channel: /path/to/genome.gff
    nextclade_db // channel: /path/to/nextclade_db/

    main:

    ch_versions = Channel.empty()


    //
    // Filter variants by allele frequency
    //
    BCFTOOLS_FILTER(
        vcf
    )
    ch_versions = ch_versions.mix(BCFTOOLS_FILTER.out.versions.first())

    //
    // Create BED file with consensus regions to mask
    //
    MAKE_BED_MASK (
        bam.join(BCFTOOLS_FILTER.out.vcf, by: [0]),
        fasta,
        params.save_mpileup
    )
    ch_versions = ch_versions.mix(MAKE_BED_MASK.out.versions.first())

    //
    // Merge intervals with BEDTools
    //
    BEDTOOLS_MERGE (
        MAKE_BED_MASK.out.bed
    )
    ch_versions = ch_versions.mix(BEDTOOLS_MERGE.out.versions.first())

    //
    // Mask regions in consensus with BEDTools
    //
    BEDTOOLS_MASKFASTA (
        BEDTOOLS_MERGE.out.bed,
        fasta
    )
    ch_versions = ch_versions.mix(BEDTOOLS_MASKFASTA.out.versions.first())

    //
    // Call consensus sequence with BCFTools
    //
    BCFTOOLS_CONSENSUS (
        vcf.join(tbi, by: [0]).join(BEDTOOLS_MASKFASTA.out.fasta, by: [0])
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions.first())

    //
    // Consensus sequence QC
    //
    CONSENSUS_QC (
        BCFTOOLS_CONSENSUS.out.fasta,
        fasta,
        gff,
        nextclade_db
    )
    ch_versions = ch_versions.mix(CONSENSUS_QC.out.versions.first())

    emit:
    consensus        = BCFTOOLS_CONSENSUS.out.fasta      // channel: [ val(meta), [ fasta ] ]

    quast_results    = CONSENSUS_QC.out.quast_results    // channel: [ val(meta), [ results ] ]
    quast_tsv        = CONSENSUS_QC.out.quast_tsv        // channel: [ val(meta), [ tsv ] ]

    pangolin_report  = CONSENSUS_QC.out.pangolin_report  // channel: [ val(meta), [ csv ] ]

    nextclade_report = CONSENSUS_QC.out.nextclade_report // channel: [ val(meta), [ csv ] ]

    bases_tsv        = CONSENSUS_QC.out.bases_tsv        // channel: [ val(meta), [ tsv ] ]
    bases_pdf        = CONSENSUS_QC.out.bases_pdf        // channel: [ val(meta), [ pdf ] ]

    versions         = ch_versions                       // channel: [ versions.yml ]
}
