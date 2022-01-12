//
// Run various tools to generate a masked genome consensus sequence
//

include { BEDTOOLS_GENOMECOV } from '../../modules/nf-core/modules/bedtools/genomecov/main'
include { BEDTOOLS_MERGE     } from '../../modules/nf-core/modules/bedtools/merge/main'
include { BEDTOOLS_MASKFASTA } from '../../modules/nf-core/modules/bedtools/maskfasta/main'
include { BCFTOOLS_CONSENSUS } from '../../modules/nf-core/modules/bcftools/consensus/main'
include { MAKE_BED_MASK      } from '../../modules/local/make_bed_mask'
include { PLOT_BASE_DENSITY  } from '../../modules/local/plot_base_density'

workflow MAKE_CONSENSUS {
    take:
    bam_vcf // channel: [ val(meta), [ bam ], [ vcf ], [ tbi ] ]
    fasta

    main:

    ch_versions = Channel.empty()

    BEDTOOLS_GENOMECOV (
        bam_vcf.map { meta, bam, vcf, tbi -> [ meta, bam, 1 ] },
        [],
        'bed',
    )
    ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions.first())

    BEDTOOLS_MERGE (
        BEDTOOLS_GENOMECOV.out.genomecov
    )
    ch_versions = ch_versions.mix(BEDTOOLS_MERGE.out.versions.first())

    MAKE_BED_MASK (
        bam_vcf.map { meta, bam, vcf, tbi -> [ meta, vcf ] }.join( BEDTOOLS_MERGE.out.bed, by: [0] ),
        fasta
    )
    ch_versions = ch_versions.mix(MAKE_BED_MASK.out.versions.first())

    BEDTOOLS_MASKFASTA (
        MAKE_BED_MASK.out.bed,
        MAKE_BED_MASK.out.fasta.map{it[1]}
    )
    ch_versions = ch_versions.mix(BEDTOOLS_MASKFASTA.out.versions.first())

    BCFTOOLS_CONSENSUS (
        bam_vcf.map { meta, bam, vcf, tbi -> [ meta, vcf, tbi ] }.join( BEDTOOLS_MASKFASTA.out.fasta, by: [0] )
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions.first())

    PLOT_BASE_DENSITY (
        BCFTOOLS_CONSENSUS.out.fasta
    )
    ch_versions = ch_versions.mix(PLOT_BASE_DENSITY.out.versions.first())

    emit:
    fasta    = BCFTOOLS_CONSENSUS.out.fasta // channel: [ val(meta), [ fasta ] ]
    tsv      = PLOT_BASE_DENSITY.out.tsv    // channel: [ val(meta), [ tsv ] ]
    pdf      = PLOT_BASE_DENSITY.out.pdf    // channel: [ val(meta), [ pdf ] ]

    versions = ch_versions                  // channel: [ versions.yml ]
}
