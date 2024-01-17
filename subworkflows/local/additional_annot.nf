//
// Run snpEff, bgzip, tabix, stats and SnpSift commands
//

include { SNPEFF_BUILD                                                    } from '../../modules/local/snpeff_build'
include { SNPEFF_ANN                                                      } from '../../modules/local/snpeff_ann'
include { SNPSIFT_EXTRACTFIELDS                                           } from '../../modules/local/snpsift_extractfields'
include { VCF_BGZIP_TABIX_STATS                                           } from './vcf_bgzip_tabix_stats'
include { BCFTOOLS_QUERY                                                  } from '../../modules/nf-core/bcftools/query/main'
include { MAKE_VARIANTS_LONG_TABLE as MAKE_VARIANTS_LONG_TABLE_ADDITIONAL } from '../../modules/local/make_variants_long_table'


workflow ADDITIONAL_ANNOT {
    take:
    vcf      // channel: [ val(meta), [ vcf ] ]
    tbi      // channel: [ val(meta), [ tbi ] ]
    fasta    // path   : genome.fasta
    annot    // path   : additional_annot
    pangolin // channel: [ val(meta), [ csv ] ]

    main:

    ch_versions = Channel.empty()

    //
    // Make snpEff database
    //
    ch_snpeff_db     = Channel.empty()
    ch_snpeff_config = Channel.empty()

    SNPEFF_BUILD (
        fasta,
        annot
    )
    ch_snpeff_db     = SNPEFF_BUILD.out.db
    ch_snpeff_config = SNPEFF_BUILD.out.config
    ch_versions      = ch_versions.mix(SNPEFF_BUILD.out.versions)

    SNPEFF_ANN (
        vcf,
        ch_snpeff_db,
        ch_snpeff_config,
        fasta
    )
    ch_versions = ch_versions.mix(SNPEFF_ANN.out.versions.first())

    VCF_BGZIP_TABIX_STATS (
        SNPEFF_ANN.out.vcf,
        [ [:], [] ],
        [ [:], [] ],
        [ [:], [] ]
    )
    ch_versions = ch_versions.mix(VCF_BGZIP_TABIX_STATS.out.versions)

    SNPSIFT_EXTRACTFIELDS (
        VCF_BGZIP_TABIX_STATS.out.vcf
    )
    ch_versions = ch_versions.mix(SNPSIFT_EXTRACTFIELDS.out.versions.first())

    BCFTOOLS_QUERY (
        vcf.join(tbi, by: [0]),
        [],
        [],
        []
    )
    ch_versions = ch_versions.mix(BCFTOOLS_QUERY.out.versions.first())

    MAKE_VARIANTS_LONG_TABLE_ADDITIONAL (
        BCFTOOLS_QUERY.out.output.collect{it[1]},
        SNPSIFT_EXTRACTFIELDS.out.txt.collect{it[1]}.ifEmpty([]),
        pangolin.collect{it[1]}.ifEmpty([])
    )
    ch_versions = ch_versions.mix(MAKE_VARIANTS_LONG_TABLE_ADDITIONAL.out.versions)

    emit:
    long_table  = MAKE_VARIANTS_LONG_TABLE_ADDITIONAL.out.csv // channel: [ val(meta), [ csv ] ]

    versions    = ch_versions    // channel: [ versions.yml ]
}
