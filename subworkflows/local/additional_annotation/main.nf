//
// Run snpEff, bgzip, tabix, stats and SnpSift commands
//

include { SNPEFF_BUILD                                                    } from '../../../modules/local/snpeff/build'
include { SNPEFF_ANN                                                      } from '../../../modules/local/snpeff/ann'
include { SNPSIFT_EXTRACTFIELDS                                           } from '../../../modules/local/snpsift/extractfields'
include { VCF_BGZIP_TABIX_STATS                                           } from '../vcf_bgzip_tabix_stats'
include { BCFTOOLS_QUERY                                                  } from '../../../modules/nf-core/bcftools/query/main'
include { MAKE_VARIANTS_LONG_TABLE as MAKE_VARIANTS_LONG_TABLE_ADDITIONAL } from '../../../modules/local/make_variants_long_table'


workflow ADDITIONAL_ANNOTATION {
    take:
    vcf      // channel: [ val(meta), [ vcf ] ]
    tbi      // channel: [ val(meta), [ tbi ] ]
    fasta    // path   : genome.fasta
    annot    // path   : additional_annotation
    pangolin // channel: [ val(meta), [ csv ] ]

    main:

    //
    // Make snpEff database
    //
    ch_snpeff_db     = channel.empty()
    ch_snpeff_config = channel.empty()

    SNPEFF_BUILD (
        fasta,
        annot
    )
    ch_snpeff_db     = SNPEFF_BUILD.out.db
    ch_snpeff_config = SNPEFF_BUILD.out.config

    SNPEFF_ANN (
        vcf,
        ch_snpeff_db,
        ch_snpeff_config,
        fasta
    )

    VCF_BGZIP_TABIX_STATS (
        SNPEFF_ANN.out.vcf,
        [ [:], [] ],
        [ [:], [] ],
        [ [:], [] ]
    )

    SNPSIFT_EXTRACTFIELDS (
        VCF_BGZIP_TABIX_STATS.out.vcf
    )

    BCFTOOLS_QUERY (
        vcf.join(tbi, by: [0]),
        [],
        [],
        []
    )

    MAKE_VARIANTS_LONG_TABLE_ADDITIONAL (
        BCFTOOLS_QUERY.out.output.collect{output -> output[1]},
        SNPSIFT_EXTRACTFIELDS.out.txt.collect{txt -> txt[1]}.ifEmpty([]),
        pangolin.collect{pgl -> pgl[1]}.ifEmpty([])
    )

    emit:
    long_table  = MAKE_VARIANTS_LONG_TABLE_ADDITIONAL.out.csv // channel: [ val(meta), [ csv ] ]
}
