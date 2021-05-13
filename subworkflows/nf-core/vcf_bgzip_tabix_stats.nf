//
// Run BCFTools bgzip, tabix and stats commands
//

params.bgzip_options = [:]
params.tabix_options = [:]
params.stats_options = [:]

include { TABIX_BGZIP     } from '../../modules/nf-core/software/tabix/bgzip/main' addParams( options: params.bgzip_options                                            )
include { VCF_TABIX_STATS } from './vcf_tabix_stats'                               addParams( tabix_options: params.tabix_options, stats_options: params.stats_options )

workflow VCF_BGZIP_TABIX_STATS {
    take:
    vcf // channel: [ val(meta), [ vcf ] ]

    main:
    TABIX_BGZIP  ( vcf )
    VCF_TABIX_STATS ( TABIX_BGZIP.out.gz )

    emit:
    vcf              = TABIX_BGZIP.out.gz                   // channel: [ val(meta), [ vcf.gz ] ]
    tabix_version    = TABIX_BGZIP.out.version              //    path: *.version.txt

    tbi              = VCF_TABIX_STATS.out.tbi              // channel: [ val(meta), [ tbi ] ]
    stats            = VCF_TABIX_STATS.out.stats            // channel: [ val(meta), [ txt ] ]
    bcftools_version = VCF_TABIX_STATS.out.bcftools_version //    path: *.version.txt
}
