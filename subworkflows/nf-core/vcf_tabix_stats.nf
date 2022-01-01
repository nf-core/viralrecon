//
// Run BCFTools tabix and stats commands
//

params.tabix_options = [:]
params.stats_options = [:]

include { TABIX_TABIX    } from '../../modules/nf-core/modules/tabix/tabix/main'    addParams( options: params.tabix_options )
include { BCFTOOLS_STATS } from '../../modules/nf-core/modules/bcftools/stats/main' addParams( options: params.stats_options )

workflow VCF_TABIX_STATS {
    take:
    vcf // channel: [ val(meta), [ vcf ] ]

    main:
    TABIX_TABIX    ( vcf )
    BCFTOOLS_STATS ( vcf )

    emit:
    tbi              = TABIX_TABIX.out.tbi        // channel: [ val(meta), [ tbi ] ]
    tabix_version    = TABIX_TABIX.out.version    //    path: *.version.txt

    stats            = BCFTOOLS_STATS.out.stats   // channel: [ val(meta), [ txt ] ]
    bcftools_version = BCFTOOLS_STATS.out.version //    path: *.version.txt
}
