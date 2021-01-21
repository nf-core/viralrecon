/*
 * Run BCFTools tabix and stats commands
 */

params.tabix_options = [:]
params.stats_options = [:]

include { BCFTOOLS_TABIX } from '../../modules/local/bcftools_tabix' addParams( options: params.tabix_options )
include { BCFTOOLS_STATS } from '../../modules/local/bcftools_stats' addParams( options: params.stats_options )

workflow VCF_TABIX_STATS {
    take:
    vcf // channel: [ val(meta), [ vcf ] ]
    
    main:
    BCFTOOLS_TABIX ( vcf )
    BCFTOOLS_STATS ( vcf )

    emit:
    tbi     = BCFTOOLS_TABIX.out.tbi     // channel: [ val(meta), [ tbi ] ]
    stats   = BCFTOOLS_STATS.out.stats   // channel: [ val(meta), [ txt ] ]
    version = BCFTOOLS_STATS.out.version //    path: *.version.txt
}
