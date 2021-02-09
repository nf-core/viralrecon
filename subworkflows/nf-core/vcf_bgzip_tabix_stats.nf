/*
 * Run BCFTools bgzip, tabix and stats commands
 */

params.bgzip_options = [:]
params.tabix_options = [:]
params.stats_options = [:]

include { BCFTOOLS_BGZIP  } from '../../modules/nf-core/software/bcftools/bgzip/main' addParams( options: params.bgzip_options )
include { VCF_TABIX_STATS } from './vcf_tabix_stats'                                  addParams( tabix_options: params.tabix_options, stats_options: params.stats_options )

workflow VCF_BGZIP_TABIX_STATS {
    take:
    vcf // channel: [ val(meta), [ vcf ] ]
    
    main:
    BCFTOOLS_BGZIP  ( vcf )
    VCF_TABIX_STATS ( BCFTOOLS_BGZIP.out.vcf )

    emit:
    vcf     = BCFTOOLS_BGZIP.out.vcf     // channel: [ val(meta), [ vcf.gz ] ]
    tbi     = VCF_TABIX_STATS.out.tbi    // channel: [ val(meta), [ tbi ] ]
    stats   = VCF_TABIX_STATS.out.stats  // channel: [ val(meta), [ txt ] ]
    version = BCFTOOLS_BGZIP.out.version //    path: *.version.txt
}
