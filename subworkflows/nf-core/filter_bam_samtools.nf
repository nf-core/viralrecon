//
// Filter co-ordinate sorted BAM, index and run samtools stats, flagstat and idxstats
//

include { SAMTOOLS_VIEW      } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_INDEX     } from '../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from './bam_stats_samtools/main'

workflow FILTER_BAM_SAMTOOLS {
    take:
    ch_bam_bai // channel: [ val(meta), [ bam ], [ bai ] ]
    ch_fasta   // path   : fasta

    main:

    ch_versions = Channel.empty()

    //
    // Filter BAM using Samtools view
    //
    SAMTOOLS_VIEW ( ch_bam_bai, ch_fasta )
    ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions.first())

    //
    // Index BAM file and run samtools stats, flagstat and idxstats
    //
    SAMTOOLS_INDEX ( SAMTOOLS_VIEW.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    SAMTOOLS_VIEW.out.bam
        .join(SAMTOOLS_INDEX.out.bai, by: [0], remainder: true)
        .set { ch_bam_bai }

    BAM_STATS_SAMTOOLS ( ch_bam_bai, ch_fasta )

    ch_versions = ch_versions.mix(BAM_STATS_SAMTOOLS.out.versions)

    emit:
    bam      = SAMTOOLS_VIEW.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    stats    = BAM_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}
