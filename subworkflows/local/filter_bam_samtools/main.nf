//
// Filter co-ordinate sorted BAM, index and run samtools stats, flagstat and idxstats
//

include { SAMTOOLS_VIEW      } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../../nf-core/bam_stats_samtools/main'

workflow FILTER_BAM_SAMTOOLS {
    take:
    bam_bai // channel: [ val(meta), [ bam ], [ bai ] ]
    fasta   // channel: [ val(meta), path(fasta), path(fai) ]

    main:

    //
    // Filter BAM using Samtools view
    //
    SAMTOOLS_VIEW (
        bam_bai,
        fasta,
        [ [:], [] ],
        [ [:], [] ],
        []
    )

    //
    // Index BAM file and run samtools stats, flagstat and idxstats
    //
    SAMTOOLS_INDEX (
        SAMTOOLS_VIEW.out.bam
    )

    BAM_STATS_SAMTOOLS (
        SAMTOOLS_VIEW.out.bam.join(SAMTOOLS_INDEX.out.index, by: [0]),
        fasta
    )

    emit:
    bam      = SAMTOOLS_VIEW.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.index        // channel: [ val(meta), [ bai ] ]
    stats    = BAM_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
}
