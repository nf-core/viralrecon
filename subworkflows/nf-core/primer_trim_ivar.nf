//
// iVar trim, sort, index BAM file and run samtools stats, flagstat and idxstats
//

include { IVAR_TRIM               } from '../../modules/nf-core/ivar/trim/main'
include { BAM_SORT_STATS_SAMTOOLS } from './bam_sort_stats_samtools/main'

workflow PRIMER_TRIM_IVAR {
    take:
    bam   // channel: [ val(meta), [ bam ], [bai] ]
    fasta // path   : fasta
    bed   // path   : bed

    main:

    ch_versions = Channel.empty()

    //
    // iVar trim primers
    //
    IVAR_TRIM (
        bam,
        bed
    )
    ch_versions = ch_versions.mix(IVAR_TRIM.out.versions.first())

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS ( IVAR_TRIM.out.bam, fasta )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:
    bam_orig = IVAR_TRIM.out.bam              // channel: [ val(meta), bam   ]
    log_out  = IVAR_TRIM.out.log              // channel: [ val(meta), log   ]

    bam      = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai      = BAM_SORT_STATS_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions = ch_versions                    // channel: [ versions.yml ]
}
