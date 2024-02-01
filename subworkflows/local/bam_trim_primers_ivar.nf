//
// iVar trim, sort, index BAM file and run samtools stats, flagstat and idxstats
//

include { IVAR_TRIM               } from '../../modules/nf-core/ivar/trim/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../nf-core/bam_sort_stats_samtools/main'
include { SAMTOOLS_COLLATE        } from '../../modules/nf-core/samtools/collate/main'
include { SAMTOOLS_FIXMATE        } from '../../modules/nf-core/samtools/fixmate/main'

workflow BAM_TRIM_PRIMERS_IVAR {
    take:
    bam   // channel: [ val(meta), [ bam ], [bai] ]
    bed   // path   : bed
    fasta // channel: reference.fasta

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
    // samtools fixmate fills in mate coordinates and insert size fields
    //
    SAMTOOLS_COLLATE ( IVAR_TRIM.out.bam, fasta )
    ch_versions = ch_versions.mix(SAMTOOLS_COLLATE.out.versions.first())

    SAMTOOLS_FIXMATE ( SAMTOOLS_COLLATE.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_FIXMATE.out.versions.first())

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    BAM_SORT_STATS_SAMTOOLS (
        SAMTOOLS_FIXMATE.out.bam,
        fasta
    )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    emit:
    bam_orig = SAMTOOLS_FIXMATE.out.bam             // channel: [ val(meta), bam   ]
    log_out  = IVAR_TRIM.out.log                    // channel: [ val(meta), log   ]

    bam      = BAM_SORT_STATS_SAMTOOLS.out.bam      // channel: [ val(meta), [ bam ] ]
    bai      = BAM_SORT_STATS_SAMTOOLS.out.bai      // channel: [ val(meta), [ bai ] ]
    csi      = BAM_SORT_STATS_SAMTOOLS.out.csi      // channel: [ val(meta), [ csi ] ]
    stats    = BAM_SORT_STATS_SAMTOOLS.out.stats    // channel: [ val(meta), [ stats ] ]
    flagstat = BAM_SORT_STATS_SAMTOOLS.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = BAM_SORT_STATS_SAMTOOLS.out.idxstats // channel: [ val(meta), [ idxstats ] ]

    versions = ch_versions                          // channel: [ versions.yml ]
}
