//
// Run snpEff, bgzip, tabix, stats and SnpSift commands
//
include { MINIMAP2_INDEX                               } from '../../../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN                               } from '../../../modules/nf-core/minimap2/align/main'
include { ARTIC_ALIGNTRIM                              } from '../../../modules/nf-core/artic/aligntrim/main'
include { CLAIR3                                       } from '../../../modules/nf-core/clair3/main'
include { BCFTOOLS_FILTER                              } from '../../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_NORM                                } from '../../../modules/nf-core/bcftools/norm/main'
include { BCFTOOLS_INDEX                               } from '../../../modules/nf-core/bcftools/index/main'
include { BAM_SORT_STATS_SAMTOOLS                      } from '../../nf-core/bam_sort_stats_samtools/main'
include { BCFTOOLS_FILTER as BCFTOOLS_CONSENSUS_FILTER } from '../../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_FILTER      } from '../../../modules/nf-core/bcftools/index/main'
include { BEDTOOLS_MERGE                               } from '../../../modules/nf-core/bedtools/merge/main'
include { BEDTOOLS_MASKFASTA                           } from '../../../modules/nf-core/bedtools/maskfasta/main'
include { BCFTOOLS_CONSENSUS                           } from '../../../modules/nf-core/bcftools/consensus/main'
include { MAKE_BED_MASK                                } from '../../../modules/local/make_bed_mask'
include { RENAME_FASTA_HEADER                          } from '../../../modules/local/rename_fasta_header'

workflow MINIMAP2_MAPPING {
    take:
    reads    // channel: [ val(meta), [ fastq ] ]
    fasta    // path   : genome.fasta
    fai      // path   : fai
    bed_file // path   : bed

    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    def ch_fasta_fai = fasta
        .combine(fai)
        .map { fasta_file, fai_file -> tuple([:], fasta_file, fai_file) }
        .collect(flat: false)
        .map { fasta_fai -> fasta_fai[0] }

    MINIMAP2_INDEX(
        fasta.map { fa -> tuple([:], fa) }
    )

    MINIMAP2_ALIGN(
        reads,
        MINIMAP2_INDEX.out.index.collect(flat: false).map { indexes -> indexes[0] },
        true,
        'bai',
        false,
        false)

    ch_minimap_bam = channel.empty()

    if (params.protocol == 'amplicon'){

        ch_input_bam_aligntrim = MINIMAP2_ALIGN.out.bam
            .join(MINIMAP2_ALIGN.out.index, by: [0])
            .combine(bed_file)
            .map { bam_bai_tuple, bed ->
                def meta = bam_bai_tuple[0]
                def bam  = bam_bai_tuple[1]

                tuple(meta, bam, bed, null)
            }

        ARTIC_ALIGNTRIM (
            ch_input_bam_aligntrim,
            false
        )

        ch_minimap_bam   = ARTIC_ALIGNTRIM.out.primertrimmed_bam
        ch_versions      = ch_versions.mix(ARTIC_ALIGNTRIM.out.versions)

        ch_multiqc_files = ch_multiqc_files.mix(ARTIC_ALIGNTRIM.out.align_trim_report)
        ch_multiqc_files = ch_multiqc_files.mix(ARTIC_ALIGNTRIM.out.amp_depth_report)

    } else {
        ch_minimap_bam = MINIMAP2_ALIGN.out.bam
    }

    BAM_SORT_STATS_SAMTOOLS (
        ch_minimap_bam,
        ch_fasta_fai
    )

    ch_multiqc_files = ch_multiqc_files.mix(BAM_SORT_STATS_SAMTOOLS.out.flagstat)

    if (params.clair3_model_dir){
        ch_input_bam_clair3 = BAM_SORT_STATS_SAMTOOLS.out.bam.join(BAM_SORT_STATS_SAMTOOLS.out.index, by: [0]).map { meta, bam, bai ->
            tuple(
                meta,
                bam,
                bai,
                null,
                file("${params.clair3_model_dir}/${params.clair3_model}", checkIfExists: true),
                'ont'
            )
        }
    } else {
        ch_input_bam_clair3 = BAM_SORT_STATS_SAMTOOLS.out.bam.join(BAM_SORT_STATS_SAMTOOLS.out.index, by: [0]).map { meta, bam, bai ->
            tuple(
                meta,
                bam,
                bai,
                params.clair3_model,
                [],
                'ont'
            )
        }
    }

    // Run CLAIR3
    CLAIR3(
        ch_input_bam_clair3,
        fasta.map { fa -> tuple([:], fa) },
        fai.map   { idx -> tuple([:], idx) }
    )

    ch_versions = ch_versions.mix(CLAIR3.out.versions.first())

    //
    // Filter variants by allele frequency, zip and index
    //
    BCFTOOLS_FILTER (
        CLAIR3.out.vcf.join(CLAIR3.out.tbi, by: [0])
    )

    BCFTOOLS_INDEX (
        BCFTOOLS_FILTER.out.vcf
    )

    //
    // Split multi-allelic positions and normalize
    //
    BCFTOOLS_NORM (
        BCFTOOLS_FILTER.out.vcf.join(BCFTOOLS_INDEX.out.tbi, by: [0]),
        fasta.map { fa -> tuple([:], fa) },
    )

    //
    // Filter variants by allele frequency, zip and index
    //

    BCFTOOLS_CONSENSUS_FILTER (
        BCFTOOLS_NORM.out.vcf.join(BCFTOOLS_NORM.out.tbi, by: [0])
    )

    BCFTOOLS_INDEX_FILTER (
        BCFTOOLS_CONSENSUS_FILTER.out.vcf
    )

    //
    // Create BED file with consensus regions to mask
    //
    MAKE_BED_MASK (
        BAM_SORT_STATS_SAMTOOLS.out.bam.join(BCFTOOLS_CONSENSUS_FILTER.out.vcf, by: [0]),
        fasta,
        params.save_mpileup
    )

    //
    // Merge intervals with BEDTools
    //
    BEDTOOLS_MERGE (
        MAKE_BED_MASK.out.bed
    )

    //
    // Mask regions in consensus with BEDTools
    //
    BEDTOOLS_MASKFASTA (
        BEDTOOLS_MERGE.out.bed,
        fasta
    )

    //
    // Call consensus sequence with BCFTools
    //
    BCFTOOLS_CONSENSUS (
        BCFTOOLS_CONSENSUS_FILTER.out.vcf
            .join(BCFTOOLS_INDEX_FILTER.out.tbi, by: [0])
            .join(BEDTOOLS_MASKFASTA.out.fasta, by: [0])
            .map { meta, vcf, tbi, mask_fasta -> tuple(meta, vcf, tbi, mask_fasta, []) }
    )

    //
    // Rename consensus header adding sample name
    //
    RENAME_FASTA_HEADER (
        BCFTOOLS_CONSENSUS.out.fasta
    )

    emit:
    bam          = BAM_SORT_STATS_SAMTOOLS.out.bam
    bai          = BAM_SORT_STATS_SAMTOOLS.out.index

    vcf           = BCFTOOLS_NORM.out.vcf
    tbi           = BCFTOOLS_NORM.out.tbi

    consensus     = RENAME_FASTA_HEADER.out.fasta

    multiqc_files = ch_multiqc_files

    versions      = ch_versions    // channel: [ versions.yml ]
}
