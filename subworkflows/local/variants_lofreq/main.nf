//
// Variant calling with LoFreq, downstream processing and QC
//
include { LOFREQ_INDELQUAL                 } from '../../../modules/nf-core/lofreq/indelqual/main'
include { LOFREQ_CALLPARALLEL              } from '../../../modules/nf-core/lofreq/callparallel/main'
include { LOFREQ_FILTER                    } from '../../../modules/nf-core/lofreq/filter/main'
include { SAMTOOLS_INDEX                   } from '../../../modules/nf-core/samtools/index/main'
include { BCFTOOLS_STATS                   } from '../../../modules/nf-core/bcftools/stats/main'
include { getNumVariantsFromBCFToolsStats  } from '../../../subworkflows/local/utils_nfcore_viralrecon_pipeline'
include { VCF_BGZIP_TABIX_STATS            } from '../vcf_bgzip_tabix_stats'
include { VARIANTS_QC                      } from '../variants_qc'

workflow VARIANTS_LOFREQ {
    take:
    bam                 // channel: [ val(meta), [ bam ] ]
    fasta               // channel: /path/to/genome.fasta
    fai                 // channel: /path/to/genome.fai
    gff                 // channel: /path/to/genome.gff
    snpeff_db           // channel: /path/to/snpeff_db/
    snpeff_config       // channel: /path/to/snpeff.config

    main:

    // Add indel qualities to BAM files
    LOFREQ_INDELQUAL (
        bam,
        fasta.map { fasta_file -> [ [:], fasta_file ] }
    )

    // Index the bam files with indel qualities
    SAMTOOLS_INDEX (
        LOFREQ_INDELQUAL.out.bam
    )

    // Call variants in parallel
    ch_lofreq_bam = LOFREQ_INDELQUAL.out.bam
        .join(SAMTOOLS_INDEX.out.index)
        .map { meta, bam_file, bai -> [meta, bam_file, bai, []] }

    LOFREQ_CALLPARALLEL (
        ch_lofreq_bam,
        fasta.map { fasta_file -> [ [:], fasta_file ] },
        fai.map { fai_file -> [ [:], fai_file ] }
    )

    //Filter out samples with 0 variants
    ch_raw_vcf_tbi =  LOFREQ_CALLPARALLEL.out.vcf
        .join(LOFREQ_CALLPARALLEL.out.tbi)
        .map { meta, vcf, tbi -> [meta, vcf, tbi] }

    BCFTOOLS_STATS (
        ch_raw_vcf_tbi,
        [[],[]],
        [[],[]],
        [[],[]],
        [[],[]],
        fasta.map { fasta_file -> [ [:], fasta_file ] }
    )

    ch_vcf_tbi_stats = ch_raw_vcf_tbi
        .join(BCFTOOLS_STATS.out.stats)
        .filter { _meta, _vcf, _tbi, stats ->
            getNumVariantsFromBCFToolsStats(stats) > 0
        }

    // Filter variants by minimum allele frequency, which is not included in the default filters applied by lofreq call.
    LOFREQ_FILTER (
        ch_vcf_tbi_stats.map { meta, vcf, _tbi, _stats -> [ meta, vcf]}
    )

    // Index the filtered VCF files
    VCF_BGZIP_TABIX_STATS (
        LOFREQ_FILTER.out.vcf,
        [ [:], [] ],
        [ [:], [] ],
        [ [:], [] ]
    )

    // Keep matching VCF, index and stats only for samples with records.
    ch_filtered_vcf_tbi_stats = VCF_BGZIP_TABIX_STATS.out.vcf
        .join(VCF_BGZIP_TABIX_STATS.out.tbi, by: [0])
        .join(VCF_BGZIP_TABIX_STATS.out.stats, by: [0])
        .filter { _meta, _vcf, _tbi, stats ->
            getNumVariantsFromBCFToolsStats(stats) > 0
        }

    ch_lofreq_vcf = ch_filtered_vcf_tbi_stats
        .map { meta, vcf, _tbi, _stats -> [meta, vcf] }

    ch_lofreq_tbi = ch_filtered_vcf_tbi_stats
        .map { meta, _vcf, tbi, _stats -> [meta, tbi] }

    ch_lofreq_stats = ch_filtered_vcf_tbi_stats
        .map { meta, _vcf, _tbi, stats -> [meta, stats] }

    // Run downstream tools for variants QC
    VARIANTS_QC (
        ch_lofreq_vcf,
        fasta,
        gff,
        snpeff_db,
        snpeff_config
    )

    emit:
    vcf_orig        = ch_vcf_tbi_stats.map{meta, vcf, _tbi, _stats -> [ meta, vcf]}                          // channel: [ val(meta), [ vcf ] ]
    tbi_orig        = ch_vcf_tbi_stats.map{meta, _vcf, tbi, _stats -> [ meta, tbi]}                          // channel: [ val(meta), [ tbi ] ]
    stats_orig      = ch_vcf_tbi_stats.map{meta, _vcf, _tbi, stats -> [ meta, stats]}                          // channel: [ val(meta), [ txt ] ]

    vcf             = ch_lofreq_vcf                   // channel: [ val(meta), [ vcf.gz ] ]
    tbi             = ch_lofreq_tbi                   // channel: [ val(meta), [ tbi ] ]
    stats           = ch_lofreq_stats                 // channel: [ val(meta), [ txt ] ]

    snpeff_vcf      = VARIANTS_QC.out.snpeff_vcf       // channel: [ val(meta), [ vcf.gz ] ]
    snpeff_tbi      = VARIANTS_QC.out.snpeff_tbi       // channel: [ val(meta), [ tbi ] ]
    snpeff_stats    = VARIANTS_QC.out.snpeff_stats     // channel: [ val(meta), [ txt ] ]
    snpeff_csv      = VARIANTS_QC.out.snpeff_csv       // channel: [ val(meta), [ csv ] ]
    snpeff_txt      = VARIANTS_QC.out.snpeff_txt       // channel: [ val(meta), [ txt ] ]
    snpeff_html     = VARIANTS_QC.out.snpeff_html      // channel: [ val(meta), [ html ] ]
    snpsift_txt     = VARIANTS_QC.out.snpsift_txt      // channel: [ val(meta), [ txt ] ]
}
