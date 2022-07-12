//
// Variant calling with BCFTools, downstream processing and QC
//

include { BCFTOOLS_MPILEUP } from '../../modules/nf-core/modules/bcftools/mpileup/main'
include { BCFTOOLS_NORM    } from '../../modules/nf-core/modules/bcftools/norm/main'
include { VCF_TABIX_STATS  } from '../nf-core/vcf_tabix_stats'
include { VARIANTS_QC      } from './variants_qc'

workflow VARIANTS_BCFTOOLS {
    take:
    bam           // channel: [ val(meta), [ bam ] ]
    fasta         // channel: /path/to/genome.fasta
    sizes         // channel: /path/to/genome.sizes
    gff           // channel: /path/to/genome.gff
    bed           // channel: /path/to/primers.bed
    snpeff_db     // channel: /path/to/snpeff_db/
    snpeff_config // channel: /path/to/snpeff.config

    main:

    ch_versions = Channel.empty()

    //
    // Call variants
    //
    BCFTOOLS_MPILEUP (
        bam,
        fasta,
        params.save_mpileup
    )
    ch_versions = ch_versions.mix(BCFTOOLS_MPILEUP.out.versions.first())

    // Filter out samples with 0 variants
    BCFTOOLS_MPILEUP
        .out
        .vcf
        .join(BCFTOOLS_MPILEUP.out.tbi)
        .join(BCFTOOLS_MPILEUP.out.stats)
        .filter { meta, vcf, tbi, stats -> WorkflowCommons.getNumVariantsFromBCFToolsStats(stats) > 0 }
        .set { ch_vcf_tbi_stats }

    ch_vcf_tbi_stats
        .map { meta, vcf, tbi, stats -> [ meta, vcf ] }
        .set { ch_vcf }

    ch_vcf_tbi_stats
        .map { meta, vcf, tbi, stats -> [ meta, tbi ] }
        .set { ch_tbi }

    ch_vcf_tbi_stats
        .map { meta, vcf, tbi, stats -> [ meta, stats ] }
        .set { ch_stats }

    //
    // Split multi-allelic positions
    //
    BCFTOOLS_NORM (
        ch_vcf.join(ch_tbi, by: [0]),
        fasta
    )
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions.first())

    VCF_TABIX_STATS (
        BCFTOOLS_NORM.out.vcf
    )
    ch_versions = ch_versions.mix(VCF_TABIX_STATS.out.versions)

    //
    // Run downstream tools for variants QC
    //
    VARIANTS_QC (
        bam,
        BCFTOOLS_NORM.out.vcf,
        VCF_TABIX_STATS.out.stats,
        fasta,
        sizes,
        gff,
        bed,
        snpeff_db,
        snpeff_config
    )
    ch_versions = ch_versions.mix(VARIANTS_QC.out.versions)

    emit:
    vcf_orig        = ch_vcf                          // channel: [ val(meta), [ vcf ] ]
    tbi_orig        = ch_tbi                          // channel: [ val(meta), [ tbi ] ]
    stats_orig      = ch_stats                        // channel: [ val(meta), [ txt ] ]

    vcf             = BCFTOOLS_NORM.out.vcf           // channel: [ val(meta), [ vcf ] ]
    tbi             = VCF_TABIX_STATS.out.tbi         // channel: [ val(meta), [ tbi ] ]
    stats           = VCF_TABIX_STATS.out.stats       // channel: [ val(meta), [ txt ] ]

    snpeff_vcf      = VARIANTS_QC.out.snpeff_vcf      // channel: [ val(meta), [ vcf.gz ] ]
    snpeff_tbi      = VARIANTS_QC.out.snpeff_tbi      // channel: [ val(meta), [ tbi ] ]
    snpeff_stats    = VARIANTS_QC.out.snpeff_stats    // channel: [ val(meta), [ txt ] ]
    snpeff_csv      = VARIANTS_QC.out.snpeff_csv      // channel: [ val(meta), [ csv ] ]
    snpeff_txt      = VARIANTS_QC.out.snpeff_txt      // channel: [ val(meta), [ txt ] ]
    snpeff_html     = VARIANTS_QC.out.snpeff_html     // channel: [ val(meta), [ html ] ]
    snpsift_txt     = VARIANTS_QC.out.snpsift_txt     // channel: [ val(meta), [ txt ] ]

    asciigenome_pdf = VARIANTS_QC.out.asciigenome_pdf // channel: [ val(meta), [ pdf ] ]

    versions        = ch_versions                     // channel: [ versions.yml ]
}
