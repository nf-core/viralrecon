//
// Variant calling with BCFTools, downstream processing and QC
//

include { BCFTOOLS_MPILEUP                 } from '../../../modules/nf-core/bcftools/mpileup/main'
include { BCFTOOLS_NORM                    } from '../../../modules/nf-core/bcftools/norm/main'
include { VCF_TABIX_STATS                  } from '../vcf_tabix_stats'
include { VARIANTS_QC                      } from '../variants_qc'
include { getNumVariantsFromBCFToolsStats  } from '../../../subworkflows/local/utils_nfcore_viralrecon_pipeline'

workflow VARIANTS_BCFTOOLS {
    take:
    bam           // channel: [ val(meta), [ bam ] ]
    fasta_fai     // channel: [ val(meta), fasta, fai ]
    gff           // channel: /path/to/genome.gff
    snpeff_db     // channel: /path/to/snpeff_db/
    snpeff_config // channel: /path/to/snpeff.config

    main:

    //
    // Call variants
    //
    ch_fasta = fasta_fai.map { _meta, fasta_file, _fai_file -> fasta_file }

    BCFTOOLS_MPILEUP (
        bam.map{ meta, bam_file -> [ meta, bam_file, [], [] ] },
        fasta_fai,
        params.save_mpileup
    )

    // Filter out samples with 0 variants
    BCFTOOLS_MPILEUP
        .out
        .vcf
        .join(BCFTOOLS_MPILEUP.out.tbi)
        .join(BCFTOOLS_MPILEUP.out.stats)
        .filter { _meta, _vcf, _tbi, stats -> getNumVariantsFromBCFToolsStats(stats) > 0 }
        .set { ch_vcf_tbi_stats }

    ch_vcf_tbi_stats
        .map { meta, vcf, _tbi, _stats -> [ meta, vcf ] }
        .set { ch_vcf }

    ch_vcf_tbi_stats
        .map { meta, _vcf, tbi, _stats -> [ meta, tbi ] }
        .set { ch_tbi }

    ch_vcf_tbi_stats
        .map { meta, _vcf, _tbi, stats -> [ meta, stats ] }
        .set { ch_stats }

    //
    // Split multi-allelic positions
    //
    BCFTOOLS_NORM (
        ch_vcf.join(ch_tbi, by: [0]),
        ch_fasta.map { fasta_file -> [ [:], fasta_file ] }
    )

    VCF_TABIX_STATS (
        BCFTOOLS_NORM.out.vcf,
        [ [:], [] ],
        [ [:], [] ],
        [ [:], [] ]
    )

    //
    // Run downstream tools for variants QC
    //
    VARIANTS_QC (
        BCFTOOLS_NORM.out.vcf,
        ch_fasta,
        gff,
        snpeff_db,
        snpeff_config
    )

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
}
