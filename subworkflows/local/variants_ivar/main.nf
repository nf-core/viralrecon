//
// Variant calling with IVar, downstream processing and QC
//

include { IVAR_VARIANTS         } from '../../../modules/nf-core/ivar/variants/main'
include { IVAR_VARIANTS_TO_VCF  } from '../../../modules/local/ivar_variants_to_vcf'
include { BCFTOOLS_SORT         } from '../../../modules/nf-core/bcftools/sort/main'
include { VCF_TABIX_STATS       } from '../vcf_tabix_stats'
include { VARIANTS_QC           } from '../variants_qc'
include { getNumLinesInFile     } from '../../../subworkflows/local/utils_nfcore_viralrecon_pipeline'


workflow VARIANTS_IVAR {
    take:
    bam                 // channel: [ val(meta), [ bam ] ]
    fasta               // channel: /path/to/genome.fasta
    fai                 // channel: /path/to/genome.fai
    gff                 // channel: /path/to/genome.gff
    snpeff_db           // channel: /path/to/snpeff_db/
    snpeff_config       // channel: /path/to/snpeff.config
    ivar_multiqc_header // channel: /path/to/multiqc_header for ivar variants

    main:

    ch_versions = channel.empty()

    //
    // Call variants
    //
    IVAR_VARIANTS (
        bam,
        fasta,
        fai,
        gff,
        params.save_mpileup
    )
    ch_versions = ch_versions.mix(IVAR_VARIANTS.out.versions)

    // Filter out samples with 0 variants
    IVAR_VARIANTS
        .out
        .tsv
        .filter { _meta, tsv -> getNumLinesInFile(tsv) > 1 }
        .set { ch_ivar_tsv }

    //
    // Convert original iVar output to VCF, zip and index
    //
    IVAR_VARIANTS_TO_VCF (
        ch_ivar_tsv,
        fasta,
        ivar_multiqc_header
    )

    BCFTOOLS_SORT (
        IVAR_VARIANTS_TO_VCF.out.vcf
    )

    VCF_TABIX_STATS (
        BCFTOOLS_SORT.out.vcf,
        [ [:], [] ],
        [ [:], [] ],
        [ [:], [] ]
    )

    //
    // Run downstream tools for variants QC
    //
    VARIANTS_QC (
        BCFTOOLS_SORT.out.vcf,
        fasta,
        gff,
        snpeff_db,
        snpeff_config
    )

    emit:
    tsv             = ch_ivar_tsv                     // channel: [ val(meta), [ tsv ] ]

    vcf_orig        = IVAR_VARIANTS_TO_VCF.out.vcf    // channel: [ val(meta), [ vcf ] ]
    log_out         = IVAR_VARIANTS_TO_VCF.out.log    // channel: [ val(meta), [ log ] ]
    multiqc_tsv     = IVAR_VARIANTS_TO_VCF.out.tsv    // channel: [ val(meta), [ tsv ] ]

    vcf             = BCFTOOLS_SORT.out.vcf           // channel: [ val(meta), [ vcf ] ]
    tbi             = VCF_TABIX_STATS.out.tbi         // channel: [ val(meta), [ tbi ] ]
    stats           = VCF_TABIX_STATS.out.stats       // channel: [ val(meta), [ txt ] ]

    snpeff_vcf      = VARIANTS_QC.out.snpeff_vcf      // channel: [ val(meta), [ vcf.gz ] ]
    snpeff_tbi      = VARIANTS_QC.out.snpeff_tbi      // channel: [ val(meta), [ tbi ] ]
    snpeff_stats    = VARIANTS_QC.out.snpeff_stats    // channel: [ val(meta), [ txt ] ]
    snpeff_csv      = VARIANTS_QC.out.snpeff_csv      // channel: [ val(meta), [ csv ] ]
    snpeff_txt      = VARIANTS_QC.out.snpeff_txt      // channel: [ val(meta), [ txt ] ]
    snpeff_html     = VARIANTS_QC.out.snpeff_html     // channel: [ val(meta), [ html ] ]
    snpsift_txt     = VARIANTS_QC.out.snpsift_txt     // channel: [ val(meta), [ txt ] ]
    versions        = ch_versions                     // channel: versions.yml
}
