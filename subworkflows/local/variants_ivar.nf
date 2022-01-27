//
// Variant calling with IVar, downstream processing and QC
//

include { IVAR_VARIANTS         } from '../../modules/nf-core/modules/ivar/variants/main'
include { IVAR_VARIANTS_TO_VCF  } from '../../modules/local/ivar_variants_to_vcf'
include { BCFTOOLS_SORT         } from '../../modules/nf-core/modules/bcftools/sort/main'


include { VCF_BGZIP_TABIX_STATS } from '../nf-core/vcf_bgzip_tabix_stats'
include { VARIANTS_QC           } from './variants_qc'

workflow VARIANTS_IVAR {
    take:
    bam                 // channel: [ val(meta), [ bam ] ]
    fasta               // channel: /path/to/genome.fasta
    sizes               // channel: /path/to/genome.sizes
    gff                 // channel: /path/to/genome.gff
    bed                 // channel: /path/to/primers.bed
    snpeff_db           // channel: /path/to/snpeff_db/
    snpeff_config       // channel: /path/to/snpeff.config
    ivar_multiqc_header // channel: /path/to/multiqc_header for ivar variants

    main:

    ch_versions = Channel.empty()

    //
    // Call variants
    //
    IVAR_VARIANTS (
        bam,
        fasta,
        gff,
        params.save_mpileup
    )
    ch_versions = ch_versions.mix(IVAR_VARIANTS.out.versions.first())

    //
    // Convert original iVar output to VCF, zip, sort and index
    //
    IVAR_VARIANTS_TO_VCF (
        IVAR_VARIANTS.out.tsv,
        ivar_multiqc_header
    )
    ch_versions = ch_versions.mix(IVAR_VARIANTS_TO_VCF.out.versions.first())

    BCFTOOLS_SORT (
        IVAR_VARIANTS_TO_VCF.out.vcf
    )
    ch_versions = ch_versions.mix(BCFTOOLS_SORT.out.versions)

    VCF_BGZIP_TABIX_STATS (
        BCFTOOLS_SORT.out.vcf
    )
    ch_versions = ch_versions.mix(VCF_BGZIP_TABIX_STATS.out.versions)

    //
    // Run downstream tools for variants QC
    //
    VARIANTS_QC (
        bam,
        VCF_BGZIP_TABIX_STATS.out.vcf,
        VCF_BGZIP_TABIX_STATS.out.stats,
        fasta,
        sizes,
        gff,
        bed,
        snpeff_db,
        snpeff_config
    )
    ch_versions = ch_versions.mix(VARIANTS_QC.out.versions)

    emit:
    tsv             = IVAR_VARIANTS.out.tsv           // channel: [ val(meta), [ tsv ] ]

    vcf_orig        = IVAR_VARIANTS_TO_VCF.out.vcf    // channel: [ val(meta), [ vcf ] ]
    log_out         = IVAR_VARIANTS_TO_VCF.out.log    // channel: [ val(meta), [ log ] ]
    multiqc_tsv     = IVAR_VARIANTS_TO_VCF.out.tsv    // channel: [ val(meta), [ tsv ] ]

    vcf             = VCF_BGZIP_TABIX_STATS.out.vcf   // channel: [ val(meta), [ vcf ] ]
    tbi             = VCF_BGZIP_TABIX_STATS.out.tbi   // channel: [ val(meta), [ tbi ] ]
    stats           = VCF_BGZIP_TABIX_STATS.out.stats // channel: [ val(meta), [ txt ] ]

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
