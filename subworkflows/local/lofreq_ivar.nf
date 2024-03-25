



//
// SUBWORKFLOW: Call variants with Lofreq, followed by ivar
//

//include { LOFREQ_INDELQUAL      } from '../modules/nf-core/lofreq_indelqual'
include { LOFREQ_CALLPARALLEL     } from '../../modules/nf-core/lofreq/callparallel/main'
include { BCFTOOLS_MERGE          } from '../../modules/nf-core/bcftools/merge/main'
include { VARIANTS_IVAR           } from '../subworkflows/local/variants_ivar'


workflow LOFREQ_IVAR {
    take:
    lofreq_bam          // channel: [ val(meta), [ bam ], [ bai ] ]
    bam                 // channel: [ val(meta), [ bam ] ] // for ivar
    fasta               // channel: /path/to/genome.fasta
    fai                 // channel: /path/to/genome.fai
    sizes               // channel: /path/to/genome.sizes
    gff                 // channel: /path/to/genome.gff
    bed                 // channel: /path/to/primers.bed
    snpeff_db           // channel: /path/to/snpeff_db/
    snpeff_config       // channel: /path/to/snpeff.config
    ivar_multiqc_header // channel: /path/to/multiqc_header for ivar variants

    main:

    ch_lofreq_vcf                    = Channel.empty()
    ch_vcf                           = Channel.empty()  // ivar
    //ch_versions                      = Channel.empty()
    //ch_lofreq_tbi                    = Channel.empty()

    //
    // SUBWORKFLOW: Call variants with Lofreq, followed by ivar
    //
    // Call Lofreq variants
    LOFREQ_CALLPARALLEL (
        lofreq_bam,
        fasta,
        fai
    )
    //ch_lofreq_vcf               = LOFREQ_CALLPARALLEL.out.vcf
    ch_versions                 = ch_versions.mix(LOFREQ_CALLPARALLEL.out.versions)

    //
    // Call ivar variants by invoking ivar subworkflow
    VARIANTS_IVAR (
        ivar_bam,
        fasta,
        fai,
        sizes,
        gff,
        bed,
        snpeff_db,
        snpeff_config,
        ivar_multiqc_header
    )
    ch_versions                 = ch_versions.mix(VARIANTS_IVAR.out.versions)

    // Merge vcf files w BCFTOOLS_MERGE module
        BCFTOOLS_MERGE(
            vcfs,
            fasta,
            fai
        )
    ch_versions                 = ch_versions.mix(BCFTOOLS_MERGE.out.versions)




    emit:
    // from VARIANTS_IVAR workflow
    tsv             = ch_ivar_tsv                     // channel: [ val(meta), [ tsv ] ]

    vcf_orig        = IVAR_VARIANTS_TO_VCF.out.vcf    // channel: [ val(meta), [ vcf ] ]
    log_out         = IVAR_VARIANTS_TO_VCF.out.log    // channel: [ val(meta), [ log ] ]
    multiqc_tsv     = IVAR_VARIANTS_TO_VCF.out.tsv    // channel: [ val(meta), [ tsv ] ]

    vcf             = BCFTOOLS_SORT.out.vcf           // channel: [ val(meta), [ vcf ] ]
    tbi             = VCF_TABIX_STATS.out.tbi         // channel: [ val(meta), [ tbi ] ]
    csi             = VCF_TABIX_STATS.out.csi         // channel: [ val(meta), [ csi ] ]
    stats           = VCF_TABIX_STATS.out.stats       // channel: [ val(meta), [ txt ] ]

    snpeff_vcf      = VARIANTS_QC.out.snpeff_vcf      // channel: [ val(meta), [ vcf.gz ] ]
    snpeff_tbi      = VARIANTS_QC.out.snpeff_tbi      // channel: [ val(meta), [ tbi ] ]
    snpeff_stats    = VARIANTS_QC.out.snpeff_stats    // channel: [ val(meta), [ txt ] ]
    snpeff_csv      = VARIANTS_QC.out.snpeff_csv      // channel: [ val(meta), [ csv ] ]
    snpeff_txt      = VARIANTS_QC.out.snpeff_txt      // channel: [ val(meta), [ txt ] ]
    snpeff_html     = VARIANTS_QC.out.snpeff_html     // channel: [ val(meta), [ html ] ]
    snpsift_txt     = VARIANTS_QC.out.snpsift_txt     // channel: [ val(meta), [ txt ] ]

    asciigenome_pdf = VARIANTS_QC.out.asciigenome_pdf // channel: [ val(meta), [ pdf ] ]

    // added from proposed subworkflow
    lofreq_vcf      = LOFREQ_CALLPARALLEL.out.vcf     // channel: [ val(meta), [ vcf ] ]
    //lofreq_tbi      = LOFREQ_CALLPARALLEL.out.tbi     // channel: [ val(meta), [ tbi ] ]
    merged_vcf      = BCFTOOLS_MERGE.out.merged_variants  // channel: [ val(meta), [ vcf ] ]

    versions        = ch_versions                     // channel: [ versions.yml ]
}