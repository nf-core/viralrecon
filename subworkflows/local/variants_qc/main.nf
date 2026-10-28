//
// Variant calling QC
//

include { SNPEFF_SNPSIFT } from '../snpeff_snpsift'

workflow VARIANTS_QC {
    take:
    vcf           // channel: [ val(meta), [ vcf ] ]
    fasta         // channel: /path/to/genome.fasta
    gff           // channel: /path/to/genome.gff
    snpeff_db     // channel: /path/to/snpeff_db/
    snpeff_config // channel: /path/to/snpeff.config

    main:

    //
    // Annotate variants
    //
    ch_snpeff_vcf   = channel.empty()
    ch_snpeff_tbi   = channel.empty()
    ch_snpeff_stats = channel.empty()
    ch_snpeff_csv   = channel.empty()
    ch_snpeff_txt   = channel.empty()
    ch_snpeff_html  = channel.empty()
    ch_snpsift_txt  = channel.empty()
    if (gff && !params.skip_snpeff) {
        SNPEFF_SNPSIFT (
            vcf,
            snpeff_db,
            snpeff_config,
            fasta
        )
        ch_snpeff_vcf   = SNPEFF_SNPSIFT.out.vcf
        ch_snpeff_tbi   = SNPEFF_SNPSIFT.out.tbi
        ch_snpeff_stats = SNPEFF_SNPSIFT.out.stats
        ch_snpeff_csv   = SNPEFF_SNPSIFT.out.csv
        ch_snpeff_txt   = SNPEFF_SNPSIFT.out.txt
        ch_snpeff_html  = SNPEFF_SNPSIFT.out.html
        ch_snpsift_txt  = SNPEFF_SNPSIFT.out.snpsift_txt
    }

    emit:
    snpeff_vcf      = ch_snpeff_vcf      // channel: [ val(meta), [ vcf.gz ] ]
    snpeff_tbi      = ch_snpeff_tbi      // channel: [ val(meta), [ tbi ] ]
    snpeff_stats    = ch_snpeff_stats    // channel: [ val(meta), [ txt ] ]
    snpeff_csv      = ch_snpeff_csv      // channel: [ val(meta), [ csv ] ]
    snpeff_txt      = ch_snpeff_txt      // channel: [ val(meta), [ txt ] ]
    snpeff_html     = ch_snpeff_html     // channel: [ val(meta), [ html ] ]
    snpsift_txt     = ch_snpsift_txt     // channel: [ val(meta), [ txt ] ]
}
