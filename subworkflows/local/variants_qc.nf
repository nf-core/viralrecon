//
// Variant calling QC
//

include { ASCIIGENOME    } from '../../modules/local/asciigenome'
include { SNPEFF_SNPSIFT } from './snpeff_snpsift'

workflow VARIANTS_QC {
    take:
    bam           // channel: [ val(meta), [ bam ] ]
    vcf           // channel: [ val(meta), [ vcf ] ]
    stats         // channel: [ val(meta), [ bcftools_stats ] ]
    fasta         // channel: /path/to/genome.fasta
    sizes         // channel: /path/to/genome.sizes
    gff           // channel: /path/to/genome.gff
    bed           // channel: /path/to/primers.bed
    snpeff_db     // channel: /path/to/snpeff_db/
    snpeff_config // channel: /path/to/snpeff.config

    main:

    ch_versions = Channel.empty()

    //
    // Annotate variants
    //
    ch_snpeff_vcf   = Channel.empty()
    ch_snpeff_tbi   = Channel.empty()
    ch_snpeff_stats = Channel.empty()
    ch_snpeff_csv   = Channel.empty()
    ch_snpeff_txt   = Channel.empty()
    ch_snpeff_html  = Channel.empty()
    ch_snpsift_txt  = Channel.empty()
    if (params.gff && !params.skip_snpeff) {
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
        ch_versions     = ch_versions.mix(SNPEFF_SNPSIFT.out.versions)
    }

    //
    // Variant screenshots with ASCIIGenome
    //
    ch_asciigenome_pdf = Channel.empty()
    if (!params.skip_asciigenome) {
        bam
            .join(vcf, by: [0])
            .join(stats, by: [0])
            .map { meta, bam, vcf, stats ->
                if (WorkflowCommons.getNumVariantsFromBCFToolsStats(stats) > 0) {
                    return [ meta, bam, vcf ]
                }
            }
            .set { ch_asciigenome }

        ASCIIGENOME (
            ch_asciigenome,
            fasta,
            sizes,
            gff,
            bed,
            params.asciigenome_window_size,
            params.asciigenome_read_depth
        )
        ch_asciigenome_pdf = ASCIIGENOME.out.pdf
        ch_versions        = ch_versions.mix(ASCIIGENOME.out.versions.first())
    }

    emit:
    snpeff_vcf      = ch_snpeff_vcf      // channel: [ val(meta), [ vcf.gz ] ]
    snpeff_tbi      = ch_snpeff_tbi      // channel: [ val(meta), [ tbi ] ]
    snpeff_stats    = ch_snpeff_stats    // channel: [ val(meta), [ txt ] ]
    snpeff_csv      = ch_snpeff_csv      // channel: [ val(meta), [ csv ] ]
    snpeff_txt      = ch_snpeff_txt      // channel: [ val(meta), [ txt ] ]
    snpeff_html     = ch_snpeff_html     // channel: [ val(meta), [ html ] ]
    snpsift_txt     = ch_snpsift_txt     // channel: [ val(meta), [ txt ] ]

    asciigenome_pdf = ch_asciigenome_pdf // channel: [ val(meta), [ pdf ] ]

    versions        = ch_versions        // channel: [ versions.yml ]
}
