//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA      } from '../../modules/nf-core/modules/gunzip/main'
include { GUNZIP as GUNZIP_GFF        } from '../../modules/nf-core/modules/gunzip/main'
include { GUNZIP as GUNZIP_PRIMER_BED } from '../../modules/nf-core/modules/gunzip/main'
include { CUSTOM_GETCHROMSIZES        } from '../../modules/nf-core/modules/custom/getchromsizes/main'
include { COLLAPSE_PRIMERS            } from '../../modules/local/collapse_primers'
include { SNPEFF_BUILD                } from '../../modules/local/snpeff_build'

workflow PREPARE_GENOME {
    main:

    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    if (params.fasta.endsWith('.gz')) {
        GUNZIP_FASTA (
            [ [:], params.fasta ]
        )
        ch_fasta    = GUNZIP_FASTA.out.gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = file(params.fasta)
    }

    //
    // Uncompress GFF annotation file
    //
    if (params.gff) {
        if (params.gff.endsWith('.gz')) {
            GUNZIP_GFF (
                [ [:], params.gff ]
            )
            ch_gff      = GUNZIP_GFF.out.gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
        } else {
            ch_gff = file(params.gff)
        }
    } else {
        ch_gff = []
    }

    //
    // Create chromosome sizes file
    //
    ch_chrom_sizes = Channel.empty()
    if (!params.skip_asciigenome) {
        CUSTOM_GETCHROMSIZES (
            ch_fasta
        )
        ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes
        ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)
    }

    //
    // Uncompress primer BED file
    //
    ch_primer_bed = Channel.empty()
    if (params.primer_bed) {
        if (params.primer_bed.endsWith('.gz')) {
            GUNZIP_PRIMER_BED (
                [ [:], params.primer_bed ]
            )
            ch_primer_bed = GUNZIP_PRIMER_BED.out.gunzip.map { it[1] }
            ch_versions   = ch_versions.mix(GUNZIP_PRIMER_BED.out.versions)
        } else {
            ch_primer_bed = file(params.primer_bed)
        }
    }

    //
    // Generate collapsed BED file
    //
    ch_primer_collapsed_bed = Channel.empty()
    if (!params.skip_mosdepth) {
        COLLAPSE_PRIMERS (
            ch_primer_bed,
            params.primer_left_suffix,
            params.primer_right_suffix
        )
        ch_primer_collapsed_bed = COLLAPSE_PRIMERS.out.bed
        ch_versions             = ch_versions.mix(COLLAPSE_PRIMERS.out.versions)
    }

    //
    // Make snpEff database
    //
    ch_snpeff_db     = Channel.empty()
    ch_snpeff_config = Channel.empty()
    if (params.gff && !params.skip_snpeff) {
        SNPEFF_BUILD (
            ch_fasta,
            ch_gff
        )
        ch_snpeff_db     = SNPEFF_BUILD.out.db
        ch_snpeff_config = SNPEFF_BUILD.out.config
        ch_versions      = ch_versions.mix(SNPEFF_BUILD.out.versions)
    }

    emit:
    fasta                = ch_fasta                 // path: genome.fasta
    gff                  = ch_gff                   // path: genome.gff
    chrom_sizes          = ch_chrom_sizes           // path: genome.sizes
    primer_bed           = ch_primer_bed            // path: primer.bed
    primer_collapsed_bed = ch_primer_collapsed_bed  // path: primer.collapsed.bed
    snpeff_db            = ch_snpeff_db             // path: snpeff_db
    snpeff_config        = ch_snpeff_config         // path: snpeff.config

    versions             = ch_versions              // channel: [ versions.yml ]
}
