//
// Uncompress and prepare reference genome files
//

params.genome_options           = [:]
params.collapse_primers_options = [:]
params.snpeff_build_options     = [:]

include {
    GUNZIP as GUNZIP_FASTA
    GUNZIP as GUNZIP_GFF
    GUNZIP as GUNZIP_PRIMER_BED } from '../../modules/nf-core/software/gunzip/main' addParams( options: params.genome_options           )
include { GET_CHROM_SIZES       } from '../../modules/local/get_chrom_sizes'        addParams( options: params.genome_options           )
include { COLLAPSE_PRIMERS      } from '../../modules/local/collapse_primers'       addParams( options: params.collapse_primers_options )
include { SNPEFF_BUILD          } from '../../modules/local/snpeff_build'           addParams( options: params.snpeff_build_options     )

workflow PREPARE_GENOME {
    take:
    dummy_file

    main:

    //
    // Uncompress genome fasta file if required
    //
    if (params.fasta.endsWith('.gz')) {
        ch_fasta = GUNZIP_FASTA ( params.fasta ).gunzip
    } else {
        ch_fasta = file(params.fasta)
    }

    //
    // Uncompress GFF annotation file
    //
    if (params.gff) {
        if (params.gff.endsWith('.gz')) {
            ch_gff = GUNZIP_GFF ( params.gff ).gunzip
        } else {
            ch_gff = file(params.gff)
        }
    } else {
        ch_gff = dummy_file
    }

    //
    // Create chromosome sizes file
    //
    ch_chrom_sizes = Channel.empty()
    if (!params.skip_asciigenome) {
        ch_chrom_sizes = GET_CHROM_SIZES ( ch_fasta ).sizes
    }

    //
    // Uncompress primer BED file
    //
    ch_primer_bed = Channel.empty()
    if (params.primer_bed) {
        if (params.primer_bed.endsWith('.gz')) {
            ch_primer_bed = GUNZIP_PRIMER_BED ( params.primer_bed ).gunzip
        } else {
            ch_primer_bed = file(params.primer_bed)
        }
    }

    //
    // Generate collapsed BED file
    //
    ch_primer_collapsed_bed = Channel.empty()
    if (!params.skip_mosdepth) {
        ch_primer_collapsed_bed = COLLAPSE_PRIMERS ( ch_primer_bed, params.primer_left_suffix, params.primer_right_suffix )
    }

    //
    // Make snpEff database
    //
    ch_snpeff_db     = Channel.empty()
    ch_snpeff_config = Channel.empty()
    if (params.gff && !params.skip_snpeff) {
        SNPEFF_BUILD ( ch_fasta, ch_gff )
        ch_snpeff_db     = SNPEFF_BUILD.out.db
        ch_snpeff_config = SNPEFF_BUILD.out.config
    }

    emit:
    fasta                = ch_fasta                 // path: genome.fasta
    gff                  = ch_gff                   // path: genome.gff
    chrom_sizes          = ch_chrom_sizes           // path: genome.sizes
    primer_bed           = ch_primer_bed            // path: primer.bed
    primer_collapsed_bed = ch_primer_collapsed_bed  // path: primer.collapsed.bed
    snpeff_db            = ch_snpeff_db             // path: snpeff_db
    snpeff_config        = ch_snpeff_config         // path: snpeff.config
}
