/*
 * Uncompress and prepare reference genome files
*/

params.genome_options           = [:]
params.clone_scheme_options     = [:]
params.collapse_primers_options = [:]
params.snpeff_build_options     = [:]

include { GUNZIP             } from '../../modules/nf-core/software/gunzip/main' addParams( options: params.genome_options           )
include { ARTIC_CLONE_SCHEME } from '../../modules/local/artic_clone_scheme'     addParams( options: params.clone_scheme_options     )
include { COLLAPSE_PRIMERS   } from '../../modules/local/collapse_primers'       addParams( options: params.collapse_primers_options )
include { SNPEFF_BUILD       } from '../../modules/local/snpeff_build'           addParams( options: params.snpeff_build_options     )

workflow PREPARE_GENOME {
    take:
    dummy_file

    main:
    /*
     * Fetch ARTIC primer scheme and reference genome if required
     */
    if (!params.artic_scheme_dir) {
        ARTIC_CLONE_SCHEME (
            params.artic_scheme,
            params.artic_scheme_version,
            params.artic_scheme_repo_url
        )
        ch_scheme_dir = ARTIC_CLONE_SCHEME.out.scheme
        ch_primer_bed = ARTIC_CLONE_SCHEME.out.bed
        ch_fasta      = ARTIC_CLONE_SCHEME.out.fasta
    } else {
        ch_scheme_dir = file(params.artic_scheme_dir, checkIfExists: true)
        // ch_fasta
        // ch_primer_bed
    }

    /*
     * Uncompress GFF annotation file
     */
    if (params.gff) {
        if (params.gff.endsWith('.gz')) {
            ch_gff = GUNZIP ( params.gff ).gunzip
        } else {
            ch_gff = file(params.gff)
        }
    } else {
        ch_gff = dummy_file
    }

    /*
     * Prepare reference files required for variant calling
     */
    ch_primer_collapsed_bed = Channel.empty()
    if (!params.skip_variants) {
        if (!params.skip_mosdepth) {
            ch_primer_collapsed_bed = COLLAPSE_PRIMERS ( ch_primer_bed, params.primer_left_suffix, params.primer_right_suffix )
        }
    }

    /*
     * Make snpEff database
     */
    ch_snpeff_db     = Channel.empty()
    ch_snpeff_config = Channel.empty()
    if (params.gff && !params.skip_variants_snpeff) {
        SNPEFF_BUILD ( ch_fasta, ch_gff )
        ch_snpeff_db     = SNPEFF_BUILD.out.db
        ch_snpeff_config = SNPEFF_BUILD.out.config
    }
    
    emit:
    fasta                = ch_fasta                 // path: genome.fasta
    gff                  = ch_gff                   // path: genome.gff
    scheme_dir           = ch_scheme_dir            // path: scheme_dir
    primer_bed           = ch_primer_bed            // path: primer.bed
    primer_collapsed_bed = ch_primer_collapsed_bed  // path: primer.collapsed.bed
    snpeff_db            = ch_snpeff_db             // path: snpeff_db
    snpeff_config        = ch_snpeff_config         // path: snpeff.config
}
