/*
 * Downstream analysis for assembly scaffolds
 */

params.bandage_options      = [:]
params.blastn_options       = [:]
params.abacas_options       = [:]
params.plasmidid_options    = [:]
params.quast_options        = [:]
params.snpeff_options       = [:]
params.snpeff_bgzip_options = [:]
params.snpeff_tabix_options = [:]
params.snpeff_stats_options = [:]
params.snpsift_options      = [:]

include { BANDAGE        } from '../../modules/local/bandage'      addParams( options: params.bandage_options   ) 
include { BLAST_BLASTN   } from '../../modules/local/blast_blastn' addParams( options: params.blastn_options    ) 
include { ABACAS         } from '../../modules/local/abacas'       addParams( options: params.abacas_options    )
include { PLASMIDID      } from '../../modules/local/plasmidid'    addParams( options: params.plasmidid_options )
include { QUAST          } from '../../modules/local/quast'        addParams( options: params.quast_options     )
include { SNPEFF_SNPSIFT } from './snpeff_snpsift'                 addParams( snpeff_options: params.snpeff_options, snpsift_options: params.snpsift_options, bgzip_options: params.snpeff_bgzip_options, tabix_options: params.snpeff_tabix_options, stats_options: params.snpeff_stats_options )

workflow ASSEMBLY_QC {
    take:
    scaffolds     // channel: [ val(meta), [ scaffolds ] ]
    graphs        // channel: [ val(meta), [ graphs ] ]
    fasta         // channel: /path/to/genome.fasta
    gff           // channel: /path/to/genome.gff
    blast_db      // channel: /path/to/blast_db/
    snpeff_db     // channel: /path/to/snpeff_db/
    snpeff_config // channel: /path/to/snpeff_db/
    
    main:
    /*
     * Run blastn on assembly scaffolds
     */
    ch_blast_txt     = Channel.empty()
    ch_blast_version = Channel.empty()
    if (!params.skip_blast) {
        BLAST_BLASTN ( scaffolds, blast_db )
        ch_blast_txt     = BLAST_BLASTN.out.txt
        ch_blast_version = BLAST_BLASTN.out.version
    }

    /*
     * Assembly QC across all samples with QUAST
     */
    ch_quast_results = Channel.empty()
    ch_quast_tsv     = Channel.empty()
    ch_quast_version = Channel.empty()
    if (!params.skip_assembly_quast) {
        QUAST ( scaffolds.collect{ it[1] }, fasta, gff )
        ch_quast_results = QUAST.out.results
        ch_quast_tsv     = QUAST.out.tsv
        ch_quast_version = QUAST.out.version
    }

    /*
     * Generate assembly visualisation with Bandage
     */
    ch_bandage_png     = Channel.empty()
    ch_bandage_svg     = Channel.empty()
    ch_bandage_version = Channel.empty()
    if (!params.skip_bandage) {
        BANDAGE ( graphs )
        ch_bandage_version = BANDAGE.out.version
        ch_bandage_png     = BANDAGE.out.png
        ch_bandage_svg     = BANDAGE.out.svg
    }

    /*
     * Contiguate assembly with ABACAS
     */
    ch_abacas_results = Channel.empty()
    ch_abacas_version = Channel.empty()
    // if (!params.skip_abacas) {
    //     ABACAS ( scaffolds, fasta )
    //     ch_abacas_results = ABACAS.out.results
    //     ch_abacas_version = ABACAS.out.version
    // }

    /*
     * Assembly report with PlasmidID
     */
    ch_plasmidid_results = Channel.empty()
    ch_plasmidid_version = Channel.empty()
    // if (!params.skip_plasmidid) {
    //     PLASMIDID ( scaffolds, fasta )
    //     ch_plasmidid_results = PLASMIDID.out.results
    //     ch_plasmidid_version = PLASMIDID.out.version
    // }

    // if (!params.skip_variants_snpeff) {
    //     SNPEFF_SNPSIFT ( vcf, snpeff_db, snpeff_config, fasta )
    // }

    emit:
    blast_txt         = ch_blast_txt         // channel: [ val(meta), [ txt ] ]
    blast_version     = ch_blast_version     //    path: *.version.txt

    quast_results     = ch_quast_results     // channel: [ val(meta), [ results ] ]
    quast_tsv         = ch_quast_tsv         // channel: [ val(meta), [ tsv ] ]
    quast_version     = ch_quast_version     //    path: *.version.txt

    bandage_png       = ch_bandage_png       // channel: [ val(meta), [ png ] ]
    bandage_svg       = ch_bandage_svg       // channel: [ val(meta), [ svg ] ]
    bandage_version   = ch_bandage_version   //    path: *.version.txt

    abacas_results    = ch_abacas_results    // channel: [ val(meta), [ results ] ]
    abacas_version    = ch_abacas_version    //    path: *.version.txt

    plasmidid_results = ch_plasmidid_results // channel: [ val(meta), [ results ] ]
    plasmidid_version = ch_plasmidid_version //    path: *.version.txt

}

