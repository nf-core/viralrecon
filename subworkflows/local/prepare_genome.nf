/*
 * Uncompress and prepare reference genome files
*/

params.genome_options        = [:]
params.index_options         = [:]
params.bowtie2_index_options = [:]
params.snpeff_build_options  = [:]

include {
    GUNZIP as GUNZIP_FASTA
    GUNZIP as GUNZIP_GFF               
    GUNZIP as GUNZIP_AMPLICON_BED
    GUNZIP as GUNZIP_AMPLICON_FASTA    } from '../../modules/local/gunzip'        addParams( options: params.genome_options        )
include { UNTAR as UNTAR_BOWTIE2_INDEX } from '../../modules/local/untar'         addParams( options: params.bowtie2_index_options )
include { UNTAR as UNTAR_KRAKEN2_DB    } from '../../modules/local/untar'         addParams( options: params.index_options         )
include { BOWTIE2_BUILD                } from '../../modules/local/bowtie2_build' addParams( options: params.bowtie2_index_options )
include { SNPEFF_BUILD                 } from '../../modules/local/snpeff_build'  addParams( options: params.snpeff_build_options  )

workflow PREPARE_GENOME {
    take:
    ch_dummy_file

    main:
    /*
     * Uncompress genome fasta file if required
     */
    if (params.fasta.endsWith('.gz')) {
        ch_fasta = GUNZIP_FASTA ( params.fasta ).gunzip
    } else {
        ch_fasta = file(params.fasta)
    }

    /*
     * Uncompress GFF annotation file
     */
    if (params.gff) {
        if (params.gff.endsWith('.gz')) {
            ch_gff = GUNZIP_GFF ( params.gff ).gunzip
        } else {
            ch_gff = file(params.gff)
        }
    } else {
        ch_gff = ch_dummy_file
    }

    /*
     * Prepare reference files required for variant calling
     */
    ch_amplicon_bed    = Channel.empty()
    ch_bowtie2_index   = Channel.empty()
    ch_bowtie2_version = Channel.empty()
    if (!params.skip_variants) {

        if (params.amplicon_bed) {
            if (params.amplicon_bed.endsWith('.gz')) {
                ch_amplicon_bed = GUNZIP_AMPLICON_BED ( params.amplicon_bed ).gunzip
            } else {
                ch_amplicon_bed = file(params.amplicon_bed)
            }
        }

        if (params.bowtie2_index) {
            if (params.bowtie2_index.endsWith('.tar.gz')) {
                ch_bowtie2_index = UNTAR_BOWTIE2_INDEX ( params.bowtie2_index ).untar
            } else {
                ch_bowtie2_index = file(params.bowtie2_index)
            }
        } else {
            ch_bowtie2_index   = BOWTIE2_BUILD ( ch_fasta ).index
            ch_bowtie2_version = BOWTIE2_BUILD.out.version
        }
    }

    /*
     * Prepare reference files required for de novo assembly
     */
    ch_amplicon_fasta = Channel.empty()
    ch_kraken2_db     = Channel.empty()
    if (!params.skip_assembly) {

        if (params.amplicon_fasta) {
            if (params.amplicon_fasta.endsWith('.gz')) {
                ch_amplicon_fasta = GUNZIP_AMPLICON_FASTA ( params.amplicon_fasta ).gunzip
            } else {
                ch_amplicon_fasta = file(params.amplicon_fasta)
            }
        }

        if (!params.skip_kraken2) {
            if (params.kraken2_db) {
                if (params.kraken2_db.endsWith('.tar.gz')) {
                    ch_kraken2_db = UNTAR_KRAKEN2_DB ( params.kraken2_db ).untar
                } else {
                    ch_kraken2_db = file(params.kraken2_db)
                }
            } else {
                //BUILD KRAKEN2 DATABASE HERE
                ch_kraken2_db = ch_dummy_file
            }
        }
    }

    /*
     * Make snpEff database
     */
    ch_snpeff_db     = Channel.empty()
    ch_snpeff_config = Channel.empty()
    if ((!params.skip_variants || !params.skip_assembly) && params.gff && !params.skip_snpeff) {
        SNPEFF_BUILD ( ch_fasta, ch_gff )
        ch_snpeff_db     = SNPEFF_BUILD.out.db
        ch_snpeff_config = SNPEFF_BUILD.out.config
    }
    
    emit:
    fasta           = ch_fasta            // path: genome.fasta
    gff             = ch_gff              // path: genome.gff
    amplicon_bed    = ch_amplicon_bed     // path: amplicon.bed
    amplicon_fasta  = ch_amplicon_fasta   // path: amplicon.fasta
    bowtie2_index   = ch_bowtie2_index    // path: bowtie2/index/
    kraken2_db      = ch_kraken2_db       // path: kraken2_db/
    snpeff_db       = ch_snpeff_db        // path: snpeff_db
    snpeff_config   = ch_snpeff_config    // path: snpeff.config
    bowtie2_version = ch_bowtie2_version  // path: *.version.txt
}
