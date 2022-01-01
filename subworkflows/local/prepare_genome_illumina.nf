//
// Uncompress and prepare reference genome files
//

params.genome_options            = [:]
params.index_options             = [:]
params.db_options                = [:]
params.bowtie2_build_options     = [:]
params.collapse_primers_options  = [:]
params.bedtools_getfasta_options = [:]
params.snpeff_build_options      = [:]
params.makeblastdb_options       = [:]
params.kraken2_build_options     = [:]

include {
    GUNZIP as GUNZIP_FASTA
    GUNZIP as GUNZIP_GFF
    GUNZIP as GUNZIP_PRIMER_BED
    GUNZIP as GUNZIP_PRIMER_FASTA      } from '../../modules/nf-core/modules/gunzip/main'            addParams( options: params.genome_options            )
include { UNTAR as UNTAR_BOWTIE2_INDEX } from '../../modules/nf-core/modules/untar/main'             addParams( options: params.index_options             )
include { UNTAR as UNTAR_KRAKEN2_DB    } from '../../modules/nf-core/modules/untar/main'             addParams( options: params.db_options                )
include { UNTAR as UNTAR_BLAST_DB      } from '../../modules/nf-core/modules/untar/main'             addParams( options: params.db_options                )
include { BOWTIE2_BUILD                } from '../../modules/nf-core/modules/bowtie2/build/main'     addParams( options: params.bowtie2_build_options     )
include { BLAST_MAKEBLASTDB            } from '../../modules/nf-core/modules/blast/makeblastdb/main' addParams( options: params.makeblastdb_options       )
include { BEDTOOLS_GETFASTA            } from '../../modules/nf-core/modules/bedtools/getfasta/main' addParams( options: params.bedtools_getfasta_options )
include { GET_CHROM_SIZES              } from '../../modules/local/get_chrom_sizes'                  addParams( options: params.genome_options            )
include { COLLAPSE_PRIMERS             } from '../../modules/local/collapse_primers'                 addParams( options: params.collapse_primers_options  )
include { KRAKEN2_BUILD                } from '../../modules/local/kraken2_build'                    addParams( options: params.kraken2_build_options     )
include { SNPEFF_BUILD                 } from '../../modules/local/snpeff_build'                     addParams( options: params.snpeff_build_options      )

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
    // Prepare reference files required for variant calling
    //
    ch_kraken2_db = Channel.empty()
    if (!params.skip_kraken2) {
        if (params.kraken2_db) {
            if (params.kraken2_db.endsWith('.tar.gz')) {
                ch_kraken2_db = UNTAR_KRAKEN2_DB ( params.kraken2_db ).untar
            } else {
                ch_kraken2_db = file(params.kraken2_db)
            }
        } else {
            ch_kraken2_db = KRAKEN2_BUILD ( params.kraken2_db_name ).db
        }
    }

    //
    // Prepare files required for amplicon data
    //
    ch_primer_bed           = Channel.empty()
    ch_primer_fasta         = Channel.empty()
    ch_primer_collapsed_bed = Channel.empty()
    if (params.protocol == 'amplicon') {
        if (params.primer_bed) {
            if (params.primer_bed.endsWith('.gz')) {
                ch_primer_bed = GUNZIP_PRIMER_BED ( params.primer_bed ).gunzip
            } else {
                ch_primer_bed = file(params.primer_bed)
            }
        }

        if (!params.skip_variants && !params.skip_mosdepth) {
            ch_primer_collapsed_bed = COLLAPSE_PRIMERS ( ch_primer_bed, params.primer_left_suffix, params.primer_right_suffix )
        }

        if (!params.skip_assembly && !params.skip_cutadapt) {
            if (params.primer_fasta) {
                if (params.primer_fasta.endsWith('.gz')) {
                    ch_primer_fasta = GUNZIP_PRIMER_FASTA ( params.primer_fasta ).gunzip
                } else {
                    ch_primer_fasta = file(params.primer_fasta)
                }
            } else {
                ch_primer_fasta = BEDTOOLS_GETFASTA ( ch_primer_bed, ch_fasta ).fasta
            }
        }
    }

    //
    // Prepare reference files required for variant calling
    //
    ch_bowtie2_index = Channel.empty()
    if (!params.skip_variants) {
        if (params.bowtie2_index) {
            if (params.bowtie2_index.endsWith('.tar.gz')) {
                ch_bowtie2_index = UNTAR_BOWTIE2_INDEX ( params.bowtie2_index ).untar
            } else {
                ch_bowtie2_index = file(params.bowtie2_index)
            }
        } else {
            ch_bowtie2_index = BOWTIE2_BUILD ( ch_fasta ).index
        }
    }

    //
    // Prepare reference files required for de novo assembly
    //
    ch_blast_db = Channel.empty()
    if (!params.skip_assembly) {
        if (!params.skip_blast) {
            if (params.blast_db) {
                if (params.blast_db.endsWith('.tar.gz')) {
                    ch_blast_db = UNTAR_BLAST_DB ( params.blast_db ).untar
                } else {
                    ch_blast_db = file(params.blast_db)
                }
            } else {
                ch_blast_db = BLAST_MAKEBLASTDB ( ch_fasta ).db
            }
        }
    }

    //
    // Make snpEff database
    //
    ch_snpeff_db     = Channel.empty()
    ch_snpeff_config = Channel.empty()
    if (!params.skip_variants && params.gff && !params.skip_snpeff) {
        SNPEFF_BUILD ( ch_fasta, ch_gff )
        ch_snpeff_db     = SNPEFF_BUILD.out.db
        ch_snpeff_config = SNPEFF_BUILD.out.config
    }

    emit:
    fasta                = ch_fasta                 // path: genome.fasta
    gff                  = ch_gff                   // path: genome.gff
    chrom_sizes          = ch_chrom_sizes           // path: genome.sizes
    bowtie2_index        = ch_bowtie2_index         // path: bowtie2/index/
    primer_bed           = ch_primer_bed            // path: primer.bed
    primer_collapsed_bed = ch_primer_collapsed_bed  // path: primer.collapsed.bed
    primer_fasta         = ch_primer_fasta          // path: primer.fasta
    blast_db             = ch_blast_db              // path: blast_db/
    kraken2_db           = ch_kraken2_db            // path: kraken2_db/
    snpeff_db            = ch_snpeff_db             // path: snpeff_db
    snpeff_config        = ch_snpeff_config         // path: snpeff.config
}
