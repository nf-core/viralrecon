//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA        } from '../../modules/nf-core/modules/gunzip/main'
include { GUNZIP as GUNZIP_GFF          } from '../../modules/nf-core/modules/gunzip/main'
include { GUNZIP as GUNZIP_PRIMER_BED   } from '../../modules/nf-core/modules/gunzip/main'
include { GUNZIP as GUNZIP_PRIMER_FASTA } from '../../modules/nf-core/modules/gunzip/main'
include { UNTAR as UNTAR_BOWTIE2_INDEX  } from '../../modules/nf-core/modules/untar/main'
include { UNTAR as UNTAR_NEXTCLADE_DB   } from '../../modules/nf-core/modules/untar/main'
include { UNTAR as UNTAR_KRAKEN2_DB     } from '../../modules/nf-core/modules/untar/main'
include { UNTAR as UNTAR_BLAST_DB       } from '../../modules/nf-core/modules/untar/main'
include { BOWTIE2_BUILD                 } from '../../modules/nf-core/modules/bowtie2/build/main'
include { BLAST_MAKEBLASTDB             } from '../../modules/nf-core/modules/blast/makeblastdb/main'
include { BEDTOOLS_GETFASTA             } from '../../modules/nf-core/modules/bedtools/getfasta/main'
include { CUSTOM_GETCHROMSIZES          } from '../../modules/nf-core/modules/custom/getchromsizes/main'
include { NEXTCLADE_DATASETGET          } from '../../modules/nf-core/modules/nextclade/datasetget/main'
include { COLLAPSE_PRIMERS              } from '../../modules/local/collapse_primers'
include { KRAKEN2_BUILD                 } from '../../modules/local/kraken2_build'
include { SNPEFF_BUILD                  } from '../../modules/local/snpeff_build'

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
    ch_fai         = Channel.empty()
    ch_chrom_sizes = Channel.empty()
    if (params.protocol == 'amplicon' || !params.skip_asciigenome) {
        CUSTOM_GETCHROMSIZES (
            ch_fasta
        )
        ch_fai         = CUSTOM_GETCHROMSIZES.out.fai
        ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes
        ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)
    }

    //
    // Prepare reference files required for variant calling
    //
    ch_kraken2_db = Channel.empty()
    if (!params.skip_kraken2) {
        if (params.kraken2_db) {
            if (params.kraken2_db.endsWith('.tar.gz')) {
                UNTAR_KRAKEN2_DB (
                    [ [:], params.kraken2_db ]
                )
                ch_kraken2_db = UNTAR_KRAKEN2_DB.out.untar.map { it[1] }
                ch_versions   = ch_versions.mix(UNTAR_KRAKEN2_DB.out.versions)
            } else {
                ch_kraken2_db = file(params.kraken2_db)
            }
        } else {
            KRAKEN2_BUILD (
                params.kraken2_db_name
            )
            ch_kraken2_db = KRAKEN2_BUILD.out.db
            ch_versions   = ch_versions.mix(KRAKEN2_BUILD.out.versions)
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
                GUNZIP_PRIMER_BED (
                    [ [:], params.primer_bed ]
                )
                ch_primer_bed = GUNZIP_PRIMER_BED.out.gunzip.map { it[1] }
                ch_versions   = ch_versions.mix(GUNZIP_PRIMER_BED.out.versions)
            } else {
                ch_primer_bed = file(params.primer_bed)
            }
        }

        if (!params.skip_variants && !params.skip_mosdepth) {
            COLLAPSE_PRIMERS (
                ch_primer_bed,
                params.primer_left_suffix,
                params.primer_right_suffix
            )
            ch_primer_collapsed_bed = COLLAPSE_PRIMERS.out.bed
            ch_versions             = ch_versions.mix(COLLAPSE_PRIMERS.out.versions)
        }

        if (!params.skip_assembly && !params.skip_cutadapt) {
            if (params.primer_fasta) {
                if (params.primer_fasta.endsWith('.gz')) {
                    GUNZIP_PRIMER_FASTA (
                        [ [:], params.primer_fasta ]
                    )
                    ch_primer_fasta = GUNZIP_PRIMER_FASTA.out.gunzip.map { it[1] }
                    ch_versions     = ch_versions.mix(GUNZIP_PRIMER_FASTA.out.versions)
                } else {
                    ch_primer_fasta = file(params.primer_fasta)
                }
            } else {
                BEDTOOLS_GETFASTA (
                    ch_primer_bed,
                    ch_fasta
                )
                ch_primer_fasta = BEDTOOLS_GETFASTA.out.fasta
                ch_versions     = ch_versions.mix(BEDTOOLS_GETFASTA.out.versions)
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
                UNTAR_BOWTIE2_INDEX (
                    [ [:], params.bowtie2_index ]
                )
                ch_bowtie2_index = UNTAR_BOWTIE2_INDEX.out.untar.map { it[1] }
                ch_versions      = ch_versions.mix(UNTAR_BOWTIE2_INDEX.out.versions)
            } else {
                ch_bowtie2_index = file(params.bowtie2_index)
            }
        } else {
            BOWTIE2_BUILD (
                ch_fasta
            )
            ch_bowtie2_index = BOWTIE2_BUILD.out.index
            ch_versions      = ch_versions.mix(BOWTIE2_BUILD.out.versions)
        }
    }

    //
    // Prepare Nextclade dataset
    //
    ch_nextclade_db = Channel.empty()
    if (!params.skip_consensus && !params.skip_nextclade) {
        if (params.nextclade_dataset) {
            if (params.nextclade_dataset.endsWith('.tar.gz')) {
                UNTAR_NEXTCLADE_DB (
                    [ [:], params.nextclade_dataset ]
                )
                ch_nextclade_db = UNTAR_NEXTCLADE_DB.out.untar.map { it[1] }
                ch_versions     = ch_versions.mix(UNTAR_NEXTCLADE_DB.out.versions)
            } else {
                ch_nextclade_db = file(params.nextclade_dataset)
            }
        } else if (params.nextclade_dataset_name) {
            NEXTCLADE_DATASETGET (
                params.nextclade_dataset_name,
                params.nextclade_dataset_reference,
                params.nextclade_dataset_tag
            )
            ch_nextclade_db = NEXTCLADE_DATASETGET.out.dataset
            ch_versions     = ch_versions.mix(NEXTCLADE_DATASETGET.out.versions)
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
                    UNTAR_BLAST_DB (
                        [ [:], params.blast_db ]
                    )
                    ch_blast_db = UNTAR_BLAST_DB.out.untar.map { it[1] }
                    ch_versions = ch_versions.mix(UNTAR_BLAST_DB.out.versions)
                } else {
                    ch_blast_db = file(params.blast_db)
                }
            } else {
                BLAST_MAKEBLASTDB (
                    ch_fasta
                )
                ch_blast_db = BLAST_MAKEBLASTDB.out.db
                ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)
            }
        }
    }

    //
    // Make snpEff database
    //
    ch_snpeff_db     = Channel.empty()
    ch_snpeff_config = Channel.empty()
    if (!params.skip_variants && params.gff && !params.skip_snpeff) {
        SNPEFF_BUILD (
            ch_fasta,
            ch_gff
        )
        ch_snpeff_db     = SNPEFF_BUILD.out.db
        ch_snpeff_config = SNPEFF_BUILD.out.config
        ch_versions      = ch_versions.mix(SNPEFF_BUILD.out.versions)
    }

    emit:
    fasta                = ch_fasta                // path: genome.fasta
    gff                  = ch_gff                  // path: genome.gff
    fai                  = ch_fai                  // path: genome.fai
    chrom_sizes          = ch_chrom_sizes          // path: genome.sizes
    bowtie2_index        = ch_bowtie2_index        // path: bowtie2/index/
    primer_bed           = ch_primer_bed           // path: primer.bed
    primer_collapsed_bed = ch_primer_collapsed_bed // path: primer.collapsed.bed
    primer_fasta         = ch_primer_fasta         // path: primer.fasta
    nextclade_db         = ch_nextclade_db         // path: nextclade_db
    blast_db             = ch_blast_db             // path: blast_db/
    kraken2_db           = ch_kraken2_db           // path: kraken2_db/
    snpeff_db            = ch_snpeff_db            // path: snpeff_db
    snpeff_config        = ch_snpeff_config        // path: snpeff.config

    versions             = ch_versions             // channel: [ versions.yml ]
}
