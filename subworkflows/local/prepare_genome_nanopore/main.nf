//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA      } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFF        } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_PRIMER_BED } from '../../../modules/nf-core/gunzip/main'
include { UNTAR                       } from '../../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_KRAKEN2_DB   } from '../../../modules/nf-core/untar/main'
include { SAMTOOLS_FAIDX              } from '../../../modules/nf-core/samtools/faidx/main'
include { NEXTCLADE_DATASETGET        } from '../../../modules/nf-core/nextclade/datasetget/main'
include { COLLAPSE_PRIMERS            } from '../../../modules/local/collapse_primers'
include { KRAKEN2_BUILD               } from '../../../modules/local/kraken2/build'
include { SNPEFF_BUILD                } from '../../../modules/local/snpeff/build'

workflow PREPARE_GENOME_NANOPORE {

    take:
    fasta
    gff
    primer_bed
    nextclade_dataset
    nextclade_dataset_name
    nextclade_dataset_tag

    main:

    ch_versions = channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    if (fasta.endsWith('.gz')) {
        GUNZIP_FASTA (
            [ [:], fasta ]
        )
        ch_fasta    = GUNZIP_FASTA.out.gunzip.map { _meta, gunzip -> gunzip }
    } else {
        ch_fasta = channel.value(file(fasta))
    }

    //
    // Uncompress GFF annotation file
    //
    ch_gff = channel.empty()
    if (gff) {
        if (gff.endsWith('.gz')) {
            GUNZIP_GFF (
                [ [:], gff ]
            )
            ch_gff      = GUNZIP_GFF.out.gunzip.map { _meta, gunzip -> gunzip }
        } else {
            ch_gff = channel.value(file(gff))
        }
    }

    //
    // Create chromosome sizes file
    //
    SAMTOOLS_FAIDX (
        ch_fasta.map { fasta_file ->
            [ [:], fasta_file, [] ] },
        true
    )
    ch_fai         = SAMTOOLS_FAIDX.out.fai.map { _meta, fai -> fai  }
    ch_chrom_sizes = SAMTOOLS_FAIDX.out.sizes.map { _meta, sizes -> sizes }

    //
    // Prepare reference files required for variant calling
    //
    ch_kraken2_db = channel.empty()
    if (!params.skip_kraken2) {
        if (params.kraken2_db) {
            if (params.kraken2_db.endsWith('.tar.gz')) {
                UNTAR_KRAKEN2_DB (
                    [ [:], params.kraken2_db ]
                )
                ch_kraken2_db = UNTAR_KRAKEN2_DB.out.untar.map { _meta, kraken2_db -> kraken2_db }
            } else {
                ch_kraken2_db = channel.value(file(params.kraken2_db))
            }
        } else {
            KRAKEN2_BUILD (
                params.kraken2_db_name
            )
            ch_kraken2_db = KRAKEN2_BUILD.out.db.first()
        }
    }

    //
    // Uncompress primer BED file
    //
    ch_primer_bed = channel.empty()
    if (primer_bed) {
        if (primer_bed.endsWith('.gz')) {
            GUNZIP_PRIMER_BED (
                [ [:], primer_bed ]
            )
            ch_primer_bed = GUNZIP_PRIMER_BED.out.gunzip.map { _meta, gunzip -> gunzip }
        } else {
            ch_primer_bed = channel.value(file(primer_bed))
        }
    }

    //
    // Generate collapsed BED file
    //
    ch_primer_collapsed_bed = channel.empty()
    if (!params.skip_mosdepth) {
        COLLAPSE_PRIMERS (
            ch_primer_bed,
            params.primer_left_suffix,
            params.primer_right_suffix
        )
        ch_primer_collapsed_bed = COLLAPSE_PRIMERS.out.bed
    }

    //
    // Prepare Nextclade dataset
    //
    ch_nextclade_db = channel.empty()
    ch_versions = channel.empty()
    if (!params.skip_consensus && !params.skip_nextclade) {
        if (nextclade_dataset) {
            if (nextclade_dataset.endsWith('.tar.gz')) {
                UNTAR (
                    [ [:], nextclade_dataset ]
                )
                ch_nextclade_db = UNTAR.out.untar.map { _meta, untar -> untar }
            } else {
                ch_nextclade_db = channel.value(file(nextclade_dataset))
            }
        } else if (nextclade_dataset_name) {
            NEXTCLADE_DATASETGET (
                nextclade_dataset_name,
                nextclade_dataset_tag
            )
            ch_nextclade_db = NEXTCLADE_DATASETGET.out.dataset
            ch_versions     = ch_versions.mix(NEXTCLADE_DATASETGET.out.versions)
        }
    }

    //
    // Make snpEff database
    //
    ch_snpeff_db     = channel.empty()
    ch_snpeff_config = channel.empty()
    if (!params.skip_snpeff) {
        SNPEFF_BUILD (
            ch_fasta,
            ch_gff
        )
        ch_snpeff_db     = SNPEFF_BUILD.out.db
        ch_snpeff_config = SNPEFF_BUILD.out.config
    }

    //
    // Materialize reference channels so they can be reused by multiple consumers
    //
    def ch_reference_fasta                = ch_fasta.collect(flat: false).map { files -> files[0] }
    def ch_reference_fai                  = ch_fai.collect(flat: false).map { files -> files[0] }
    def ch_reference_gff                  = gff ? ch_gff.collect(flat: false).map { files -> files[0] } : []
    def ch_reference_primer_bed           = ch_primer_bed.collect(flat: false).map { files -> files[0] }
    def ch_reference_primer_collapsed_bed = !params.skip_mosdepth ? ch_primer_collapsed_bed.collect(flat: false).map { files -> files[0] } : []
    def ch_reference_nextclade_db         = !params.skip_nextclade ? ch_nextclade_db.collect(flat: false).map { dirs -> dirs[0] } : []
    def ch_reference_kraken2_db           = !params.skip_kraken2 ? ch_kraken2_db.collect(flat: false).map { dirs -> dirs[0] } : []
    def ch_reference_snpeff_db            = (gff && !params.skip_snpeff) ? ch_snpeff_db.collect(flat: false).map { dirs -> dirs[0] } : []
    def ch_reference_snpeff_config        = (gff && !params.skip_snpeff) ? ch_snpeff_config.collect(flat: false).map { files -> files[0] } : []

    emit:
    fasta                = ch_reference_fasta                // path: genome.fasta
    gff                  = ch_reference_gff                  // path: genome.gff
    fai                  = ch_reference_fai                  // path: genome.fai
    chrom_sizes          = ch_chrom_sizes                    // path: genome.sizes
    primer_bed           = ch_reference_primer_bed           // path: primer.bed
    primer_collapsed_bed = ch_reference_primer_collapsed_bed // path: primer.collapsed.bed
    nextclade_db         = ch_reference_nextclade_db         // path: nextclade_db
    kraken2_db           = ch_reference_kraken2_db           // path: kraken2_db/
    snpeff_db            = ch_reference_snpeff_db            // path: snpeff_db
    snpeff_config        = ch_reference_snpeff_config        // path: snpeff.config

    versions             = ch_versions                       // channel: versions.yml
}
