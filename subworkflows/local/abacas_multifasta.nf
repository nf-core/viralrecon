//
// Process multiple multifasta files with abacas
//
include { ABACAS        } from '../../modules/nf-core/abacas/main'

workflow ABACAS_MULTI {
    take:
    scaffold   // channel: [ val(meta), path(scaffold) ]
    multifasta // channel: /path/to/genome.fasta
    assembler  // string: assembler name

    main:

    ch_versions = Channel.empty()

    //
    // Split multifasta file into individual fasta files
    //
    multifasta
        .splitFasta( by: 1, file: true )
        .set { ch_fasta }

    //
    // Run abacas on each fasta file
    //
    ch_abacas = Channel.empty()

    scaffold
        .combine (ch_fasta)
        .set { ch_scaffold_fasta }

    ABACAS (
        ch_scaffold_fasta.map { meta, scaffold, fasta -> tuple( meta, scaffold ) },
        ch_scaffold_fasta.map { meta, scaffold, fasta -> fasta }
    )
    ch_abacas = ABACAS.out.results
    ch_versions = ch_versions.mix(ABACAS.out.versions)

    //
    // Concatenate abacas results
    //

    ch_abacas
        .map { meta, files -> tuple(meta, files[3]) }
        .multiMap{ meta, fasta ->
            metadata: [meta.id, meta]
            fasta: [meta.id, fasta]
        }
        .set { ch_abacas_split }

        ch_abacas_split.fasta
        .collectFile (storeDir: "${params.outdir}/assembly/${assembler}/abacas_multi") { id, fasta ->
            ["${id}.fa",fasta]
        }
        .map { file -> [file.simpleName, file] }
        .join(ch_abacas_split.metadata)
        .map { id, fasta, meta -> tuple(meta, fasta) }
        .set { ch_abacas_results }

    emit:
    abacas_results     = ch_abacas_results     // channel: [ val(meta), path('*.abacas*') ]
    versions           = ch_versions           // channel: [ versions.yml ]
}
