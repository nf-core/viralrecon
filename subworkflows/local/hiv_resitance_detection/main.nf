//
// HIV resistance detecton
//

include { SIERRALOCAL                                        } from '../../../modules/local/sierralocal'
include { LIFTOFF                                            } from '../../../modules/nf-core/liftoff'
include { GFF2JSON                                           } from '../../../modules/local/gff2json'
include { BAM2CODFREQ                                        } from '../../../modules/local/bam2codfreq'
include { RESISTANCE_TABLES                                  } from '../../../modules/local/resistance_tables'
include { RESISTANCE_REPORT                                  } from '../../../modules/local/resistance_report'
include { ADDITIONAL_ANNOTATION as HIV_RESISTANCE_ANNOTATION } from '../additional_annotation'
include { LIFTOFF as CONSENSUS_LIFTOFF                       } from '../../../modules/nf-core/liftoff'

workflow HIV_RESISTANCE {
    take:
    consensus        // channel: [ val(meta), [ consensus ] ]
    bam              // channel: [ val(meta), [ bam ], [bai] ]
    fasta            // path   : genome.fasta
    gff              // path   : genome.gff
    vcf              // channel: [ val(meta), [ vcf ] ]
    tbi              // channel: [ val(meta), [ tbi ] ]
    pangolin         // channel: [ val(meta), [ csv ] ]
    nextclade_report // channel: [ val(meta), [ csv ] ]

    main:

    //
    // HIV resistance detection
    //

    SIERRALOCAL (
        consensus,
        params.hivdb_xml    ? file(params.hivdb_xml, checkIfExists: true)    : [],
        params.apobec_drm   ? file(params.apobec_drm, checkIfExists: true)   : [],
        params.apobec_csv   ? file(params.apobec_csv, checkIfExists: true)   : [],
        params.unusual_csv  ? file(params.unusual_csv, checkIfExists: true)  : [],
        params.sdrms_csv    ? file(params.sdrms_csv, checkIfExists: true)    : [],
        params.mutation_csv ? file(params.mutation_csv, checkIfExists: true) : []
    )

    codfreq_sequence   = file("$projectDir/assets/codfreq.fasta", checkIfExists: true)
    codfreq_annotation = file("$projectDir/assets/codfreq.gff", checkIfExists: true)

    if (params.genome != 'codfreq') {

        // Not setted as param because we can't changed these files for codfreq and sierralocal to work toghether.
        // Using the ones in assets instead of the ones in test-datasets to avoid download issues when running offline
        LIFTOFF (
            fasta.map { f -> [ [id: f.baseName], f ] },
            codfreq_sequence,
            codfreq_annotation,
            []
        )

        GFF2JSON (
            fasta,
            LIFTOFF.out.gff3
        )

        HIV_RESISTANCE_ANNOTATION (
            vcf,
            tbi,
            fasta,
            LIFTOFF.out.gff3.map { _meta, annotation -> annotation },
            pangolin
        )

    } else {
        GFF2JSON (
            fasta,
            gff.map { f -> [ [id: f.baseName], f ] }
        )
    }

    BAM2CODFREQ (
        bam,
        GFF2JSON.out.profile_json
    )

    CONSENSUS_LIFTOFF (
        consensus,
        codfreq_sequence,
        codfreq_annotation,
        []
    )

    RESISTANCE_TABLES(
        SIERRALOCAL.out.json.join(BAM2CODFREQ.out.codfreq, by: [0])
    )

    SIERRALOCAL.out.json.dump(tag:"SIERRALOCAL.out.json")
    RESISTANCE_TABLES.out.mutation_csv.dump(tag:"RESISTANCE_TABLES.out.mutation_csv")
    RESISTANCE_TABLES.out.resistance_csv.dump(tag:"RESISTANCE_TABLES.out.resistance_csv")
    nextclade_report.dump(tag:"nextclade_report")
    consensus.dump(tag:"consensus")
    CONSENSUS_LIFTOFF.out.gff3.dump(tag:"CONSENSUS_LIFTOFF.out.gff3")

    ch_resistance_report_input = SIERRALOCAL.out.json
        .join(RESISTANCE_TABLES.out.mutation_csv, by: [0])
        .join(RESISTANCE_TABLES.out.resistance_csv, by: [0])
        .join(nextclade_report, by: [0])
        .join(consensus, by: [0])
        .join(CONSENSUS_LIFTOFF.out.gff3, by: [0])

    RESISTANCE_REPORT(ch_resistance_report_input)

    emit:
    sierralocal_results  = SIERRALOCAL.out.json                      // channel: [ val(meta), [ json ] ]
    bam2codfreq_results  = BAM2CODFREQ.out.codfreq                   // channel: [ val(meta), [ codfreq ] ]
    mutation_table       = RESISTANCE_TABLES.out.mutation_csv        // channel: [ val(meta), [ mutation_csv ] ]
    mutation_table_short = RESISTANCE_TABLES.out.mutation_short_csv  // channel: [ val(meta), [ mutation_csv ] ]
    resistance_table     = RESISTANCE_TABLES.out.resistance_csv      // channel: [ val(meta), [ resistance_csv ] ]
    resistance_report    = RESISTANCE_REPORT.out.html
}
