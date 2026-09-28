#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/viralrecon
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/viralrecon
    Website: https://nf-co.re/viralrecon
    Slack  : https://nfcore.slack.com/channels/viralrecon
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta         = getGenomeAttribute('fasta')
params.gff           = getGenomeAttribute('gff')
params.bowtie2_index = getGenomeAttribute('bowtie2')

params.nextclade_dataset           = getGenomeAttribute('nextclade_dataset_v3pl')
params.nextclade_dataset_name      = getGenomeAttribute('nextclade_dataset_name')
params.nextclade_dataset_tag       = getGenomeAttribute('nextclade_dataset_tag_v3pl')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { VIRALRECON              } from './workflows/viralrecon'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_viralrecon_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_viralrecon_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main nf-core/viralrecon analysis pipeline depending on type of input
//
workflow NFCORE_VIRALRECON {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    def primer_set         = ''
    def primer_set_version = ''

    // Make sure platform is defined
    if( !params.platform ) {
        error "Parameter --platform is required (illumina / nanopore). Please specify."
    }

    // Check that platform value is valid
    def valid_platforms = ["illumina","nanopore"]
    if( !(params.platform in valid_platforms) ) {
        error "Invalid value for --platform: '${params.platform}'. Allowed values: ${valid_platforms.join(', ')}"
    }

    if (params.platform == 'illumina' && params.trim_primers) {
        primer_set         = params.primer_set
        primer_set_version = params.primer_set_version
    } else if (params.platform == 'nanopore') {
        primer_set          = params.primer_set
        primer_set_version  = params.primer_set_version
    }

    def primer_bed   = params.primer_bed ?: getGenomeAttribute('primer_bed', primer_set, primer_set_version)
    def genome_gff   = (params.gff == false || params.gff == 'false') ? null : params.gff

    //
    // WORKFLOW: Run pipeline
    //
    VIRALRECON (
        samplesheet,
        params.outdir,
        params.fasta,
        genome_gff,
        primer_bed,
        params.bowtie2_index,
        params.nextclade_dataset,
        params.nextclade_dataset_name,
        params.nextclade_dataset_tag
    )

    emit:
    multiqc_report = VIRALRECON.out.multiqc_report // channel: /path/to/multiqc_report.html
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_VIRALRECON (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFCORE_VIRALRECON.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def getGenomeAttribute(attribute, primer_set='', primer_set_version='') {
        def val = null
        def support_link =  " The default genome config used by the pipeline can be found here:\n" +
                            "   - https://github.com/nf-core/configs/blob/master/conf/pipeline/viralrecon/genomes.config\n\n" +
                            " If you would still like to blame us please come and find us on nf-core Slack:\n" +
                            "   - https://nf-co.re/viralrecon#contributions-and-support\n" +
                            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

        if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
            def genome_map = params.genomes[ params.genome ]
            if (primer_set) {
                if (genome_map.containsKey('primer_sets')) {
                    genome_map = genome_map[ 'primer_sets' ]
                    if (genome_map.containsKey(primer_set)) {
                        genome_map = genome_map[ primer_set ]
                        primer_set_version = primer_set_version.toString()
                        if (genome_map.containsKey(primer_set_version)) {
                            genome_map = genome_map[ primer_set_version ]
                        } else {
                            error("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                                " --primer_set_version '${primer_set_version}' not found!\n\n" +
                                " Currently, the available primer set version keys are: ${genome_map.keySet().join(", ")}\n\n" +
                                " Please check:\n" +
                                "   - The value provided to --primer_set_version (currently '${primer_set_version}')\n" +
                                "   - The value provided to --primer_set (currently '${primer_set}')\n" +
                                "   - The value provided to --genome (currently '${params.genome}')\n" +
                                "   - Any custom config files provided to the pipeline.\n\n" + support_link)
                        }
                    } else {
                        error("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                            " --primer_set '${primer_set}' not found!\n\n" +
                            " Currently, the available primer set keys are: ${genome_map.keySet().join(", ")}\n\n" +
                            " Please check:\n" +
                            "   - The value provided to --primer_set (currently '${primer_set}')\n" +
                            "   - The value provided to --genome (currently '${params.genome}')\n" +
                            "   - Any custom config files provided to the pipeline.\n\n" + support_link)
                    }
                } else {
                    error("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                        " Genome '${params.genome}' does not contain any primer sets!\n\n" +
                        " Please check:\n" +
                        "   - The value provided to --genome (currently '${params.genome}')\n" +
                        "   - Any custom config files provided to the pipeline.\n\n" + support_link)
                }
            }
            if (genome_map.containsKey(attribute)) {
                val = genome_map[ attribute ]
            } else if (params.genomes[ params.genome ].containsKey(attribute)) {
                val = params.genomes[ params.genome ][ attribute ]
            }
        }
        return val
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
