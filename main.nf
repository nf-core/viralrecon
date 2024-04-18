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

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.platform == 'illumina') {
    include { ILLUMINA as VIRALRECON } from './workflows/illumina'
} else if (params.platform == 'nanopore') {
    include { NANOPORE as VIRALRECON } from './workflows/nanopore'
}

include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_viralrecon_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_viralrecon_pipeline'

include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_viralrecon_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def primer_set         = ''
def primer_set_version = 0
if (params.platform == 'illumina' && params.protocol == 'amplicon') {
    primer_set         = params.primer_set
    primer_set_version = params.primer_set_version
} else if (params.platform == 'nanopore') {
    primer_set          = params.primer_set
    primer_set_version  = params.primer_set_version
    params.artic_scheme = getGenomeAttribute('scheme')
}

params.fasta         = getGenomeAttribute('fasta')
params.gff           = getGenomeAttribute('gff')
params.bowtie2_index = getGenomeAttribute('bowtie2')
params.primer_bed    = getGenomeAttribute('primer_bed', primer_set, primer_set_version)

params.nextclade_dataset           = getGenomeAttribute('nextclade_dataset')
params.nextclade_dataset_name      = getGenomeAttribute('nextclade_dataset_name')
params.nextclade_dataset_reference = getGenomeAttribute('nextclade_dataset_reference')
params.nextclade_dataset_tag       = getGenomeAttribute('nextclade_dataset_tag')

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

    //
    // WORKFLOW: Run pipeline
    //
    VIRALRECON (
        samplesheet
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
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
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
        params.hook_url,
        NFCORE_VIRALRECON.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
