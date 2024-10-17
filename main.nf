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

def primer_set         = ''
def primer_set_version = 0
if (params.platform == 'illumina' && params.protocol == 'amplicon') {
    primer_set         = params.primer_set
    primer_set_version = params.primer_set_version
} else if (params.platform == 'nanopore') {
    primer_set          = params.primer_set
    primer_set_version  = params.primer_set_version
    params.artic_scheme = getGenomeAttribute('scheme', primer_set, primer_set_version)
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
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (params.platform == 'illumina') {
    include { ILLUMINA } from './workflows/illumina'
} else if (params.platform == 'nanopore') {
    include { NANOPORE } from './workflows/nanopore'
}

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

    //
    // WORKFLOW: Run pipeline
    //
    multiqc_report   = Channel.empty()

    if (params.platform == 'illumina') {
        ILLUMINA (
            samplesheet,
            params.fasta,
            params.gff,
            params.primer_bed,
            params.bowtie2_index,
            params.nextclade_dataset,
            params.nextclade_dataset_name,
            params.nextclade_dataset_reference,
            params.nextclade_dataset_tag
        )

        multiqc_report = ILLUMINA.out.multiqc_report

    } else if (params.platform == 'nanopore') {
        NANOPORE (
            samplesheet,
            params.fasta,
            params.gff,
            params.primer_bed,
            params.artic_scheme,
            params.bowtie2_index,
            params.nextclade_dataset,
            params.nextclade_dataset_name,
            params.nextclade_dataset_reference,
            params.nextclade_dataset_tag
        )

        multiqc_report = NANOPORE.out.multiqc_report
    }

    emit:
    multiqc_report // channel: /path/to/multiqc_report.html

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
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def getGenomeAttribute(attribute, primer_set='', primer_set_version=0) {
        def val = ''
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
                            Nextflow.error("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                                " --primer_set_version '${primer_set_version}' not found!\n\n" +
                                " Currently, the available primer set version keys are: ${genome_map.keySet().join(", ")}\n\n" +
                                " Please check:\n" +
                                "   - The value provided to --primer_set_version (currently '${primer_set_version}')\n" +
                                "   - The value provided to --primer_set (currently '${primer_set}')\n" +
                                "   - The value provided to --genome (currently '${params.genome}')\n" +
                                "   - Any custom config files provided to the pipeline.\n\n" + support_link)
                        }
                    } else {
                        Nextflow.error("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                            " --primer_set '${primer_set}' not found!\n\n" +
                            " Currently, the available primer set keys are: ${genome_map.keySet().join(", ")}\n\n" +
                            " Please check:\n" +
                            "   - The value provided to --primer_set (currently '${primer_set}')\n" +
                            "   - The value provided to --genome (currently '${params.genome}')\n" +
                            "   - Any custom config files provided to the pipeline.\n\n" + support_link)
                    }
                } else {
                    Nextflow.error("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
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
