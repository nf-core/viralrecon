#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/viralrecon
========================================================================================
 nf-core/viralrecon Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/viralrecon
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

////////////////////////////////////////////////////
/* --               PRINT HELP                 -- */
////////////////////////////////////////////////////

log.info Utils.logo(workflow, params.monochrome_logs)

def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run nf-core/viralrecon --input samplesheet.csv --genome 'MN908947.3' -profile docker"
    log.info NfcoreSchema.paramsHelp(workflow, params, json_schema, command)
    log.info Workflow.citation(workflow)
    log.info Utils.dashedLine(params.monochrome_logs)
    exit 0
}

////////////////////////////////////////////////////
/* --        GENOME PARAMETER VALUES           -- */
////////////////////////////////////////////////////

def primer_set         = ''
def primer_set_version = 0
if (!params.public_data_ids && params.platform == 'illumina' && params.protocol == 'amplicon') {
    primer_set         = params.primer_set
    primer_set_version = params.primer_set_version
} else if (!params.public_data_ids && params.platform == 'nanopore') {
    primer_set          = 'artic'
    primer_set_version  = params.primer_set_version
    params.artic_scheme = Workflow.getGenomeAttribute(params, 'scheme', log, primer_set, primer_set_version)
}

params.fasta         = Workflow.getGenomeAttribute(params, 'fasta'     , log, primer_set, primer_set_version)
params.gff           = Workflow.getGenomeAttribute(params, 'gff'       , log, primer_set, primer_set_version)
params.bowtie2_index = Workflow.getGenomeAttribute(params, 'bowtie2'   , log, primer_set, primer_set_version)
params.primer_bed    = Workflow.getGenomeAttribute(params, 'primer_bed', log, primer_set, primer_set_version)

////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params, json_schema)
log.info NfcoreSchema.paramsSummaryLog(workflow, params, json_schema)
log.info Workflow.citation(workflow)
log.info Utils.dashedLine(params.monochrome_logs)

////////////////////////////////////////////////////
/* --          PARAMETER CHECKS                -- */
////////////////////////////////////////////////////

Workflow.validateMainParams(workflow, params, log)

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow {
    /*
     * WORKFLOW: Get SRA run information for public database ids, download and md5sum check FastQ files, auto-create samplesheet
     */
    if (params.public_data_ids) {
        include { SRA_DOWNLOAD } from './workflows/sra_download' addParams( summary_params: summary_params )
        SRA_DOWNLOAD ()
    
    /*
     * WORKFLOW: Variant and de novo assembly analysis for Illumina data
     */
    } else if (params.platform == 'illumina') {
        include { ILLUMINA } from './workflows/illumina' addParams( summary_params: summary_params )
        ILLUMINA ()

    /*
     * WORKFLOW: Variant analysis for Nanopore data
     */ 
    } else if (params.platform == 'nanopore') {
        include { NANOPORE } from './workflows/nanopore' addParams( summary_params: summary_params )
        NANOPORE ()
    }
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
