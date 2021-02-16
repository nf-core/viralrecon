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

def json_schema = "$projectDir/nextflow_schema.json"
if (params.help) {
    def command = "nextflow run nf-core/viralrecon --input samplesheet.csv --genome 'MN908947.3' -profile docker"
    log.info Schema.params_help(workflow, params, json_schema, command)
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
    params.artic_scheme = Workflow.get_genome_attribute(params, 'scheme', log, primer_set, primer_set_version)
}

params.fasta         = Workflow.get_genome_attribute(params, 'fasta'     , log, primer_set, primer_set_version)
params.gff           = Workflow.get_genome_attribute(params, 'gff'       , log, primer_set, primer_set_version)
params.bowtie2_index = Workflow.get_genome_attribute(params, 'bowtie2'   , log, primer_set, primer_set_version)
params.primer_bed    = Workflow.get_genome_attribute(params, 'primer_bed', log, primer_set, primer_set_version)

////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////

def summary_params = Schema.params_summary_map(workflow, params, json_schema)
log.info Schema.params_summary_log(workflow, params, json_schema)

////////////////////////////////////////////////////
/* --          PARAMETER CHECKS                -- */
////////////////////////////////////////////////////

Workflow.validate_main_params(workflow, params, log)

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
