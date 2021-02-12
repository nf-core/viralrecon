/*
 * This file holds several functions specific to the pipeline.
 */

class Workflow {

    // Citation string
    private static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
               "* The pipeline\n" + 
               "  https://doi.org/10.5281/zenodo.1400710\n\n" +
               "* The nf-core framework\n" +
               "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
               "* Software dependencies\n" +
               "  https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md"
    }

    static void validate_params(params, log, valid_params) {
        genome_exists(params, log)

        // Generic parameter validation
        if (!params.fasta) { 
            log.error "Genome fasta file not specified!"
            System.exit(0)
        }

        if (!params.skip_kraken2 && !params.kraken2_db) {
            if (!params.kraken2_db_name) { 
                log.error "Please specify a valid name to build Kraken2 database for host e.g. 'human'!"
                System.exit(0)
            }
        }
        
        if (!valid_params['protocols'].contains(params.protocol)) {
            log.error "Invalid protocol option: ${params.protocol}. Valid options: ${valid_params['protocols'].join(', ')}"
            System.exit(0)
        }

        // Variant calling parameter validation
        def callers = params.callers ? params.callers.split(',').collect{ it.trim().toLowerCase() } : []
        if ((valid_params['callers'] + callers).unique().size() != valid_params['callers'].size()) {
            log.error "Invalid variant calller option: ${params.callers}. Valid options: ${valid_params['callers'].join(', ')}"
            System.exit(0)
        }

        if (params.protocol == 'amplicon' && !params.skip_variants && !params.primer_bed) {
            log.error "To perform variant calling in 'amplicon' mode please provide a valid primer BED file!"
            System.exit(0)
        }

        // Assembly parameter validation
        def assemblers = params.assemblers ? params.assemblers.split(',').collect{ it.trim().toLowerCase() } : []
        if ((valid_params['assemblers'] + assemblers).unique().size() != valid_params['assemblers'].size()) {
            log.error "Invalid assembler option: ${params.assemblers}. Valid options: ${valid_params['assemblers'].join(', ')}"
            System.exit(0)
        }

        if (!valid_params['spades_modes'].contains(params.spades_mode)) {
            log.error "Invalid spades mode option: ${params.spades_mode}. Valid options: ${valid_params['spades_modes'].join(', ')}"
            System.exit(0)
        }
    }

    // Print a warning after SRA download has completed
    static void sra_download(log) {
        log.warn "=============================================================================\n" +
                 "  THIS IS AN EXPERIMENTAL FEATURE!\n\n" + 
                 "  Please double-check the samplesheet that has been auto-created using the\n" +
                 "  public database ids provided via the '--public_data_ids' parameter.\n\n" +
                 "  All of the sample metadata obtained from the ENA has been appended\n" +
                 "  as additional columns to help you manually curate the samplesheet before\n" +
                 "  you run the main branch of the pipeline.\n" +
                 "==================================================================================="
    }

    // Exit pipeline if incorrect --genome key provided
    static void genome_exists(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            log.error "=============================================================================\n" +
                      "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                      "  Currently, the available genome keys are:\n" +
                      "  ${params.genomes.keySet().join(", ")}\n" +
                      "==================================================================================="
            System.exit(0)
        }
    }

    // Get attribute from genome config file e.g. fasta
    static String get_genome_attribute(params, attribute) {
        def val = ''
        if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
            if (params.genomes[ params.genome ].containsKey(attribute)) {
                val = params.genomes[ params.genome ][ attribute ]
            }
        }
        return val
    }  

    // Print warning if genome fasta has more than one sequence
    static void is_multifasta(fasta, log) {
        def count = 0
        def line  = null
        fasta.withReader { reader ->
            while (line = reader.readLine()) {
                if (line.contains('>')) {
                    count++
                    if (count > 1) {
                        log.warn "=============================================================================\n" +
                                "  This pipeline does not officially support multi-fasta genome files!\n\n" + 
                                "  The parameters and processes are tailored for viral genome analysis.\n" +
                                "  Please amend the '--fasta' parameter.\n" +
                                "==================================================================================="
                        break
                    }
                }
            }
        }
    }

    // Function that parses and returns the number of mapped reasds from flagstat files
    static ArrayList get_flagstat_mapped_reads(workflow, params, log, flagstat) {
        def mapped_reads = 0
        flagstat.eachLine { line ->
            if (line.contains(' mapped (')) {
                mapped_reads = line.tokenize().first().toInteger()
            }
        }
        
        def pass = false
        def logname = flagstat.getBaseName() - 'flagstat'
        Map colors = Headers.log_colours(params.monochrome_logs)
        if (mapped_reads <= params.min_mapped_reads.toInteger()) {
            log.info "-${colors.purple}[$workflow.manifest.name]${colors.red} [FAIL] Mapped read threshold >= ${params.min_mapped_reads}. IGNORING FOR FURTHER DOWNSTREAM ANALYSIS: ${mapped_reads} - $logname${colors.reset}."
        } else {
            pass = true
            //log.info "-${colors.purple}[$workflow.manifest.name]${colors.green} [PASS] Mapped read threshold >=${params.min_mapped_reads}: ${mapped_reads} - $logname${colors.reset}."
        }
        return [ mapped_reads, pass ]
    }

    static void has_primer_suffixes(primer_bed_file, primer_left_suffix, primer_right_suffix, log) {
        def total = 0
        def left  = 0
        def right = 0
        primer_bed_file.eachLine { 
            line ->
                total += 1
                def name = line.split('\t')[3]
                if (name.endsWith(primer_left_suffix)) {
                    left += 1
                } else if (name.endsWith(primer_right_suffix)) (
                    right += 1
                )
        }
        if (total != (left + right)) {
            log.warn "=============================================================================\n" +
                    "  Please check the name field (column 4) in the file supplied via --primer_bed.\n\n" +
                    "  All of the values in that column do not end with those supplied by:\n" +
                    "      --primer_left_suffix : $primer_left_suffix\n" +
                    "      --primer_right_suffix: $primer_right_suffix\n\n" +
                    "  This information is required to collapse the primer intervals into amplicons\n" + 
                    "  for the coverage plots generated by the pipeline.\n" +
                    "==================================================================================="
        }
    }
}
