/*
 *This file holds several functions specific to the main.nf workflow in the nf-core/viralrecon pipeline
 */

class WorkflowMain {

    /*
     * Citation string for pipeline
     */
    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
               "* The pipeline\n" + 
               "  https://doi.org/10.5281/zenodo.3901628\n\n" +
               "* The nf-core framework\n" +
               "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
               "* Software dependencies\n" +
               "  https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md"
    }

    /*
     * Print help to screen if required
     */
    public static String help(workflow, params, log) {
        def command = "nextflow run nf-core/viralrecon --input samplesheet.csv --genome 'MN908947.3' -profile docker"
        def help_string = ''
        help_string += NfcoreTemplate.logo(workflow, params.monochrome_logs)
        help_string += NfcoreSchema.paramsHelp(workflow, params, command)
        help_string += '\n' + citation(workflow) + '\n'
        help_string += NfcoreTemplate.dashedLine(params.monochrome_logs)
        return help_string
    }

    /*
     * Print parameter summary log to screen
     */
    public static String paramsSummaryLog(workflow, params, log) {
        def summary_log = ''
        summary_log += NfcoreTemplate.logo(workflow, params.monochrome_logs)
        summary_log += NfcoreSchema.paramsSummaryLog(workflow, params)
        summary_log += '\n' + citation(workflow) + '\n'
        summary_log += NfcoreTemplate.dashedLine(params.monochrome_logs)
        return summary_log
    }

    /*
     * Validate parameters and print summary to screen
     */
    public static void initialise(workflow, params, log) {
        // Print help to screen if required
        if (params.help) {
            log.info help(workflow, params, log)
            System.exit(0)
        }

        // Validate workflow parameters via the JSON schema
        if (params.validate_params) {
            NfcoreSchema.validateParameters(workflow, params, log)
        }

        // Print parameter summary log to screen
        log.info paramsSummaryLog(workflow, params, log)

        // Check that conda channels are set-up correctly
        if (params.enable_conda) {
            Utils.checkCondaChannels(log)
        }

        // Check AWS batch settings
        NfcoreTemplate.awsBatch(workflow, params)

        // Check the hostnames against configured profiles
        NfcoreTemplate.hostName(workflow, params, log)

        // Check sequencing platform
        def platformList = ['illumina', 'nanopore']
        if (!params.public_data_ids) {
            if (!params.platform) {
                log.error "Platform not specified with e.g. '--platform illumina'. Valid options: ${platformList.join(', ')}."
                System.exit(1)
            } else if (!platformList.contains(params.platform)) {
                log.error "Invalid platform option: '${params.platform}'. Valid options: ${platformList.join(', ')}."
                System.exit(1)
            }
        }
    }

    /*
     * Get attribute from genome config file e.g. fasta
     */
    public static String getGenomeAttribute(params, attribute, log, primer_set='', primer_set_version=0) {
        def val = ''
        def support_str = " The default genome config used by the pipeline can be found here:\n" +
                          "   - https://github.com/nf-core/configs/blob/master/conf/pipeline/viralrecon/genomes.config\n\n" +
                          " If you would still like to blame us please come and find us on nf-core Slack:\n" + 
                          "   - https://nf-co.re/viralrecon#contributions-and-support\n" +
                          "============================================================================="
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
                            log.error "=============================================================================\n" +
                                      " --primer_set_version '${primer_set_version}' not found!\n\n" +
                                      " Currently, the available primer set version keys are: ${genome_map.keySet().join(", ")}\n\n" +
                                      " Please check:\n" +
                                      "   - The value provided to --primer_set_version (currently '${primer_set_version}')\n" +
                                      "   - The value provided to --primer_set (currently '${primer_set}')\n" +
                                      "   - The value provided to --genome (currently '${params.genome}')\n" +
                                      "   - Any custom config files provided to the pipeline.\n\n" + support_str
                            System.exit(1)
                        }
                    } else {
                        log.error "=============================================================================\n" +
                                  " --primer_set '${primer_set}' not found!\n\n" +
                                  " Currently, the available primer set keys are: ${genome_map.keySet().join(", ")}\n\n" +
                                  " Please check:\n" +
                                  "   - The value provided to --primer_set (currently '${primer_set}')\n" +
                                  "   - The value provided to --genome (currently '${params.genome}')\n" +
                                  "   - Any custom config files provided to the pipeline.\n\n" + support_str
                        System.exit(1)
                    }
                } else {
                    log.error "=============================================================================\n" +
                              " Genome '${params.genome}' does not contain any primer sets!\n\n" +
                              " Please check:\n" +
                              "   - The value provided to --genome (currently '${params.genome}')\n" +
                              "   - Any custom config files provided to the pipeline.\n\n" + support_str
                    System.exit(1)
                }
            }
            if (genome_map.containsKey(attribute)) {
                val = genome_map[ attribute ]
            }
        }
        return val
    }

    /*
     * Exit pipeline if incorrect --genome key provided
     */
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            log.error "=============================================================================\n" +
                      "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                      "  Currently, the available genome keys are:\n" +
                      "  ${params.genomes.keySet().join(", ")}\n" +
                      "==================================================================================="
            System.exit(1)
        }
    }

    /*
     * Get workflow summary for MultiQC
     */
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    /*
     * Function to check whether primer BED file has the correct suffixes as provided to the pipeline
     */
    public static void checkPrimerSuffixes(primer_bed_file, primer_left_suffix, primer_right_suffix, log) {
        def total = 0
        def left  = 0
        def right = 0
        primer_bed_file.eachLine { line ->
            total += 1
            def name = line.split('\t')[3]
            if (name.contains(primer_left_suffix)) {
                left += 1
            } else if (name.contains(primer_right_suffix)) (
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

    
