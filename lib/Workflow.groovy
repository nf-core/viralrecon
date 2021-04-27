/*
 * This file holds several functions specific to the pipeline.
 */

import groovy.json.JsonSlurper

class Workflow {

    // Citation string
    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
               "* The pipeline\n" + 
               "  https://doi.org/10.5281/zenodo.3901628\n\n" +
               "* The nf-core framework\n" +
               "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
               "* Software dependencies\n" +
               "  https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md"
    }

    public static void validateMainParams(workflow, params, json_schema, log) {
        // Validate workflow parameters via the JSON schema
        if (params.validate_params) {
            NfcoreSchema.validateParameters(params, json_schema, log)
        }

        // Check that conda channels are set-up correctly
        if (params.enable_conda) {
            Checks.checkCondaChannels(log)
        }

        // Check AWS batch settings
        Checks.awsBatch(workflow, params)

        // Check the hostnames against configured profiles
        Checks.hostName(workflow, params, log)

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

    public static void validateIlluminaParams(params, log, valid_params) {
        genomeExists(params, log)

        // Generic parameter validation
        if (!valid_params['protocols'].contains(params.protocol)) {
            log.error "Invalid option: '${params.protocol}'. Valid options for '--protocol': ${valid_params['protocols'].join(', ')}."
            System.exit(1)
        }

        if (!params.fasta) { 
            log.error "Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file."
            System.exit(1)
        }

        if (!params.skip_kraken2 && !params.kraken2_db) {
            if (!params.kraken2_db_name) { 
                log.error "Please specify a valid name to build Kraken2 database for host e.g. '--kraken2_db_name human'."
                System.exit(1)
            }
        }
        
        // Variant calling parameter validation
        def callers = params.callers ? params.callers.split(',').collect{ it.trim().toLowerCase() } : []
        if ((valid_params['callers'] + callers).unique().size() != valid_params['callers'].size()) {
            log.error "Invalid option: ${params.callers}. Valid options for '--callers': ${valid_params['callers'].join(', ')}."
            System.exit(1)
        }

        if (params.protocol == 'amplicon' && !params.skip_variants && !params.primer_bed) {
            log.error "To perform variant calling in amplicon mode please provide a valid primer BED file e.g. '--primer_bed primers.bed'."
            System.exit(1)
        }

        // Assembly parameter validation
        def assemblers = params.assemblers ? params.assemblers.split(',').collect{ it.trim().toLowerCase() } : []
        if ((valid_params['assemblers'] + assemblers).unique().size() != valid_params['assemblers'].size()) {
            log.error "Invalid option: ${params.assemblers}. Valid options for '--assemblers': ${valid_params['assemblers'].join(', ')}."
            System.exit(1)
        }

        if (!valid_params['spades_modes'].contains(params.spades_mode)) {
            log.error "Invalid option: ${params.spades_mode}. Valid options for '--spades_modes': ${valid_params['spades_modes'].join(', ')}."
            System.exit(1)
        }
    }

    public static void validateNanoporeParams(params, log, valid_params) {
        genomeExists(params, log)

        // Generic parameter validation
        if (!params.fasta) { 
            log.error "Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file."
            System.exit(1)
        }

        if (!params.primer_bed) {
            log.error "Primer BED file not specified with e.g. '--primer_bed primers.bed' or via a detectable config file."
            System.exit(1)
        }

        if (!params.artic_scheme) {
            log.error "ARTIC scheme not specified with e.g. --artic_scheme 'nCoV-2019' or via a detectable config file."
            System.exit(1)
        }

        if (!valid_params['artic_minion_caller'].contains(params.artic_minion_caller)) {
            log.error "Invalid option: ${params.artic_minion_caller}. Valid options for '--artic_minion_caller': ${valid_params['artic_minion_caller'].join(', ')}."
            System.exit(1)
        }

        if (!valid_params['artic_minion_aligner'].contains(params.artic_minion_aligner)) {
            log.error "Invalid option: ${params.artic_minion_aligner}. Valid options for '--artic_minion_aligner': ${valid_params['artic_minion_aligner'].join(', ')}."
            System.exit(1)
        }

        if (!params.fastq_dir) {
            log.error "Please specify a valid folder containing ONT basecalled fastq files generated by guppy_barcoder or guppy_basecaller e.g. '--fastq_dir ./20191023_1522_MC-110615_0_FAO93606_12bf9b4f/fastq_pass/"
            System.exit(1)
        }

        if (params.artic_minion_caller == 'nanopolish') {
            if (!params.fast5_dir) {
                log.error "Please specify a valid folder containing ONT fast5 files e.g. '--fast5_dir ./20191023_1522_MC-110615_0_FAO93606_12bf9b4f/fast5_pass/"
                System.exit(1)
            }
            if (!params.sequencing_summary) {
                log.error "Please specify a valid ONT sequencing summary file e.g. '--sequencing_summary ./20191023_1522_MC-110615_0_FAO93606_12bf9b4f/sequencing_summary.txt"
                System.exit(1)
            }
        }

        if (params.artic_minion_caller == 'medaka') {
            if (!params.artic_minion_medaka_model) {
                log.error "Please specify the '--artic_minion_medaka_model' parameter too if using the '--artic_minion_caller medaka' workflow.\nSee https://github.com/nanoporetech/medaka"
                System.exit(1)
            }
        }
    }

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

    // Print warning if genome fasta has more than one sequence
    public static void isMultiFasta(fasta_file, log) {
        def count = 0
        def line  = null
        fasta_file.withReader { reader ->
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
    public static ArrayList getFlagstatMappedReads(flagstat_file, params) {
        def mapped_reads = 0
        flagstat_file.eachLine { line ->
            if (line.contains(' mapped (')) {
                mapped_reads = line.tokenize().first().toInteger()
            }
        }
        
        def pass = false
        def logname = flagstat_file.getBaseName() - 'flagstat'
        if (mapped_reads > params.min_mapped_reads.toInteger()) {
            pass = true
        }
        return [ mapped_reads, pass ]
    }

    // Function to check whether primer BED file has the correct suffixes as provided to the pipeline
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

    // Check if the primer BED file supplied to the pipeline is from the SWIFT/SNAP protocol
    public static void checkIfSwiftProtocol(primer_bed_file, log) {
        def count = 0
        def line  = null
        primer_bed_file.withReader { reader ->
            while (line = reader.readLine()) {
                def name = line.split('\t')[3]
                if (name.contains('covid19genome')) {
                    count++
                    if (count > 1) {
                        log.warn "=============================================================================\n" +
                                 "  Found 'covid19genome' in the name field of the primer BED file!\n" + 
                                 "  This suggests that you have used the SWIFT/SNAP protocol to prep your samples.\n" + 
                                 "  If so, please set '--ivar_trim_offset 5' as suggested in the issue below:\n" +
                                 "  https://github.com/nf-core/viralrecon/issues/170\n" +
                                 "==================================================================================="
                        break
                    }
                }
            }
        }
    }

    // Function that parses fastp json output file to get total number of reads after trimming
    public static Integer getFastpReadsAfterFiltering(json_file) {
        def Map json = (Map) new JsonSlurper().parseText(json_file.text).get('summary')
        return json['after_filtering']['total_reads'].toInteger()
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

    // Print a warning after SRA download has completed
    public static void sraDownload(log) {
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
    private static void genomeExists(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            log.error "=============================================================================\n" +
                      "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                      "  Currently, the available genome keys are:\n" +
                      "  ${params.genomes.keySet().join(", ")}\n" +
                      "==================================================================================="
            System.exit(1)
        }
    }
}
