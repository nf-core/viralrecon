//
// This file holds several functions specific to the main.nf workflow in the nf-core/viralrecon pipeline
//

import nextflow.Nextflow

class WorkflowMain {

    //
    // Citation string for pipeline
    //
    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
            "* The pipeline\n" +
            "  https://doi.org/10.5281/zenodo.3901628\n\n" +
            "* The nf-core framework\n" +
            "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
            "* Software dependencies\n" +
            "  https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md"
    }

    //
    // Validate parameters and print summary to screen
    //
    public static void initialise(workflow, params, log, args) {

        // Print workflow version and exit on --version
        if (params.version) {
            String workflow_version = NfcoreTemplate.version(workflow)
            log.info "${workflow.manifest.name} ${workflow_version}"
            System.exit(0)
        }

        // Warn about using custom configs to provide pipeline parameters
        NfcoreTemplate.warnParamsProvidedInConfig(workflow, log)

        // Check that a -profile or Nextflow config has been provided to run the pipeline
        NfcoreTemplate.checkConfigProvided(workflow, log)
        // Check that the profile doesn't contain spaces and doesn't end with a trailing comma
        checkProfile(workflow.profile, args, log)

        // Check that conda channels are set-up correctly
        if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
            Utils.checkCondaChannels(log)
        }

        // Check AWS batch settings
        NfcoreTemplate.awsBatch(workflow, params)

        // Check sequencing platform
        def platformList = ['illumina', 'nanopore']
        if (!params.platform) {
            Nextflow.error("Platform not specified with e.g. '--platform illumina'. Valid options: ${platformList.join(', ')}.")
        } else if (!platformList.contains(params.platform)) {
            Nextflow.error("Invalid platform option: '${params.platform}'. Valid options: ${platformList.join(', ')}.")
        }

        // Check Nextclade dataset parameters
        if (!params.skip_consensus && !params.skip_nextclade) {
            if (!params.nextclade_dataset && !params.nextclade_dataset_name) {
                Nextflow.error("Nextclade dataset not specified with '--nextclade_dataset' or '--nextclade_dataset_name'. A list of available datasets can be obtained using the Nextclade 'nextclade dataset list' command.")
            }
        }
    }
    //
    // Get attribute from genome config file e.g. fasta
    //
    public static String getGenomeAttribute(params, attribute, log, primer_set='', primer_set_version=0) {
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

    //
    // Exit pipeline if --profile contains spaces
    //
    private static void checkProfile(profile, args, log) {
        if (profile.endsWith(',')) {
            Nextflow.error "Profile cannot end with a trailing comma. Please remove the comma from the end of the profile string.\nHint: A common mistake is to provide multiple values to `-profile` separated by spaces. Please use commas to separate profiles instead,e.g., `-profile docker,test`."
        }
        if (args[0]) {
            log.warn "nf-core pipelines do not accept positional arguments. The positional argument `${args[0]}` has been detected.\n      Hint: A common mistake is to provide multiple values to `-profile` separated by spaces. Please use commas to separate profiles instead,e.g., `-profile docker,test`."
        }
    }
}
