//
// This file holds several functions specific to the workflow/illumina.nf in the nf-core/viralrecon pipeline
//

import groovy.json.JsonSlurper

class WorkflowIllumina {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log, valid_params) {
        WorkflowCommons.genomeExistsError(params, log)

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

    //
    // Print warning if genome fasta has more than one sequence
    //
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

    //
    // Function that parses and returns the number of mapped reasds from flagstat files
    //
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

    //
    // Check if the primer BED file supplied to the pipeline is from the SWIFT/SNAP protocol
    //
    public static void checkIfSwiftProtocol(primer_bed_file, name_prefix, log) {
        def count = 0
        def line  = null
        primer_bed_file.withReader { reader ->
            while (line = reader.readLine()) {
                def name = line.split('\t')[3]
                if (name.contains(name_prefix)) {
                    count++
                    if (count > 1) {
                        log.warn "=============================================================================\n" +
                            "  Found '${name_prefix}' in the name field of the primer BED file!\n" +
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

    //
    // Function that parses fastp json output file to get total number of reads after trimming
    //
    public static Integer getFastpReadsAfterFiltering(json_file) {
        def Map json = (Map) new JsonSlurper().parseText(json_file.text).get('summary')
        return json['after_filtering']['total_reads'].toInteger()
    }

    //
    // Function that parses fastp json output file to get total number of reads before trimming
    //
    public static Integer getFastpReadsBeforeFiltering(json_file) {
        def Map json = (Map) new JsonSlurper().parseText(json_file.text).get('summary')
        return json['before_filtering']['total_reads'].toInteger()
    }
}
