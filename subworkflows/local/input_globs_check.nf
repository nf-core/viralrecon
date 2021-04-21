params.options = [:]

workflow INPUT_GLOBS_CHECK {
    take:
    input // file  : /path/to/samplesheet.csv
    platform    // string: sequencing platform. Accepted values: 'illumina', 'nanopore'
    
    main:

    input
        .map { create_fastq_channels(it) }
        .set { reads }

    emit:
    reads // channel: [ val(meta), [ reads ] ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channels(ArrayList tuple) {
    def meta = [:]
    meta.id           = tuple[0]
    meta.single_end   = params.single_end

    array = [ meta, tuple[1] ]
    return array    
}