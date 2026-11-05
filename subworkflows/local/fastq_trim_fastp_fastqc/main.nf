//
// Read QC and trimming
//

include { FASTQC as FASTQC_RAW  } from '../../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIM } from '../../../modules/nf-core/fastqc/main'
include { FASTP                 } from '../../../modules/nf-core/fastp/main'

//
// Function that parses fastp json output file to get total number of reads after trimming
//
def getFastpReadsAfterFiltering(json_file) {
    def Map json = new groovy.json.JsonSlurper().parseText(json_file.text).get('summary')
    return json['after_filtering']['total_reads'].toInteger()
}

workflow FASTQ_TRIM_FASTP_FASTQC {
    take:
    reads                // channel: [ val(meta), [ reads ] ]
    adapter_fasta        //    file: adapter.fasta
    discard_trimmed_pass //   value: boolean
    save_trimmed_fail    //   value: boolean
    save_merged          //   value: boolean

    main:

    fastqc_raw_html = channel.empty()
    fastqc_raw_zip  = channel.empty()
    if (!params.skip_fastqc) {
        FASTQC_RAW (
            reads
        )
        fastqc_raw_html = FASTQC_RAW.out.html
        fastqc_raw_zip  = FASTQC_RAW.out.zip
    }

    trim_reads        = reads
    trim_json         = channel.empty()
    trim_html         = channel.empty()
    trim_log          = channel.empty()
    trim_reads_fail   = channel.empty()
    trim_reads_merged = channel.empty()
    fastqc_trim_html  = channel.empty()
    fastqc_trim_zip   = channel.empty()
    if (!params.skip_fastp) {
        FASTP (
            reads.map { meta, reads_ -> tuple(meta, reads_, adapter_fasta) },
            discard_trimmed_pass,
            save_trimmed_fail,
            save_merged
        )
        trim_reads        = FASTP.out.reads
        trim_json         = FASTP.out.json
        trim_html         = FASTP.out.html
        trim_log          = FASTP.out.log
        trim_reads_fail   = FASTP.out.reads_fail
        trim_reads_merged = FASTP.out.reads_merged

        //
        // Filter empty FastQ files after adapter trimming so FastQC doesn't fail
        //
        trim_reads
            .join(trim_json)
            .map {
                meta, reads_, json ->
                    if (getFastpReadsAfterFiltering(json) > 0) {
                        [ meta, reads_ ]
                    }
            }
            .set { trim_reads }

        if (!params.skip_fastqc) {
            FASTQC_TRIM (
                trim_reads
            )
            fastqc_trim_html = FASTQC_TRIM.out.html
            fastqc_trim_zip  = FASTQC_TRIM.out.zip
        }
    }

    emit:
    reads    = trim_reads  // channel: [ val(meta), [ reads ] ]
    trim_json              // channel: [ val(meta), [ json ] ]
    trim_html              // channel: [ val(meta), [ html ] ]
    trim_log               // channel: [ val(meta), [ log ] ]
    trim_reads_fail        // channel: [ val(meta), [ fastq.gz ] ]
    trim_reads_merged      // channel: [ val(meta), [ fastq.gz ] ]

    fastqc_raw_html        // channel: [ val(meta), [ html ] ]
    fastqc_raw_zip         // channel: [ val(meta), [ zip ] ]
    fastqc_trim_html       // channel: [ val(meta), [ html ] ]
    fastqc_trim_zip        // channel: [ val(meta), [ zip ] ]
}
