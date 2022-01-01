//
// Read QC and trimming
//

params.fastqc_raw_options  = [:]
params.fastqc_trim_options = [:]
params.fastp_options       = [:]

include { FASTQC as FASTQC_RAW  } from '../../modules/nf-core/modules/fastqc/main' addParams( options: params.fastqc_raw_options  )
include { FASTQC as FASTQC_TRIM } from '../../modules/nf-core/modules/fastqc/main' addParams( options: params.fastqc_trim_options )
include { FASTP                 } from '../../modules/nf-core/modules/fastp/main'  addParams( options: params.fastp_options       )

workflow FASTQC_FASTP {
    take:
    reads // channel: [ val(meta), [ reads ] ]

    main:
    fastqc_raw_html = Channel.empty()
    fastqc_raw_zip  = Channel.empty()
    fastqc_version  = Channel.empty()
    if (!params.skip_fastqc) {
        FASTQC_RAW ( reads ).html.set { fastqc_raw_html }
        fastqc_raw_zip = FASTQC_RAW.out.zip
        fastqc_version = FASTQC_RAW.out.version
    }

    trim_reads       = reads
    trim_json        = Channel.empty()
    trim_html        = Channel.empty()
    trim_log         = Channel.empty()
    trim_reads_fail  = Channel.empty()
    fastp_version    = Channel.empty()
    fastqc_trim_html = Channel.empty()
    fastqc_trim_zip  = Channel.empty()
    if (!params.skip_fastp) {
        FASTP ( reads ).reads.set { trim_reads }
        trim_json       = FASTP.out.json
        trim_html       = FASTP.out.html
        trim_log        = FASTP.out.log
        trim_reads_fail = FASTP.out.reads_fail
        fastp_version   = FASTP.out.version

        if (!params.skip_fastqc) {
            FASTQC_TRIM ( trim_reads ).html.set { fastqc_trim_html }
            fastqc_trim_zip = FASTQC_TRIM.out.zip
        }
    }

    emit:
    reads = trim_reads // channel: [ val(meta), [ reads ] ]
    trim_json          // channel: [ val(meta), [ json ] ]
    trim_html          // channel: [ val(meta), [ html ] ]
    trim_log           // channel: [ val(meta), [ log ] ]
    trim_reads_fail    // channel: [ val(meta), [ fastq.gz ] ]
    fastp_version      //    path: *.version.txt

    fastqc_raw_html    // channel: [ val(meta), [ html ] ]
    fastqc_raw_zip     // channel: [ val(meta), [ zip ] ]
    fastqc_trim_html   // channel: [ val(meta), [ html ] ]
    fastqc_trim_zip    // channel: [ val(meta), [ zip ] ]
    fastqc_version     //    path: *.version.txt
}
