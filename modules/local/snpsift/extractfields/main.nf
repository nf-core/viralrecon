process SNPSIFT_EXTRACTFIELDS {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::snpsift=4.3.1t"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snpsift:4.3.1t--hdfd78af_3' :
        'quay.io/biocontainers/snpsift:4.3.1t--hdfd78af_3' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.snpsift.txt"), emit: txt
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def avail_mem = 4
    if (!task.memory) {
        log.info '[SnpSift] Available memory not known - defaulting to 4GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    """
    SnpSift \\
        -Xmx${avail_mem}g \\
        extractFields \\
        -s "," \\
        -e "." \\
        $args \\
        $vcf \\
        CHROM POS REF ALT \\
        "ANN[*].GENE" "ANN[*].GENEID" \\
        "ANN[*].IMPACT" "ANN[*].EFFECT" \\
        "ANN[*].FEATURE" "ANN[*].FEATUREID" \\
        "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" \\
        "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" \\
        "ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" \\
        "ANN[*].AA_LEN" "ANN[*].DISTANCE" "EFF[*].EFFECT" \\
        "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" \\
        > ${prefix}.snpsift.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpsift: \$( echo \$(SnpSift split -h 2>&1) | sed 's/^.*version //' | sed 's/(.*//' | sed 's/t//g' )
    END_VERSIONS
    """
}
