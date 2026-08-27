process BLAST_REPORT {
    label 'process_single'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/63/63c2848a35d4087421472ea90c04148fc97c4ccc2839c76f7cb3919458bb10ef/data' :
        'community.wave.seqera.io/library/jinja2_matplotlib_pandas_python_seaborn:35dc011346333319' }"

    input:
    tuple val(meta), path(blast), path(fasta)

    output:
    tuple val(meta), path("*.html")          , emit: blast_report
    tuple val(meta), path("*.fa")            , emit: reversed_contigs
    tuple val(meta), path("*_genotype.csv")  , emit: genotype
    tuple val("${task.process}"), val("python"), eval("python --version 2>&1 | sed 's/Python //g'"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val("matplotlib"), eval("python -c 'import matplotlib; print(matplotlib.__version__)'"), topic: versions, emit: versions_matplotlib
    tuple val("${task.process}"), val("seaborn"), eval("python -c 'import seaborn; print(seaborn.__version__)'"), topic: versions, emit: versions_seaborn
    tuple val("${task.process}"), val("pandas"), eval("python -c 'import pandas; print(pandas.__version__)'"), topic: versions, emit: versions_pandas
    tuple val("${task.process}"), val("jinja2"), eval("python -c 'import jinja2; print(jinja2.__version__)'"), topic: versions, emit: versions_jinja2

    when:
    task.ext.when == null || task.ext.when

    script:  // This script is bundled with the pipeline, in nf-core/viralrecon/bin/
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id

    """
    blast_report.py \\
        --blast_file ${blast} \\
        --fasta_file ${fasta} \\
        --sample_name ${meta.id} \\
        --output_html ${prefix}_blast_report.html \\
        --output_fasta ${prefix}_reversed_contigs.fa \\
        --output_genotype ${prefix}_genotype.csv \\
        $args
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.id
    """
    touch  ${prefix}_blast_report.html
    touch  ${prefix}_reversed_contigs.fa
    touch  ${prefix}_genotype.csv
    """
}
