//why
nextflow.enable.dsl=2

process ExpressionSimilarity {
    publishDir 'results', mode: 'copy'

    input: 
    path expression_csv
    path python_script

    output:
    path "summary.tsv"

    script:
    """
    python3 ${python_script} \
        --input ${expression_csv} \
        --output summary.tsv 
    """
}

workflow {
    expression = channel.fromPath(
        params.expression_file,
        checkIfExists: true
    )
    ExpressionSimilarity(
        expression,
        file("${projectDir}/python/main.py")
    )
}