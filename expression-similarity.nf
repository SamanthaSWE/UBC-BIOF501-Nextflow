//why
nextflow.enable.dsl=2

process PreprocessExpression {
    publishDir 'results', mode: 'copy'

    input: 
    path expression_csv_input
    path python_script

    output:
    path "processed_expression_protein_coding_genes.csv", emit: data
    path "figures", emit: figures

    script:
    """
    python3 ${python_script} \
        --expression_input ${expression_csv_input} \
        --expression_output processed_expression_protein_coding_genes.csv \
        --figures figures
    """
}

process PCAAnalysis {
    publishDir 'results/pca', mode: 'copy'

    input:
    path processed_expression_csv_input
    path python_script

    output:
    tuple path("pca_results.csv"),
        path("pca_cosine_results.csv"),
        emit: pca

    script:
    """
    python3 ${python_script} \
        --processed_expression_input ${processed_expression_csv_input} \
        --pca_output pca_results.csv \
        --pca_cosine_output pca_cosine_results.csv 
    """
}

process PearsonAnalysis {
    publishDir 'results/pearson', mode: 'copy'

    input:
    path processed_expression_csv_input
    path human_protein_tsv_input
    path python_script

    output:
    tuple path("pearson_results.csv"),
        path("spearman_results.csv"),
        emit: correlations

    script:
    """
    python3 ${python_script} \
        --processed_expression_input ${processed_expression_csv_input} \
        --human_protein_atlas_input ${human_protein_tsv_input} \
        --pearson_output pearson_results.csv \
        --spearman_output spearman_results.csv 
    """
}

process SimilarityEvaluation {
    publishDir 'results/final_similarity', mode: 'copy'

    input:
    tuple path(pearson_results),
        path(spearman_results)
    tuple path(pca_results),
        path(pca_cosine_results)
    path metadata
    path python_script

    output:
    path "figures"

    script:
    """
    python3 ${python_script} \
        --pearson_results ${pearson_results} \
        --spearman_results ${spearman_results} \
        --pca_results ${pca_results} \
        --pca_cosine_results ${pca_cosine_results} \
        --metadata ${metadata} \
        --figures figures 
    """
}

workflow {
    expression = channel.fromPath(
        params.expression_file,
        checkIfExists: true
    )

    human_protein_atlas = channel.fromPath(
        params.human_protein_atlas_input_file,
        checkIfExists: true
    )

    metadata = channel.fromPath(
        params.metadata_file,
        checkIfExists: true
    )

    preprocessed_expression = PreprocessExpression(
        expression,
        file("${projectDir}/python/preprocess_expression.py")
    )

    pca = PCAAnalysis(
        preprocessed_expression.data,
        file("${projectDir}/python/pca_analysis.py")
    )

    pearson_spearman = PearsonAnalysis(
        preprocessed_expression.data,
        human_protein_atlas,
        file("${projectDir}/python/pearson_analysis.py")
    )

    SimilarityEvaluation(
        pearson_spearman.correlations,
        pca.pca,
        metadata,
        file("${projectDir}/python/similarity_evaluation.py")
    )

}