import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.metrics.pairwise import cosine_similarity


def principal_component_analysis(
        gene_expression_by_cell_lines: pd.DataFrame,  
        args: argparse.Namespace) -> None:
    """
    Perform PCA on the gene expression data and plot the first two principal components.
    """
    # Tested with 50, 200, and 500 components, more components led to better clustering of cell lines by lineage
    pca = PCA(n_components=500, svd_solver="randomized", random_state=42)
    scores = pca.fit_transform(gene_expression_by_cell_lines.to_numpy(dtype=pd.Float32Dtype))

    # Save embeddings for later similarity calculations
    embeddings = pd.DataFrame(
        scores,
        index=gene_expression_by_cell_lines.index,
        columns=[f"PC{i+1}" for i in range(scores.shape[1])]
    )
    embeddings.to_csv(args.pca_output, index=True)

    # Compute cosine similarity between cell lines based on PCA embeddings
    similarity: np.ndarray = cosine_similarity(scores)
    pd.DataFrame(
        similarity,
        index=gene_expression_by_cell_lines.index,
        columns=gene_expression_by_cell_lines.index
    ).to_csv(args.pca_cosine_output, index=True)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed_expression_input", required=True)
    parser.add_argument("--pca_output", required=True)
    parser.add_argument("--pca_cosine_output", required=True)
    args = parser.parse_args()

    # read the cleaned gene expression data and metadata
    gene_expression_by_cell_lines = pd.read_csv(args.processed_expression_input, index_col="ModelID")

    #vp.plot_gene_variance(gene_expression_by_cell_lines, args)
    principal_component_analysis(gene_expression_by_cell_lines, args)

if __name__ == "__main__":
    main()