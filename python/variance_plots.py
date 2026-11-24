import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

def plot_gene_variance(gene_expression_by_cell_lines: pd.DataFrame, args: argparse.Namespace) -> None:
    """
    param gene_expression_by_cell_lines: DataFrame
        Rows = cell lines / Model ID (1 per)
        Columns = genes
        Values = log2(TPM + 1)
    output to results:
    """
    output_dir = Path(args.figures)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Calculate statistics across cell lines
    gene_means = gene_expression_by_cell_lines.mean(axis=0)
    gene_variances = gene_expression_by_cell_lines.var(axis=0)

    # 1. Histogram of gene variances
    plt.figure(figsize=(9, 5))
    plt.hist(gene_variances, bins=200) # reasonable bin size to see distribution of variances (~19000 genes)
    plt.xlabel("Gene expression variance")
    plt.ylabel("Number of genes")
    plt.title("Distribution of gene variances")
    plt.tight_layout()
    plt.savefig(output_dir / "gene_variance_histogram.png", dpi=150) # image resolution
    plt.close()

    # 2. Mean-variance relationship
    # All genes in grey
    plt.figure(figsize=(9, 5))
    plt.scatter(
        gene_means,
        gene_variances,
        s=5, # plot configuration: size of points
        alpha=0.15, # plot configuration: transparency of points
        color="grey",
        label="All genes"
    )
    # Selected HVGs in red
    top_genes = gene_variances.nlargest(5000).index # TODO: why 5000?
    plt.scatter(
        gene_means.loc[top_genes],
        gene_variances.loc[top_genes],
        s=5, # plot configuration: size of points
        alpha=0.15, # plot configuration: transparency of points
        color="red",
        label="Top 5,000 HVGs"
    )

    plt.xlabel("Mean log2(TPM + 1)")
    plt.ylabel("Gene expression variance")
    plt.title("Mean–variance relationship")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir / "mean_variance_scatter.png", dpi=150)
    plt.close()

    # Save statistics for later HVG selection
    stats = gene_expression_by_cell_lines.columns.to_frame(index=False)
    stats.columns = ["gene"]
    stats["mean"] = gene_means.to_numpy()
    stats["variance"] = gene_variances.to_numpy()
    stats.to_csv(output_dir / "gene_statistics.csv", index=False)