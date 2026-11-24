import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt

def clustermap_pearson_correlation(pearson_matrix: pd.DataFrame, output_dir: Path) -> None:
    sns.clustermap(
        pearson_matrix,
        cmap="vlag",
        vmin=-1,
        vmax=1,
        figsize=(12, 12)
    )
    plt.savefig(output_dir / f"pearson_clustermap_filtered_tissue_categories.png", dpi=150)
    plt.close()

def clustermap_spearman_correlation(spearman_matrix: pd.DataFrame, output_dir: Path)-> None:
    # visualize spearman
    sns.clustermap(
        spearman_matrix,
        cmap="vlag",
        vmin=-1,
        vmax=1,
        figsize=(12, 12)
    )
    plt.savefig(output_dir / f"spearman_clustermap_filtered_tissue_categories.png", dpi=150)
    plt.close()

def scatter_pearson_vs_spearman_correlation(
        pearson_matrix: pd.DataFrame, 
        spearman_matrix: pd.DataFrame,
        output_dir: Path)-> None:
    # Use upper triangle to avoid duplicates and self-correlations
    idx = np.triu_indices(len(pearson_matrix), k=1)

    plt.scatter(
        pearson_matrix.to_numpy()[idx],
        spearman_matrix.to_numpy()[idx],
        alpha=0.2,
        s=5
    )

    plt.plot([-1, 1], [-1, 1], "--", color="gray")
    plt.xlabel("Pearson correlation")
    plt.ylabel("Spearman correlation")
    plt.xlim(-1, 1)
    plt.ylim(-1, 1)
    plt.savefig(output_dir / f"pearson_vs_spearman_filtered_tissue_categories.png", dpi=150)
    plt.close()

def histogram_pearson_vs_spearman_correlation(
        pearson_matrix: pd.DataFrame, 
        spearman_matrix: pd.DataFrame,
        output_dir: Path)-> None:
    idx = np.triu_indices(len(pearson_matrix), k=1)

    sns.histplot(
    pearson_matrix.to_numpy()[idx],
    label="Pearson",
    stat="density",
    element="step",
    fill=False
    )
    sns.histplot(
        spearman_matrix.to_numpy()[idx],
        label="Spearman",
        stat="density",
        element="step",
        fill=False
    )
    plt.legend()
    plt.savefig(output_dir / f"pearson_vs_spearman_hist_filtered_tissue_categories.png", dpi=150)
    plt.close()

def clustermap_pca_cosine_clustermap(pca_cosine_data: pd.DataFrame, output_dir: Path) -> None:
    pca_cosine_data = sns.clustermap(
        pca_cosine_data,
        cmap="viridis",
        vmin=-1,
        vmax=1,
        figsize=(10, 8),
        xticklabels=True,
        yticklabels=True,
    )
    plt.savefig(
        output_dir / f"cosine_similarity_clustered.png",
        dpi=150
    )
    plt.close() 

def scatter_pca(pca_data: pd.DataFrame, output_dir: Path) -> None:
    # Plot first two principal components
    plot_data = pd.DataFrame(
        {
            "PC1": pca_data[:, 0],
            "PC2": pca_data[:, 1],
        },
        index=pca_data.index
    )

    plt.figure(figsize=(9, 5))
    sns.scatterplot(
        data=plot_data,
        x="PC1",
        y="PC2",
        palette="husl",
        s=35,
        alpha=0.7,
        legend="full",
    )
    plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    plt.xlabel(f"PC1 ({pca_data.explained_variance_ratio_[0]:.1%})")
    plt.ylabel(f"PC2 ({pca_data.explained_variance_ratio_[1]:.1%})")
    plt.title(f"PCA")
    plt.tight_layout()
    plt.savefig(output_dir / f"pca1_pca2.png", dpi=150)
    plt.close()

def 

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pearson_results", required=True)
    parser.add_argument("--spearman_results", required=True)
    parser.add_argument("--pca_results", required=True)
    parser.add_argument("--pca_cosine_results", required=True)
    parser.add_argument("--protein_atlas", required=True)
    parser.add_argument("--figures", required=True)
    args = parser.parse_args()

    output_dir :Path = Path(args.figures)
    output_dir.mkdir(parents=True, exist_ok=True)

    # pearson and spearman visualizations
    pearson_csv = pd.read_csv(args.pearson_results, index_col="ModelID")
    clustermap_pearson_correlation(pearson_csv, output_dir)

    spearman_csv = pd.read_csv(args.spearman_results, index_col="ModelID")
    clustermap_spearman_correlation(spearman_csv, output_dir)

    # pearson vs spearman visualizations
    scatter_pearson_vs_spearman_correlation(
        pearson_matrix=pearson_csv,
        spearman_matrix=spearman_csv,
        output_dir=output_dir)
    scatter_pearson_vs_spearman_correlation(
            pearson_matrix=pearson_csv,
            spearman_matrix=spearman_csv,
            output_dir=output_dir)

    # pca visualizations
    pca_results = pd.read_csv(args.pca_results, index_col="ModelID")
    scatter_pca(pca_results, output_dir)

    pca_cosine_results = pd.read_csv(args.pca_cosine_results, index_col="ModelID")
    clustermap_pca_cosine_clustermap(pca_cosine_results, output_dir)

    # REAL WORK
    # get model information

if __name__ == "__main__":
    main()