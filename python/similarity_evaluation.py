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
    plot_data = pca_data[["PC1", "PC2"]].copy()

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
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title(f"PCA")
    plt.tight_layout()
    plt.savefig(output_dir / f"pca1_pca2.png", dpi=150)
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pearson_results", required=True)
    parser.add_argument("--spearman_results", required=True)
    parser.add_argument("--pca_results", required=True)
    parser.add_argument("--pca_cosine_results", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--figures", required=True)
    args = parser.parse_args()

    output_dir :Path = Path(args.figures)
    output_dir.mkdir(parents=True, exist_ok=True)

    # read all csv
    pearson_csv = pd.read_csv(args.pearson_results, index_col="ModelID")
    spearman_csv = pd.read_csv(args.spearman_results, index_col="ModelID")
    pca_results = pd.read_csv(args.pca_results, index_col="ModelID")
    pca_cosine_results = pd.read_csv(args.pca_cosine_results, index_col="ModelID")
    metadata = pd.read_csv(args.metadata)

    # pearson and spearman visualizations
    # clustermap_pearson_correlation(pearson_csv, output_dir)

    # clustermap_spearman_correlation(spearman_csv, output_dir)

    # # pearson vs spearman visualizations
    # scatter_pearson_vs_spearman_correlation(
    #     pearson_matrix=pearson_csv,
    #     spearman_matrix=spearman_csv,
    #     output_dir=output_dir)
    # scatter_pearson_vs_spearman_correlation(
    #         pearson_matrix=pearson_csv,
    #         spearman_matrix=spearman_csv,
    #         output_dir=output_dir)

    # pca visualizations
    #scatter_pca(pca_results, output_dir)
    #clustermap_pca_cosine_clustermap(pca_cosine_results, output_dir)

    # REAL WORK
    model_to_lineage: dict[str, str] = (
        metadata.set_index("ModelID")["OncotreeLineage"]
        .fillna("Unknown")
        .to_dict()
    )
    # what is random lineage 
    lineages = pd.Series(model_to_lineage)
    lineages = lineages[lineages != "Unknown"]

    counts = lineages.value_counts()
    n = counts.sum()
    random_baseline = (counts * (counts - 1)).sum() / (n * (n - 1))

    results = {
        "Pearson": pearson_csv,
        "Spearman": spearman_csv,
        "PCA Cosine": pca_cosine_results,
    }

    similarity_analysis = {}

    for method, data in results.items():
        # Keep only cell lines with known lineage, in matching row/column order
        ids = [
            model_id for model_id in data.index
            if model_to_lineage.get(model_id, "Unknown") != "Unknown"
        ]
        scores = data.loc[ids, ids].to_numpy()
        lineages = np.array([model_to_lineage[i] for i in ids])

        # Exclude self-comparisons
        np.fill_diagonal(
            scores,
            -np.inf
        )
        # Find top 10 neighbors without sorting the entire matrix
        neighbors = np.argpartition(scores, -10, axis=1)[:, -10:]
        agreement = (lineages[neighbors] == lineages[:, None]).mean()
        similarity_analysis[method] = agreement

    plt.figure(figsize=(8, 5))
    plt.bar(similarity_analysis.keys(), similarity_analysis.values())
    plt.ylabel("Top-10 lineage agreement")
    plt.title("Similarity method performance by lineage")
    plt.ylim(0, 1)
    plt.xticks(rotation=20)
    plt.axhline(
        y=random_baseline,
        color="red",
        linestyle="--",
        label=f"Random baseline ({random_baseline:.1%})"
    )
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_dir / f"method_lineage_agreement.png", dpi=150)
    plt.close()


if __name__ == "__main__":
    main()