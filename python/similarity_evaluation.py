import argparse
import pandas as pd
import numpy as np
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt
import umap

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

def scatter_pca_basic(pca_data: pd.DataFrame, output_dir: Path) -> None:
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

def umap_scatter_plot_spearman(spearman: pd.DataFrame, output_dir: Path, model_to_lineage: dict[str, str]) -> None:
    
    distance = np.clip(1 - spearman.to_numpy(), 0, 2)
    np.fill_diagonal(distance, 0)

    reducer = umap.UMAP(
        n_components=2,
        metric="precomputed",
        n_neighbors=30,
        min_dist=0.1,
        random_state=42
    )
    embedding = reducer.fit_transform(distance)
    lineages = [
        model_to_lineage.get(mid, "Unknown")
        for mid in spearman.index
    ]

    plt.figure(figsize=(14, 9))
    sns.scatterplot(
        x=embedding[:, 0],
        y=embedding[:, 1],
        hue=lineages,
        palette="husl",
        s=8,
        alpha=0.7
    )

    rng = np.random.default_rng(42)
    for lineage in set(lineages):
        indices = np.where(np.array(lineages) == lineage)[0]
        selected = rng.choice(indices, size=min(2, len(indices)), replace=False)

        for i in selected:
            plt.annotate(
                lineage[:3],
                (embedding[i, 0], embedding[i, 1]),
                fontsize=8,
                fontweight="bold",
                xytext=(4, 4),
                textcoords="offset points"
            )

    plt.legend(
        title="Lineage",
        bbox_to_anchor=(0.5, -0.15),
        loc="upper center",
        ncol=5,
        fontsize=8,
        markerscale=2
    )
    plt.tight_layout()
    plt.xlabel("UMAP 1")
    plt.ylabel("UMAP 2")
    plt.title("UMAP of cell-line spearman similarity")
    plt.savefig(output_dir / f"spearman_umap_30.png", dpi=150, bbox_inches="tight")
    plt.close()

def umap_scatter_plot_cosine(pca_cosine: pd.DataFrame, output_dir: Path, model_to_lineage: dict[str, str]) -> None:
    cosine_distance = np.clip(1 - pca_cosine.to_numpy(), 0, 2)
    np.fill_diagonal(cosine_distance, 0)

    reducer = umap.UMAP(
        metric="precomputed",
        n_neighbors=30,
        min_dist=0.1,
        random_state=42
    )
    embedding = reducer.fit_transform(cosine_distance)
    lineages = [
        model_to_lineage.get(mid, "Unknown")
        for mid in pca_cosine.index
    ]

    plt.figure(figsize=(14, 9))
    sns.scatterplot(
        x=embedding[:, 0],
        y=embedding[:, 1],
        hue=lineages,
        palette="husl",
        s=8,
        alpha=0.7
    )

    rng = np.random.default_rng(42)

    for lineage in set(lineages):
        indices = np.where(np.array(lineages) == lineage)[0]
        selected = rng.choice(indices, size=min(2, len(indices)), replace=False)

        for i in selected:
            plt.annotate(
                lineage[:3],
                (embedding[i, 0], embedding[i, 1]),
                fontsize=8,
                fontweight="bold",
                xytext=(4, 4),
                textcoords="offset points"
            )

    plt.legend(
        title="Lineage",
        bbox_to_anchor=(0.5, -0.15),
        loc="upper center",
        ncol=5,
        fontsize=8,
        markerscale=2
    )
    plt.tight_layout()
    plt.xlabel("UMAP 1")
    plt.ylabel("UMAP 2")
    plt.title("UMAP of cell-line PCA cosine similarity")
    plt.savefig(output_dir / f"pca_umap_30.png", dpi=150, bbox_inches="tight")
    plt.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--spearman_results", required=True)
    parser.add_argument("--pca_cosine_results", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--figures", required=True)
    args = parser.parse_args()

    output_dir :Path = Path(args.figures)
    output_dir.mkdir(parents=True, exist_ok=True)

    # read all csv
    spearman_csv = pd.read_csv(args.spearman_results, index_col="ModelID")
    pca_cosine_results = pd.read_csv(args.pca_cosine_results, index_col="ModelID")
    lineage_metadata = pd.read_csv(args.metadata)

    # clustermap_spearman_correlation(spearman_csv, output_dir)
    # scatter_pca(pca_results, output_dir)
    # clustermap_pca_cosine_clustermap(pca_cosine_results, output_dir)

    # Evaluate Lineage Agreements
    model_to_lineage: dict[str, str] = (
        lineage_metadata.set_index("ModelID")["OncotreeLineage"]
        .fillna("Unknown")
        .to_dict()
    )

    umap_scatter_plot_cosine(pca_cosine_results, output_dir, model_to_lineage)
    umap_scatter_plot_spearman(spearman_csv, output_dir, model_to_lineage)
    # what is random lineage 
    lineages = pd.Series(model_to_lineage)
    lineages = lineages[lineages != "Unknown"]

    counts = lineages.value_counts()
    n = counts.sum()
    random_baseline = (counts * (counts - 1)).sum() / (n * (n - 1))

    similarity_methods = {
        "Spearman": spearman_csv,
        "PCA Cosine": pca_cosine_results,
    }

    unknown_ids: set[str] = {
        mid for mid, lineage in model_to_lineage.items()
        if lineage == "Unknown"
    }
    similarity_analysis = {}
    records = []

    for method, data in similarity_methods.items():
        # Keep only cell lines with known lineage, in matching row/column order
        known_ids = [
            mid for mid in data.index
            if mid not in unknown_ids
        ]
        scores = data.loc[known_ids, known_ids].to_numpy()
        lineages = np.array([model_to_lineage[i] for i in known_ids])

        # Exclude self-comparisons
        np.fill_diagonal(
            scores,
            -np.inf
        )
        # Find top 10 neighbors without sorting the entire matrix
        neighbors = np.argpartition(scores, -10, axis=1)[:, -10:]
        agreement = (lineages[neighbors] == lineages[:, None]).mean()
        similarity_analysis[method] = agreement


        for mid in unknown_ids:
            if mid not in data.index:
                continue

            neighbors = data.loc[mid, known_ids].nlargest(10)

            for neighbor in neighbors.index:
                records.append({
                    "ModelID": mid,
                    "Method": method,
                    "Lineage": model_to_lineage[neighbor]
                })

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

    counts = (
        pd.DataFrame(records)
        .groupby(["ModelID", "Method", "Lineage"])
        .size()
        .unstack(fill_value=0)
    )

    colors = sns.color_palette("husl", n_colors=len(counts.columns))

    ax = counts.plot(
        kind="bar",
        stacked=True,
        figsize=(15, 7),
        width=0.8,
        color=colors
    )

    for lineage, container in zip(counts.columns, ax.containers):
        short_lineage = lineage[:3]

        labels = [
            short_lineage if v >= 1 else ""
            for v in container.datavalues
        ]

        ax.bar_label(
            container,
            labels=labels,
            label_type="center",
            fontsize=7,
            rotation=90
        )

    plt.ylabel("Number of top-10 neighbors")
    plt.xlabel("Unknown ModelID / Method")
    plt.title("Predicted lineages of unknown cell lines")
    plt.legend(title="Lineage", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(output_dir / f"unknown_lineage_predictions.png", dpi=150)
    plt.close()

if __name__ == "__main__":
    main()