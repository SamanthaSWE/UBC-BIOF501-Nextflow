

import argparse
import pandas as pd
import re
import constants
from time import perf_counter

def extract_gene_symbol(column_name: str) -> str:
    # Example: "TP53 (7157)" -> "TP53"
    return re.sub(r"\s+\(\d+\)$", "", column_name).strip()

def process_hpa_data(args: argparse.Namespace) -> tuple[set[str] ,set[str]]:
    hpa = pd.read_csv(
        args.human_protein_atlas_input, 
        sep="\t", 
        low_memory=False,
        usecols=[
            "Gene",
            "Gene synonym",
            "RNA tissue specificity"
        ])

    # Filter HPA to tissue-specific genes
    hpa = hpa[
        hpa["RNA tissue specificity"].isin(constants.TISSUE_CATEGORIES)
    ]

    gene_names = set(hpa["Gene"].dropna())

    gene_synonyms = {
        synonym.strip()
        for value in hpa["Gene synonym"].dropna()
        for synonym in re.split(r"[,;]", value)
        if synonym.strip()
    }

    return gene_names, gene_synonyms

def pearson_analysis(gene_expression_by_cell_lines: pd.DataFrame, args: argparse.Namespace) -> None:
    start = perf_counter()
    gene_names, gene_synonyms = process_hpa_data(args)

    exact_matches = []
    synonym_matches = []
    unmatched = []

    for column in gene_expression_by_cell_lines.columns:
        symbol = re.sub(r"\s+\(\d+\)$", "", column).strip()

        if symbol in gene_names:
            exact_matches.append(column)
        elif symbol in gene_synonyms:
            synonym_matches.append(column)
        else:
            unmatched.append(column)

    # Primary names vs synonyms
    overlap = gene_names & gene_synonyms
    print(f"Name-synonym overlaps: {len(overlap)}", flush=True)


    filtered_expression = gene_expression_by_cell_lines[
        exact_matches + synonym_matches
    ]

    print(f"Filter: {perf_counter() - start:.2f}s", flush=True)
    print(f"Exact matches: {len(exact_matches)}", flush=True)
    print(f"Synonym matches: {len(synonym_matches)}", flush=True)
    print(f"Excluded genes: {len(unmatched)}", flush=True)
    print(f"Retained genes: {filtered_expression.shape[1]}", flush=True)
    
    start = perf_counter()
    pearson_matrix = filtered_expression.T.corr(
        method="pearson"
    )
    pearson_matrix.to_csv(args.pearson_output)
    print(f"Pearson: {perf_counter() - start:.2f}s", flush=True)

    start = perf_counter()
    spearman_matrix = filtered_expression.T.corr(
        method="spearman"
    )
    spearman_matrix.to_csv(args.spearman_output)
    print(f"Spearman: {perf_counter() - start:.2f}s", flush=True)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed_expression_input", required=True)
    parser.add_argument("--human_protein_atlas_input", required=True)
    parser.add_argument("--pearson_output", required=True)
    parser.add_argument("--spearman_output", required=True)
    args = parser.parse_args()

    # read the cleaned gene expression data and metadata
    gene_expression_by_cell_lines = pd.read_csv(args.processed_expression_input, index_col="ModelID")

    #vp.plot_gene_variance(gene_expression_by_cell_lines, args)
    pearson_analysis(gene_expression_by_cell_lines, args)

if __name__ == "__main__":
    main()