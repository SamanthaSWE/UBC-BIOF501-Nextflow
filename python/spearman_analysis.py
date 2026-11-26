import argparse
import pandas as pd
import re
import constants

HPA_GENE_COL: str = "Gene"
HPA_GENE_SYN_COL: str = "Gene synonym"
HPA_RNA_TISSUE_SPECIFICITY_COL: str = "RNA tissue specificity"

def extract_gene_symbol(column_name: str) -> str:
    # Example: "TP53 (7157)" -> "TP53"
    return re.sub(r"\s+\(\d+\)$", "", column_name).strip().upper()

def get_tissue_specific_gene_symbols(hpa_all: pd.DataFrame) -> tuple[set[str] ,set[str]]:
    # proper gene names
    all_hpa_genes: set[str] = set(hpa_all[HPA_GENE_COL].dropna().str.strip().str.upper())

    # Filter HPA to tissue-specific genes
    hpa_filtered: pd.DataFrame = hpa_all[hpa_all[HPA_RNA_TISSUE_SPECIFICITY_COL].isin(constants.TISSUE_CATEGORIES)]
    filtered_gene_names = set(hpa_filtered[HPA_GENE_COL].dropna().str.strip().str.upper())

    # exclude synonyms that collide with another official gene symbol
    gene_synonyms = {
        synonym.strip().upper()
        for value in hpa_filtered[HPA_GENE_SYN_COL].dropna()
        for synonym in re.split(r"[,;]", value)
        if synonym.strip() and synonym.strip().upper() not in all_hpa_genes
    }

    return filtered_gene_names, gene_synonyms

def spearman_analysis(
        gene_expression_by_cell_lines: pd.DataFrame, 
        tissue_genes: tuple[set[str] ,set[str]],
        args: argparse.Namespace) -> None:
    tissue_gene_names, tissue_gene_synonyms = tissue_genes
    exact_matches = []
    synonym_matches = []
    seen_symbols = set()

    for gene_expr_column in gene_expression_by_cell_lines.columns:
        gene_expr_symbol = extract_gene_symbol(gene_expr_column)
        if gene_expr_symbol in seen_symbols: 
            continue

        if gene_expr_symbol in tissue_gene_names:
            exact_matches.append(gene_expr_column)
        elif gene_expr_symbol in tissue_gene_synonyms:
            synonym_matches.append(gene_expr_column)
        else:
            seen_symbols.add(gene_expr_column)

    filtered_tissue_gene_expression = gene_expression_by_cell_lines[
        exact_matches + synonym_matches
    ]

    spearman_matrix = filtered_tissue_gene_expression.T.corr(
        method="spearman"
    )
    spearman_matrix.to_csv(args.spearman_output)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processed_expression_input", required=True)
    parser.add_argument("--human_protein_atlas_input", required=True)
    parser.add_argument("--spearman_output", required=True)
    args = parser.parse_args()

    # read the cleaned gene expression data and metadata
    gene_expression_by_cell_lines = pd.read_csv(args.processed_expression_input, index_col="ModelID")
    # read all hpa from csv
    hpa_all: pd.DataFrame = pd.read_csv(
        args.human_protein_atlas_input, 
        sep="\t", 
        low_memory=False,
        usecols=[HPA_GENE_COL, HPA_GENE_SYN_COL, HPA_RNA_TISSUE_SPECIFICITY_COL]
    )
    # process hpa data
    gene_names, gene_synonyms = get_tissue_specific_gene_symbols(hpa_all)

    spearman_analysis(gene_expression_by_cell_lines, (gene_names, gene_synonyms), args)

if __name__ == "__main__":
    main()