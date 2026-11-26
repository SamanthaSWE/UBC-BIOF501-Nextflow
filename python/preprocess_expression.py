import argparse
import pandas as pd
import constants
import variance_plots as vp


def remove_metadata_and_multiple_profiles(raw_data_by_profile: pd.DataFrame, args: argparse.Namespace
                                          ) -> pd.DataFrame:
    """
    Remove raw data's metadata columns and filter to default expression profiles.
    param raw_data_by_profile: DataFrame
        Rows = profiles (roughly 1 per cell line, some cell lines have multiple profiles)
        Columns = a few metadata, genes
        Values = metadata values,log2(TPM + 1)
    output to results:
    """

    # 1. Filter to default expression profiles
    raw_data_by_cell_line = raw_data_by_profile[
        raw_data_by_profile["IsDefaultEntryForModel"]
        .astype(str)
        .str.lower()
        .isin(["yes", "true"])
    ].copy()

    # 2. Verify uniqueness, set ModelID as index
    assert raw_data_by_cell_line["ModelID"].is_unique, "Duplicate models remain!"
    gene_expression_by_cell_lines = raw_data_by_cell_line.set_index("ModelID")

    # 4. Drop metadata columns, keep only gene expression values (only numeric, exclude incomplete genes)
    gene_expression_by_cell_lines = gene_expression_by_cell_lines.drop(
        columns=constants.OMICSEXPRESSIONTPMLOGP1_METADATA_COLUMNS + ["Unnamed: 0"],
        errors="ignore")
    gene_expression_by_cell_lines = gene_expression_by_cell_lines.apply(pd.to_numeric, errors="coerce").dropna(axis=1)

    # 5. Log Summary Information
    num_profiles, num_columns = raw_data_by_profile.shape
    num_cell_lines, num_genes = gene_expression_by_cell_lines.shape

    print(f"Number of raw columns:{num_columns}\n", flush=True) 
    print(f"Number of raw profiles:{num_profiles}\n", flush=True)
    print(f"Number of genes:{num_genes}\n", flush=True) 
    print(f"Number of cell-line models:{num_cell_lines}\n", flush=True)

    gene_expression_by_cell_lines.to_csv(args.expression_output,index=True)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--expression_input", required=True)
    parser.add_argument("--expression_output", required=True)
    args = parser.parse_args()

    raw_data_by_profile = pd.read_csv(args.expression_input)
    remove_metadata_and_multiple_profiles(raw_data_by_profile, args)

if __name__ == "__main__":
    main()