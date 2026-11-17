import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    # 1. Read the expression dataset, print some basic information on raw input
    # The first 6 columns are metadata (below), the rest are gene expression values for each cell line 
    # '',SequencingID,ModelConditionID,ModelID,IsDefaultEntryForMC,IsDefaultEntryForModel
    # each row is an expression profile, some cell lines(ModelID) may have multiple profiles
    raw_data_by_profile = pd.read_csv(args.input)
    num_profiles, num_columns = raw_data_by_profile.shape

    # 2. Filter to default expression profiles
    raw_data_by_cell_line = raw_data_by_profile[
        raw_data_by_profile["IsDefaultEntryForModel"]
        .astype(str)
        .str.lower()
        .isin(["yes", "true"])
    ].copy()

    # 3. Verify uniqueness, set ModelID as index
    assert raw_data_by_cell_line["ModelID"].is_unique, "Duplicate models remain!"
    gene_expression_by_cell_lines = raw_data_by_cell_line.set_index("ModelID")

    # 4. Drop metadata columns, keep only gene expression values
    metadata_columns = [
        "SequencingID",
        "ModelConditionID",
        "ModelID",
        "IsDefaultEntryForMC",
        "IsDefaultEntryForModel"
    ]
    gene_expression_by_cell_lines = gene_expression_by_cell_lines.drop(
        columns=metadata_columns + ["Unnamed: 0"],
        errors="ignore")
    num_cell_lines, num_genes = gene_expression_by_cell_lines.shape

    # 5. Write summary information to output file
    with open(args.output, "w") as f:
        f.write(f"Number of raw columns:{num_columns}\n") 
        f.write(f"Number of raw profiles:{num_profiles}\n")
        f.write(f"Number of genes:{num_genes}\n") 
        f.write(f"Number of cell-line models:{num_cell_lines}\n")


if __name__ == "__main__":
    main()