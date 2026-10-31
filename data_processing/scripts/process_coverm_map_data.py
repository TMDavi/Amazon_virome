import pandas as pd
import argparse


def merge_coverm_data(filelist):

    mean_df = None
    count_df = None
    coverage_df = None

    for file in filelist:

        temp_df = pd.read_csv(file, sep="\t")

        # First column contains contig IDs
        contig_col = temp_df.columns[0]

        # Identify CoverM columns
        mean_col = next(
            col for col in temp_df.columns
            if col.endswith(" Mean")
        )

        count_col = next(
            col for col in temp_df.columns
            if col.endswith(" Read Count")
        )

        coverage_col = next(
            col for col in temp_df.columns
            if col.endswith(" Covered Fraction")
        )

        # Extract sample ID from column name
        sampleID = mean_col.removesuffix(" Mean")

        # Create temporary tables
        mean_temp = pd.DataFrame({
            "contig_id": temp_df[contig_col],
            sampleID: temp_df[mean_col]
        })

        count_temp = pd.DataFrame({
            "contig_id": temp_df[contig_col],
            sampleID: temp_df[count_col]
        })

        coverage_temp = pd.DataFrame({
            "contig_id": temp_df[contig_col],
            sampleID: temp_df[coverage_col]
        })

        # Merge mean
        if mean_df is None:
            mean_df = mean_temp
        else:
            mean_df = pd.merge(
                mean_df,
                mean_temp,
                on="contig_id",
                how="outer"
            )

        # Merge count
        if count_df is None:
            count_df = count_temp
        else:
            count_df = pd.merge(
                count_df,
                count_temp,
                on="contig_id",
                how="outer"
            )

        # Merge covered fraction
        if coverage_df is None:
            coverage_df = coverage_temp
        else:
            coverage_df = pd.merge(
                coverage_df,
                coverage_temp,
                on="contig_id",
                how="outer"
            )

    return mean_df, count_df, coverage_df


def main():

    parser = argparse.ArgumentParser(
        description="Merge CoverM results into separate matrices."
    )

    parser.add_argument(
        "-i",
        "--input",
        nargs="+",
        required=True,
        help="List of CoverM output files"
    )

    parser.add_argument(
        "-o",
        "--output",
        default="coverm",
        help="Output prefix"
    )

    args = parser.parse_args()

    mean_df, count_df, coverage_df = merge_coverm_data(args.input)

    mean_df.to_csv(
        f"{args.output}_mean.tsv",
        sep="\t",
        index=False
    )

    count_df.to_csv(
        f"{args.output}_count.tsv",
        sep="\t",
        index=False
    )

    coverage_df.to_csv(
        f"{args.output}_coverage.tsv",
        sep="\t",
        index=False
    )


if __name__ == "__main__":
    main()