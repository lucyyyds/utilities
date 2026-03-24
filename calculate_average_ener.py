#!/usr/bin/env python3
import argparse
import pandas as pd

COLS = [
    "Step",
    "Time_fs",
    "Kin",
    "Temp",
    "Pot",
    "ConsQty",
    "UsedTime_s",
]

def read_cp2k_ener(filename):
    df = pd.read_csv(
        filename,
        comment="#",
        sep=r"\s+",
        header=None,
        names=COLS,
        engine="python",
    )
    return df

def main():
    parser = argparse.ArgumentParser(
        description="Calculate average of one or more selected columns within a Time_fs range."
    )
    parser.add_argument("file", help="Input .ener file")
    parser.add_argument(
        "-y", "--ycol",
        nargs="+",
        required=True,
        help="One or more column names to average"
    )
    parser.add_argument(
        "--tmin",
        type=float,
        required=True,
        help="Minimum time in fs"
    )
    parser.add_argument(
        "--tmax",
        type=float,
        required=True,
        help="Maximum time in fs"
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List available columns and exit"
    )

    args = parser.parse_args()

    df = read_cp2k_ener(args.file)

    if args.list:
        print("Available columns:")
        for col in df.columns:
            print(col)
        return

    if args.tmin > args.tmax:
        raise ValueError("tmin must be less than or equal to tmax.")

    # Check all requested columns exist
    missing_cols = [col for col in args.ycol if col not in df.columns]
    if missing_cols:
        raise ValueError(
            f"Column(s) not found: {missing_cols}\nAvailable columns: {list(df.columns)}"
        )

    df_range = df[(df["Time_fs"] >= args.tmin) & (df["Time_fs"] <= args.tmax)]

    if df_range.empty:
        raise ValueError(
            f"No data points found in the range {args.tmin} to {args.tmax} fs."
        )

    print(f"Time range      : {args.tmin} to {args.tmax} fs")
    print(f"Points in range : {len(df_range)}")
    print("Average values:")

    for col in args.ycol:
        avg_value = df_range[col].mean()
        print(f"  {col:<12} {avg_value:.10f}")

if __name__ == "__main__":
    main()
