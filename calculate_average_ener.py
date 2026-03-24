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
        description="Calculate average of a selected column within a Time_fs range."
    )
    parser.add_argument("file", help="Input .ener file")
    parser.add_argument(
        "-y", "--ycol",
        required=True,
        help="Column name to average"
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

    if args.ycol not in df.columns:
        raise ValueError(
            f"Column '{args.ycol}' not found.\nAvailable columns: {list(df.columns)}"
        )

    if args.tmin > args.tmax:
        raise ValueError("tmin must be less than or equal to tmax.")

    df_range = df[(df["Time_fs"] >= args.tmin) & (df["Time_fs"] <= args.tmax)]

    if df_range.empty:
        raise ValueError(
            f"No data points found in the range {args.tmin} to {args.tmax} fs."
        )

    avg_value = df_range[args.ycol].mean()

    print(f"Selected column : {args.ycol}")
    print(f"Time range      : {args.tmin} to {args.tmax} fs")
    print(f"Points in range : {len(df_range)}")
    print(f"Average value   : {avg_value:.10f}")

if __name__ == "__main__":
    main()
