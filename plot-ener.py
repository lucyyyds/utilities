#!/usr/bin/env python3
import argparse
import pandas as pd
import matplotlib.pyplot as plt

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
        description="Plot CP2K .ener data using Time_fs as x-axis."
    )
    parser.add_argument("file", help="Input .ener file")
    parser.add_argument(
        "-y", "--ycols",
        nargs="+",
        required=True,
        help="Column name(s) to plot on y-axis"
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List available columns and exit"
    )
    parser.add_argument(
        "-o", "--output",
        help="Save figure to file instead of showing it"
    )
    args = parser.parse_args()

    df = read_cp2k_ener(args.file)

    if args.list:
        print("Available columns:")
        for col in df.columns:
            print(col)
        return

    xcol = "Time_fs"

    for ycol in args.ycols:
        if ycol not in df.columns:
            raise ValueError(
                f"Column '{ycol}' not found.\nAvailable columns: {list(df.columns)}"
            )

    plt.figure(figsize=(8, 5))
    for ycol in args.ycols:
        plt.plot(df[xcol], df[ycol], label=ycol)

    plt.xlabel("Time (fs)")
    plt.ylabel("Value")
    plt.title("CP2K .ener data")
    plt.legend()
    plt.tight_layout()

    if args.output:
        plt.savefig(args.output, dpi=300)
        print(f"Saved plot to {args.output}")
    else:
        plt.show()

if __name__ == "__main__":
    main()
