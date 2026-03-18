#!/usr/bin/env python3
import argparse
import pandas as pd
import matplotlib.pyplot as plt


def read_cp2k_ener(filename):
    # Read header line
    with open(filename, "r") as f:
        header_line = f.readline().strip()

    # Clean header names
    headers = header_line.lstrip("#").split()

    # Read data
    df = pd.read_csv(
        filename,
        delim_whitespace=True,
        comment="#",
        header=None,
        names=headers
    )

    return df


def main():
    parser = argparse.ArgumentParser(
        description="Plot CP2K .ener data using Time[fs] as x-axis."
    )
    parser.add_argument(
        "file",
        help="Input .ener file"
    )
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
    args = parser.parse_args()

    df = read_cp2k_ener(args.file)

    if args.list:
        print("Available columns:")
        for col in df.columns:
            print(col)
        return

    xcol = "Time[fs]"
    if xcol not in df.columns:
        raise ValueError(f"Column '{xcol}' not found in file.")

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
    plt.title("CP2K Energy Data")
    plt.legend()
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
