# example usage:
# python3 batch_plot_ener.py -y Temp Pot Kin --tmin 5000 --tmax 10000 --prefix nvt

#!/usr/bin/env python3
import argparse
import os
import pandas as pd
import matplotlib
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
        description="Batch plot CP2K .ener data from multiple folders."
    )
    parser.add_argument(
        "-y", "--ycols",
        nargs="+",
        required=True,
        help="Column name(s) to plot"
    )
    parser.add_argument(
        "--folder-prefix",
        default="inas-water-nvt-",
        help="Folder prefix, default: inas-water-nvt-"
    )
    parser.add_argument(
        "--start",
        type=int,
        default=1,
        help="Starting folder index, default: 1"
    )
    parser.add_argument(
        "--end",
        type=int,
        default=10,
        help="Ending folder index, default: 10"
    )
    parser.add_argument(
        "--file-pattern",
        default="INAS110-2-2-4-water-nvt-1.ener",
        help='Energy filename inside each folder'
    )
    parser.add_argument(
        "--tmin",
        type=float,
        help="Minimum time in fs"
    )
    parser.add_argument(
        "--tmax",
        type=float,
        help="Maximum time in fs"
    )
    parser.add_argument(
        "--prefix",
        default="batch",
        help="Prefix for output plot filenames"
    )

    args = parser.parse_args()

    for ycol in args.ycols:
        if ycol not in COLS:
            raise ValueError(f"Column '{ycol}' not found. Available columns: {COLS}")

    if args.tmin is not None and args.tmax is not None and args.tmin > args.tmax:
        raise ValueError("tmin must be less than or equal to tmax.")

    for ycol in args.ycols:
        plt.figure(figsize=(10, 6))
        plotted_any = False

        for i in range(args.start, args.end + 1):
            folder = f"{args.folder_prefix}{i}"
            ener_file = os.path.join(folder, args.file_pattern)

            if not os.path.isdir(folder):
                print(f"Skipping {folder}: missing folder")
                continue

            if not os.path.isfile(ener_file):
                print(f"Skipping {folder}: missing file {ener_file}")
                continue

            try:
                df = read_cp2k_ener(ener_file)

                if args.tmin is not None:
                    df = df[df["Time_fs"] >= args.tmin]
                if args.tmax is not None:
                    df = df[df["Time_fs"] <= args.tmax]

                if df.empty:
                    print(f"Skipping {folder}: no data in selected time range")
                    continue

                plt.plot(df["Time_fs"], df[ycol], label=folder)
                plotted_any = True

            except Exception as e:
                print(f"Skipping {folder}: error reading {ener_file}: {e}")

        if not plotted_any:
            print(f"No data plotted for {ycol}")
            plt.close()
            continue

        plt.xlabel("Time (fs)")
        plt.ylabel(ycol)
        plt.title(f"{ycol} from multiple folders")
        plt.legend(fontsize=8, ncol=2)
        plt.tight_layout()
        plt.show()

        outfile = f"{args.prefix}_{ycol}.png"
        plt.savefig(outfile, dpi=300)
        plt.close()
        print(f"Saved plot to {outfile}")

if __name__ == "__main__":
    main()