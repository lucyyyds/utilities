# example usage:
# python3 batch_avg_ener.py -y Temp Pot Kin --tmin 5000 --tmax 10000

#!/usr/bin/env python3
import argparse
import os
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
        description="Calculate average values from CP2K .ener files across multiple folders."
    )
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
        help='Energy filename inside each folder, default: "INAS110-2-2-4-water-nvt-1.ener"'
    )
    parser.add_argument(
        "--output",
        help="Optional output CSV filename"
    )

    args = parser.parse_args()

    if args.tmin > args.tmax:
        raise ValueError("tmin must be less than or equal to tmax.")

    results = []

    for i in range(args.start, args.end + 1):
        folder = f"{args.folder_prefix}{i}"
        ener_file = os.path.join(folder, args.file_pattern)

        row = {
            "Folder": folder,
            "File": ener_file,
        }

        if not os.path.isdir(folder):
            row["Status"] = "Missing folder"
            results.append(row)
            continue

        if not os.path.isfile(ener_file):
            row["Status"] = "Missing file"
            results.append(row)
            continue

        try:
            df = read_cp2k_ener(ener_file)

            missing_cols = [col for col in args.ycol if col not in df.columns]
            if missing_cols:
                row["Status"] = f"Missing columns: {missing_cols}"
                results.append(row)
                continue

            df_range = df[(df["Time_fs"] >= args.tmin) & (df["Time_fs"] <= args.tmax)]

            if df_range.empty:
                row["Status"] = "No data in range"
                results.append(row)
                continue

            row["Status"] = "OK"
            row["Points"] = len(df_range)

            for col in args.ycol:
                row[f"Avg_{col}"] = df_range[col].mean()

            results.append(row)

        except Exception as e:
            row["Status"] = f"Error: {e}"
            results.append(row)

    result_df = pd.DataFrame(results)

    avg_cols = [f"Avg_{col}" for col in args.ycol]
    ok_df = result_df[result_df["Status"] == "OK"]

    summary_row = {
        "Folder": "OVERALL_MEAN",
        "File": "-",
        "Status": f"{len(ok_df)} folders"
    }

    if not ok_df.empty:
        summary_row["Points"] = ok_df["Points"].mean() if "Points" in ok_df.columns else None
        for col in avg_cols:
            if col in ok_df.columns:
                summary_row[col] = ok_df[col].mean()
    else:
        summary_row["Points"] = None
        for col in avg_cols:
            summary_row[col] = None

    result_df = pd.concat([result_df, pd.DataFrame([summary_row])], ignore_index=True)

    pd.set_option("display.max_columns", None)
    pd.set_option("display.width", 200)
    pd.set_option("display.colheader_justify", "center")

    print(f"\nTime range: {args.tmin} to {args.tmax} fs\n")
    print(result_df.to_string(index=False))

    if args.output:
        result_df.to_csv(args.output, index=False)
        print(f"\nSaved table to {args.output}")

if __name__ == "__main__":
    main()