#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib.pyplot as plt


def convert_units(data, unit, molar_mass):
    if unit == "nm-3":
        return data, "Number density (particles / nm³)"
    elif unit == "kg/m3":
        factor = (molar_mass * 1e-3) / (6.02214076e23 * 1e-27)
        return data * factor, "Mass density (kg/m³)"
    elif unit == "g/cm3":
        factor = molar_mass / (6.02214076e23 * 1e-21)
        return data * factor, "Mass density (g/cm³)"
    else:
        raise ValueError("Unknown unit. Choose from: nm-3, kg/m3, g/cm3")


def read_xyz_trajectory(filename):
    frames = []
    with open(filename, "r") as f:
        while True:
            line = f.readline()
            if not line:
                break

            line = line.strip()
            if not line:
                continue

            natoms = int(line)
            _comment = f.readline()

            symbols = []
            coords = []
            for _ in range(natoms):
                parts = f.readline().split()
                if len(parts) < 4:
                    raise ValueError(f"Malformed XYZ file: {filename}")
                symbols.append(parts[0])
                coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

            frames.append((symbols, np.array(coords, dtype=float)))

    if not frames:
        raise ValueError(f"No frames found in {filename}")
    return frames


def replicate_xy(symbols, coords, Lx, Ly, nx_rep=3, ny_rep=2):
    replicated_coords = []
    replicated_symbols = []

    for ix in range(nx_rep):
        for iy in range(ny_rep):
            shift = np.array([ix * Lx, iy * Ly, 0.0])
            replicated_coords.append(coords + shift)
            replicated_symbols.extend(symbols)

    replicated_coords = np.vstack(replicated_coords)
    return replicated_symbols, replicated_coords


def main():
    p = argparse.ArgumentParser(
        description="Compute and plot 2D density map from multiple CP2K XYZ trajectories."
    )
    p.add_argument(
        "xyzfiles",
        nargs="+",
        help="One or more multi-frame XYZ trajectories from CP2K"
    )
    p.add_argument("--Lx", type=float, required=True, help="Original box length in x (nm)")
    p.add_argument("--Ly", type=float, required=True, help="Original box length in y (nm)")
    p.add_argument("--Lz", type=float, required=True, help="Original box length in z (nm)")

    p.add_argument("--repx", type=int, default=3, help="Replication factor in x (default: 3)")
    p.add_argument("--repy", type=int, default=2, help="Replication factor in y (default: 2)")

    p.add_argument("--aver", choices=["x", "y", "z"], required=True,
                   help="Axis to average over")
    p.add_argument("--sel", nargs="+", required=True,
                   help="Atom names/elements to include, e.g. O or O H")

    p.add_argument("--bin", type=float, default=0.02,
                   help="Bin size in nm for plotted axes (default: 0.02 nm)")
    p.add_argument("--nx", type=int, default=None,
                   help="Override number of bins along first plotted axis")
    p.add_argument("--ny", type=int, default=None,
                   help="Override number of bins along second plotted axis")

    p.add_argument("--skip", type=int, default=0,
                   help="Number of initial frames to skip in EACH file")
    p.add_argument("--stride", type=int, default=1,
                   help="Use every Nth frame after skip in EACH file")

    p.add_argument("--unit", choices=["nm-3", "kg/m3", "g/cm3"], default="nm-3",
                   help="Density unit")
    p.add_argument("--molar-mass", type=float, default=None,
                   help="Molar mass in g/mol, needed for kg/m3 or g/cm3")
    p.add_argument("--cmap", default="viridis", help="Matplotlib colormap")
    p.add_argument("--vmin", type=float, default=None, help="Minimum color scale")
    p.add_argument("--vmax", type=float, default=None, help="Maximum color scale")
    p.add_argument("--out", default="densmap.png", help="Output PNG filename")
    args = p.parse_args()

    if args.unit in ["kg/m3", "g/cm3"] and args.molar_mass is None:
        raise ValueError("--molar-mass is required for kg/m3 or g/cm3")
    if args.bin <= 0:
        raise ValueError("--bin must be > 0")

    selected = set(args.sel)
    angstrom_to_nm = 0.1

    Lx_eff = args.Lx * args.repx
    Ly_eff = args.Ly * args.repy
    Lz_eff = args.Lz

    if args.aver == "x":
        L1, L2, Lavg = Ly_eff, Lz_eff, Lx_eff
        idx1, idx2 = 1, 2
        xlabel, ylabel = "y (nm)", "z (nm)"
    elif args.aver == "y":
        L1, L2, Lavg = Lx_eff, Lz_eff, Ly_eff
        idx1, idx2 = 0, 2
        xlabel, ylabel = "x (nm)", "z (nm)"
    elif args.aver == "z":
        L1, L2, Lavg = Lx_eff, Ly_eff, Lz_eff
        idx1, idx2 = 0, 1
        xlabel, ylabel = "x (nm)", "y (nm)"
    else:
        raise ValueError("Invalid axis choice")

    nx_bins = args.nx if args.nx is not None else max(1, int(np.ceil(L1 / args.bin)))
    ny_bins = args.ny if args.ny is not None else max(1, int(np.ceil(L2 / args.bin)))

    dx = L1 / nx_bins
    dy = L2 / ny_bins
    bin_volume = dx * dy * Lavg

    hist_sum = np.zeros((nx_bins, ny_bins), dtype=float)
    nframes_used = 0

    for xyzfile in args.xyzfiles:
        frames = read_xyz_trajectory(xyzfile)
        frames = frames[args.skip::args.stride]

        if len(frames) == 0:
            print(f"Warning: no frames left after skip/stride in {xyzfile}")
            continue

        print(f"Processing {xyzfile}: using {len(frames)} frames")

        for symbols, coords_ang in frames:
            coords = coords_ang * angstrom_to_nm

            coords[:, 0] = np.mod(coords[:, 0], args.Lx)
            coords[:, 1] = np.mod(coords[:, 1], args.Ly)
            coords[:, 2] = np.mod(coords[:, 2], args.Lz)

            rep_symbols, rep_coords = replicate_xy(
                symbols, coords, args.Lx, args.Ly,
                nx_rep=args.repx, ny_rep=args.repy
            )

            mask = np.array([s in selected for s in rep_symbols], dtype=bool)
            sel_coords = rep_coords[mask]

            if sel_coords.size == 0:
                continue

            xy = sel_coords[:, [idx1, idx2]]

            H, _, _ = np.histogram2d(
                xy[:, 0], xy[:, 1],
                bins=[nx_bins, ny_bins],
                range=[[0.0, L1], [0.0, L2]]
            )

            hist_sum += H
            nframes_used += 1

    if nframes_used == 0:
        raise ValueError("No usable frames found across all input files.")

    density_nm3 = hist_sum / (nframes_used * bin_volume)

    rho, label = convert_units(
        density_nm3,
        args.unit,
        args.molar_mass if args.molar_mass is not None else 1.0
    )

    xcenters = np.linspace(0, L1, nx_bins, endpoint=False) + 0.5 * dx
    ycenters = np.linspace(0, L2, ny_bins, endpoint=False) + 0.5 * dy

    plt.figure(figsize=(8, 5))
    plt.imshow(
        rho.T,
        origin="lower",
        aspect="auto",
        extent=[xcenters.min(), xcenters.max(), ycenters.min(), ycenters.max()],
        cmap=args.cmap,
        vmin=args.vmin,
        vmax=args.vmax
    )
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.colorbar(label=label)
    plt.title(
        f"2D density map from XYZ ({len(args.xyzfiles)} files, "
        f"{nframes_used} frames total, averaged over {args.aver})"
    )
    plt.tight_layout()
    plt.savefig(args.out, dpi=200)
    plt.show()


if __name__ == "__main__":
    main()
