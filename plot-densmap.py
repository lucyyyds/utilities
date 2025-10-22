#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib.pyplot as plt

def convert_units(data, unit, molar_mass):
    """Convert number density (molecules/nm^3) to desired unit."""
    if unit == "nm-3":
        return data, "Number density (molecules / nm³)"
    elif unit == "kg/m3":
        # kg/m³ conversion
        factor = (molar_mass*1e-3) / (6.022e23 * 1e-27)
        return data * factor, "Mass density (kg/m³)"
    elif unit == "g/cm3":
        # g/cm³ conversion
        factor = molar_mass / (6.022e23 * 1e-21)
        return data * factor, "Mass density (g/cm³)"
    else:
        raise ValueError("Unknown unit. Choose from: nm-3, kg/m3, g/cm3")

def main():
    p = argparse.ArgumentParser(
        description="Visualize GROMACS densmap.dat with units and axis choice."
    )
    p.add_argument("datfile", help="densmap.dat file from gmx densmap")
    p.add_argument("--Lx", type=float, required=True, help="Box length in x (nm)")
    p.add_argument("--Ly", type=float, required=True, help="Box length in y (nm)")
    p.add_argument("--Lz", type=float, required=True, help="Box length in z (nm)")
    p.add_argument("--aver", choices=["x", "y", "z"], required=True,
                   help="Axis that was averaged over in gmx densmap")
    p.add_argument("--unit", choices=["nm-3", "kg/m3", "g/cm3"], default="g/cm3",
                   help="Output unit for density (default: g/cm3)")
    p.add_argument("--molar-mass", type=float, required=True,
                   help="Molecular mass in g/mol (e.g., 18.015 for water)")
    p.add_argument("--out", default="densmap.png", help="Output PNG filename")
    args = p.parse_args()

    # Load densmap
    dm = np.loadtxt(args.datfile)

    # Shapes
    n1, n2 = dm.shape

    # Decide which axes are present
    if args.aver == "x":
        # averaged over x → map is (y,z)
        L1, L2 = args.Ly, args.Lz
        xlabel, ylabel = "y (nm)", "z (nm)"
    elif args.aver == "y":
        # averaged over y → map is (x,z)
        L1, L2 = args.Lx, args.Lz
        xlabel, ylabel = "x (nm)", "z (nm)"
    elif args.aver == "z":
        # averaged over z → map is (x,y)
        L1, L2 = args.Lx, args.Ly
        xlabel, ylabel = "x (nm)", "y (nm)"
    else:
        raise ValueError("Invalid axis")

    # Bin centers
    axis1 = np.linspace(0, L1, n1, endpoint=False) + 0.5 * L1 / n1
    axis2 = np.linspace(0, L2, n2, endpoint=False) + 0.5 * L2 / n2

    # Convert units
    rho, label = convert_units(dm, args.unit)

    # Plot
    plt.figure(figsize=(7,5))
    plt.imshow(rho.T, origin="lower", aspect="auto",
               extent=[axis1.min(), axis1.max(), axis2.min(), axis2.max()])
    plt.xlabel(xlabel); plt.ylabel(ylabel)
    plt.colorbar(label=label)
    plt.title(f"Densmap averaged over {args.aver.upper()}")
    plt.tight_layout()
    plt.savefig(args.out, dpi=200)
    plt.show()

if __name__ == "__main__":
    main()
