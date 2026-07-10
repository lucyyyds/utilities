#!/usr/bin/env python3
"""
Reduced normal-mode analysis for a constrained GROMACS Hessian.

Inputs
------
1. GROMACS Hessian (.mtx), e.g. nm.mtx
2. GROMACS index file (.ndx), containing the frozen atom group
3. GROMACS structure file (.g96), providing atom names/order

The script:
- reads atom names directly from the G96 POSITION/POSITIONRED section,
- infers chemical elements and atomic masses,
- removes frozen coordinates from the Hessian,
- mass-weights and diagonalizes the reduced Hessian,
- generates:
    <prefix>_vibrational_modes.csv
    <prefix>_broadened_vdos.csv
    <prefix>_harmonic_vdos.png

Example
-------
python gromacs_reduced_vib_g96.py \
    nm.mtx \
    --index index.ndx \
    --structure InAs110-6-4-2.g96 \
    --freeze-group FreezeAtoms \
    --prefix gromacs_reduced \
    --sigma 3 \
    --xmax 300
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Standard atomic masses in unified atomic mass units (u).
ATOMIC_MASSES = {
    "H": 1.008,
    "He": 4.002602,
    "Li": 6.94,
    "Be": 9.0121831,
    "B": 10.81,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "F": 18.998403163,
    "Ne": 20.1797,
    "Na": 22.98976928,
    "Mg": 24.305,
    "Al": 26.9815385,
    "Si": 28.085,
    "P": 30.973761998,
    "S": 32.06,
    "Cl": 35.45,
    "Ar": 39.948,
    "K": 39.0983,
    "Ca": 40.078,
    "Sc": 44.955908,
    "Ti": 47.867,
    "V": 50.9415,
    "Cr": 51.9961,
    "Mn": 54.938044,
    "Fe": 55.845,
    "Co": 58.933194,
    "Ni": 58.6934,
    "Cu": 63.546,
    "Zn": 65.38,
    "Ga": 69.723,
    "Ge": 72.630,
    "As": 74.921595,
    "Se": 78.971,
    "Br": 79.904,
    "Kr": 83.798,
    "Rb": 85.4678,
    "Sr": 87.62,
    "Y": 88.90584,
    "Zr": 91.224,
    "Nb": 92.90637,
    "Mo": 95.95,
    "Ru": 101.07,
    "Rh": 102.90550,
    "Pd": 106.42,
    "Ag": 107.8682,
    "Cd": 112.414,
    "In": 114.818,
    "Sn": 118.710,
    "Sb": 121.760,
    "Te": 127.60,
    "I": 126.90447,
    "Xe": 131.293,
    "Cs": 132.90545196,
    "Ba": 137.327,
    "La": 138.90547,
    "Ce": 140.116,
    "Pt": 195.084,
    "Au": 196.96657,
    "Hg": 200.592,
    "Pb": 207.2,
}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract the mobile-coordinate block of a GROMACS Hessian "
            "and perform mass-weighted normal-mode analysis."
        )
    )

    parser.add_argument(
        "mtx_file",
        type=Path,
        help="GROMACS Hessian file, typically nm.mtx",
    )

    parser.add_argument(
        "--index",
        type=Path,
        required=True,
        help="GROMACS index file (.ndx)",
    )

    parser.add_argument(
        "--structure",
        type=Path,
        required=True,
        help="GROMACS structure file in .g96 format",
    )

    parser.add_argument(
        "--freeze-group",
        default="FreezeAtoms",
        help="Frozen atom group in the index file (default: FreezeAtoms)",
    )

    parser.add_argument(
        "--system-group",
        default="System",
        help="Full-system group in the index file (default: System)",
    )

    parser.add_argument(
        "--prefix",
        default=None,
        help="Output prefix (default: <mtx stem>_reduced)",
    )

    parser.add_argument(
        "--sigma",
        type=float,
        default=3.0,
        help="Gaussian broadening sigma in cm^-1 (default: 3.0)",
    )

    parser.add_argument(
        "--xmin",
        type=float,
        default=0.0,
        help="Minimum plotted wavenumber in cm^-1 (default: 0)",
    )

    parser.add_argument(
        "--xmax",
        type=float,
        default=300.0,
        help="Maximum plotted wavenumber in cm^-1 (default: 300)",
    )

    parser.add_argument(
        "--step",
        type=float,
        default=0.05,
        help="Wavenumber grid spacing in cm^-1 (default: 0.05)",
    )

    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="PNG resolution (default: 300)",
    )

    parser.add_argument(
        "--include-negative",
        action="store_true",
        help="Include negative frequencies in the broadened VDOS",
    )

    return parser.parse_args()


def parse_index_file(index_file: Path) -> dict[str, list[int]]:
    """Parse a standard GROMACS .ndx file."""
    groups: dict[str, list[int]] = {}
    current_group: str | None = None

    for raw_line in index_file.read_text(errors="replace").splitlines():
        line = raw_line.strip()

        if line.startswith("[") and line.endswith("]"):
            current_group = line[1:-1].strip()
            groups[current_group] = []
        elif line and current_group is not None:
            groups[current_group].extend(int(value) for value in line.split())

    return groups


def parse_g96_atom_names(g96_file: Path) -> list[str]:
    """
    Read atom names from a GROMACS G96 POSITION or POSITIONRED section.

    Standard record format:
        residue_number residue_name atom_name atom_number x y z
    """
    lines = g96_file.read_text(errors="replace").splitlines()

    atom_names: list[str] = []
    in_position = False
    found_position = False

    for raw_line in lines:
        line = raw_line.strip()

        if line in {"POSITION", "POSITIONRED"}:
            in_position = True
            found_position = True
            continue

        if in_position and line == "END":
            break

        if not in_position or not line:
            continue

        fields = line.split()

        if len(fields) < 7:
            raise ValueError(
                f"Invalid G96 atom record in POSITION section:\n{raw_line}"
            )

        atom_names.append(fields[2])

    if not found_position:
        raise ValueError(
            f"No POSITION or POSITIONRED section found in {g96_file}."
        )

    if not atom_names:
        raise ValueError(
            f"No atom records found in the G96 POSITION section of {g96_file}."
        )

    return atom_names


def infer_element(atom_name: str) -> str:
    """
    Infer a chemical element from a GROMACS atom name.

    Examples
    --------
    As   -> As
    As1  -> As
    IN   -> In
    In12 -> In
    OW   -> O
    HW1  -> H
    """
    letters = re.sub(r"[^A-Za-z]", "", atom_name)

    if not letters:
        raise ValueError(
            f"Cannot infer a chemical element from atom name '{atom_name}'."
        )

    # First try a valid two-letter element symbol.
    if len(letters) >= 2:
        candidate = letters[0].upper() + letters[1].lower()
        if candidate in ATOMIC_MASSES:
            return candidate

    # Then try a one-letter element.
    candidate = letters[0].upper()
    if candidate in ATOMIC_MASSES:
        return candidate

    raise ValueError(
        f"Unknown element for atom name '{atom_name}'. "
        "Add the element and mass to ATOMIC_MASSES."
    )


def read_dense_gromacs_hessian(
    mtx_file: Path,
    n_coordinates: int,
) -> np.ndarray:
    """
    Read the dense double-precision Hessian from a GROMACS .mtx file.

    The dense matrix is assumed to occupy the final:
        n_coordinates * n_coordinates * 8 bytes

    The validated file format here is big-endian float64, matching the
    double-precision GROMACS Hessian used for this workflow.
    """
    matrix_count = n_coordinates * n_coordinates
    matrix_nbytes = matrix_count * 8
    file_size = mtx_file.stat().st_size
    offset = file_size - matrix_nbytes

    if offset < 0:
        raise ValueError(
            "The .mtx file is too small for the expected Hessian dimension "
            f"{n_coordinates} x {n_coordinates}."
        )

    data = np.fromfile(
        mtx_file,
        dtype=">f8",
        offset=offset,
        count=matrix_count,
    )

    if data.size != matrix_count:
        raise ValueError(
            f"Expected {matrix_count} Hessian values but read {data.size}."
        )

    matrix = data.reshape(n_coordinates, n_coordinates)

    if not np.all(np.isfinite(matrix)):
        raise ValueError(
            "Non-finite Hessian values were read. "
            "The .mtx precision or byte order may differ."
        )

    return matrix


def calculate_element_fractions(
    mobile_elements: np.ndarray,
    mode_vectors: np.ndarray,
) -> tuple[list[str], dict[str, np.ndarray]]:
    """
    Calculate qualitative eigenvector participation fractions by element.

    These are not IR, Raman, or neutron-scattering intensities.
    """
    unique_elements = list(dict.fromkeys(mobile_elements.tolist()))
    total_weight = np.sum(mode_vectors**2, axis=(1, 2))

    fractions: dict[str, np.ndarray] = {}

    for element in unique_elements:
        mask = mobile_elements == element
        element_weight = np.sum(
            mode_vectors[:, mask, :] ** 2,
            axis=(1, 2),
        )

        fractions[element] = np.divide(
            element_weight,
            total_weight,
            out=np.zeros_like(element_weight),
            where=total_weight > 0,
        )

    return unique_elements, fractions


def gaussian_broaden(
    frequencies: np.ndarray,
    weights: np.ndarray,
    wavenumbers: np.ndarray,
    sigma: float,
) -> np.ndarray:
    """Gaussian broaden discrete normal-mode frequencies."""
    delta = wavenumbers[None, :] - frequencies[:, None]
    gaussians = np.exp(-0.5 * (delta / sigma) ** 2)
    return np.sum(weights[:, None] * gaussians, axis=0)


def main() -> int:
    args = parse_arguments()

    for path in (args.mtx_file, args.index, args.structure):
        if not path.is_file():
            print(f"Error: file not found: {path}", file=sys.stderr)
            return 1

    if args.structure.suffix.lower() != ".g96":
        print(
            "Error: --structure must be a .g96 file.",
            file=sys.stderr,
        )
        return 1

    if args.sigma <= 0:
        print("Error: --sigma must be positive.", file=sys.stderr)
        return 1

    if args.step <= 0:
        print("Error: --step must be positive.", file=sys.stderr)
        return 1

    if args.xmax <= args.xmin:
        print(
            "Error: --xmax must be greater than --xmin.",
            file=sys.stderr,
        )
        return 1

    prefix = args.prefix or f"{args.mtx_file.stem}_reduced"

    try:
        # --------------------------------------------------------
        # Read atom groups
        # --------------------------------------------------------
        groups = parse_index_file(args.index)

        if args.system_group not in groups:
            raise ValueError(
                f"System group '{args.system_group}' not found in {args.index}."
            )

        if args.freeze_group not in groups:
            raise ValueError(
                f"Freeze group '{args.freeze_group}' not found in {args.index}."
            )

        all_atoms = groups[args.system_group]
        frozen_atoms = set(groups[args.freeze_group])
        mobile_atoms = [
            atom for atom in all_atoms
            if atom not in frozen_atoms
        ]

        if not mobile_atoms:
            raise ValueError(
                "No mobile atoms remain after applying the frozen atom group."
            )

        # --------------------------------------------------------
        # Read atom names/elements/masses from G96
        # --------------------------------------------------------
        atom_names = parse_g96_atom_names(args.structure)

        if len(atom_names) != len(all_atoms):
            raise ValueError(
                "Atom count mismatch: "
                f"G96 file has {len(atom_names)} atoms, "
                f"but index group '{args.system_group}' has "
                f"{len(all_atoms)} atoms."
            )

        elements = [infer_element(name) for name in atom_names]
        masses = np.array(
            [ATOMIC_MASSES[element] for element in elements],
            dtype=float,
        )

        # --------------------------------------------------------
        # Read full Hessian
        # --------------------------------------------------------
        n_atoms = len(all_atoms)
        n_coordinates = 3 * n_atoms

        full_hessian = read_dense_gromacs_hessian(
            args.mtx_file,
            n_coordinates,
        )

        # --------------------------------------------------------
        # Extract mobile-coordinate block
        # --------------------------------------------------------
        coordinate_indices = np.array(
            [
                3 * (atom - 1) + direction
                for atom in mobile_atoms
                for direction in range(3)
            ],
            dtype=int,
        )

        reduced_hessian = full_hessian[
            np.ix_(coordinate_indices, coordinate_indices)
        ]

        # Remove tiny numerical asymmetry.
        reduced_hessian = 0.5 * (
            reduced_hessian + reduced_hessian.T
        )

        # --------------------------------------------------------
        # Mass weighting
        # --------------------------------------------------------
        mobile_indices = np.asarray(mobile_atoms, dtype=int) - 1
        mobile_elements = np.asarray(elements)[mobile_indices]
        mobile_atom_masses = masses[mobile_indices]
        coordinate_masses = np.repeat(mobile_atom_masses, 3)

        mass_weighted_hessian = reduced_hessian / np.sqrt(
            coordinate_masses[:, None]
            * coordinate_masses[None, :]
        )

        # --------------------------------------------------------
        # Diagonalize
        # --------------------------------------------------------
        eigenvalues, eigenvector_columns = np.linalg.eigh(
            mass_weighted_hessian
        )

        # Conversion validated against the GROMACS nmeig workflow:
        # eigenvalue -> frequency in cm^-1
        conversion = 1.0e12 / (
            2.0 * np.pi * 2.99792458e10
        )

        frequencies = (
            np.sign(eigenvalues)
            * np.sqrt(np.abs(eigenvalues))
            * conversion
        )

        n_modes = frequencies.size
        n_mobile_atoms = len(mobile_atoms)

        mode_vectors = eigenvector_columns.T.reshape(
            n_modes,
            n_mobile_atoms,
            3,
        )

        unique_elements, element_fractions = (
            calculate_element_fractions(
                mobile_elements,
                mode_vectors,
            )
        )

    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    # ------------------------------------------------------------
    # Output 1: mode table
    # ------------------------------------------------------------
    mode_data: dict[str, np.ndarray] = {
        "Mode": np.arange(1, n_modes + 1),
        "Frequency_cm-1": frequencies,
    }

    for element in unique_elements:
        mode_data[f"{element}_eigenvector_fraction"] = (
            element_fractions[element]
        )

    mode_df = pd.DataFrame(mode_data)
    mode_csv = Path(f"{prefix}_vibrational_modes.csv")
    mode_df.to_csv(mode_csv, index=False)

    # ------------------------------------------------------------
    # Output 2: broadened VDOS
    # ------------------------------------------------------------
    if args.include_negative:
        use_mask = np.ones(n_modes, dtype=bool)
    else:
        use_mask = frequencies > 0

    selected_frequencies = frequencies[use_mask]

    x = np.arange(
        args.xmin,
        args.xmax + 0.5 * args.step,
        args.step,
    )

    total_raw = gaussian_broaden(
        selected_frequencies,
        np.ones(selected_frequencies.size),
        x,
        args.sigma,
    )

    projected_raw: dict[str, np.ndarray] = {}

    for element in unique_elements:
        projected_raw[element] = gaussian_broaden(
            selected_frequencies,
            element_fractions[element][use_mask],
            x,
            args.sigma,
        )

    scale = total_raw.max()

    if scale > 0:
        total_vdos = total_raw / scale
        projected_vdos = {
            element: values / scale
            for element, values in projected_raw.items()
        }
    else:
        total_vdos = total_raw
        projected_vdos = projected_raw

    vdos_data: dict[str, np.ndarray] = {
        "Wavenumber_cm-1": x,
        "Total_VDOS": total_vdos,
    }

    for element in unique_elements:
        vdos_data[f"{element}_projected_qualitative"] = (
            projected_vdos[element]
        )

    vdos_df = pd.DataFrame(vdos_data)
    vdos_csv = Path(f"{prefix}_broadened_vdos.csv")
    vdos_df.to_csv(vdos_csv, index=False)

    # ------------------------------------------------------------
    # Output 3: plot
    # ------------------------------------------------------------
    plt.figure(figsize=(8, 5))
    plt.plot(
        x,
        total_vdos,
        label="Total",
        linewidth=1.8,
    )

    for element in unique_elements:
        plt.plot(
            x,
            projected_vdos[element],
            label=f"{element} contribution",
            linewidth=1.4,
        )

    plt.xlabel(r"Wavenumber (cm$^{-1}$)")
    plt.ylabel("Relative intensity")
    plt.xlim(args.xmin, args.xmax)
    plt.ylim(bottom=0)
    plt.legend(frameon=False)
    plt.tight_layout()

    plot_png = Path(f"{prefix}_harmonic_vdos.png")
    plt.savefig(
        plot_png,
        dpi=args.dpi,
        bbox_inches="tight",
    )
    plt.close()

    # ------------------------------------------------------------
    # Terminal summary
    # ------------------------------------------------------------
    print()
    print("Reduced GROMACS vibrational analysis complete")
    print("---------------------------------------------")
    print(f"Hessian file        : {args.mtx_file}")
    print(f"Index file          : {args.index}")
    print(f"G96 structure       : {args.structure}")
    print(f"Total atoms         : {len(all_atoms)}")
    print(f"Frozen atoms        : {len(frozen_atoms)}")
    print(f"Mobile atoms        : {len(mobile_atoms)}")
    print(
        f"Reduced Hessian     : "
        f"{reduced_hessian.shape[0]} x "
        f"{reduced_hessian.shape[1]}"
    )
    print(f"Modes               : {n_modes}")
    print(
        f"Frequency range     : "
        f"{frequencies.min():.4f} to "
        f"{frequencies.max():.4f} cm^-1"
    )
    print(f"Negative modes      : {int(np.sum(frequencies < 0))}")
    print(
        f"|frequency| < 10    : "
        f"{int(np.sum(np.abs(frequencies) < 10))}"
    )
    print(f"Gaussian sigma      : {args.sigma:.3f} cm^-1")
    print(f"Elements            : {', '.join(unique_elements)}")
    print()
    print("Element masses used:")
    for element in unique_elements:
        print(f"  {element:>3s} : {ATOMIC_MASSES[element]:.8f} u")
    print()
    print("Created:")
    print(f"  {mode_csv}")
    print(f"  {vdos_csv}")
    print(f"  {plot_png}")
    print()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
