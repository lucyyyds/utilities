#!/usr/bin/env python3
"""
Post-process a CP2K VIBRATIONAL_ANALYSIS MOLDEN_VIB file.

Outputs
-------
<prefix>_vibrational_modes.csv
    Mode number, frequency, and qualitative element-projected
    eigenvector fractions.

<prefix>_broadened_vdos.csv
    Wavenumber grid, total broadened VDOS, and element-projected curves.

<prefix>_harmonic_vdos.png
    Publication-ready PNG plot of the total and element-projected VDOS.

Example
-------
python cp2k_vib_postprocess.py INAS110-6-4-2-VIBRATIONS-1.mol

python cp2k_vib_postprocess.py INAS110-6-4-2-VIBRATIONS-1.mol \
    --prefix inas110 --sigma 3 --xmax 300
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Post-process a CP2K VIBRATIONAL_ANALYSIS MOLDEN_VIB file "
            "into mode tables, broadened VDOS data, and a PNG plot."
        )
    )

    parser.add_argument(
        "input_file",
        type=Path,
        help="CP2K MOLDEN_VIB file, typically *-VIBRATIONS-*.mol",
    )

    parser.add_argument(
        "--prefix",
        type=str,
        default=None,
        help=(
            "Output filename prefix. Default: input filename stem with "
            "'-VIBRATIONS-*' removed."
        ),
    )

    parser.add_argument(
        "--sigma",
        type=float,
        default=3.0,
        help="Gaussian broadening standard deviation in cm^-1 (default: 3.0)",
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
        help="PNG resolution in dots per inch (default: 300)",
    )

    parser.add_argument(
        "--include-negative",
        action="store_true",
        help="Include negative frequencies in the broadened VDOS.",
    )

    return parser.parse_args()


def find_section(lines: list[str], section_name: str) -> int:
    """Return the line index of a Molden section header."""
    for i, line in enumerate(lines):
        if line.strip().upper() == section_name.upper():
            return i
    raise ValueError(f"Required section not found: {section_name}")


def parse_molden_vibrations(
    input_file: Path,
) -> tuple[list[str], np.ndarray, np.ndarray]:
    """
    Parse atom symbols, frequencies, and normal-mode eigenvectors.

    Returns
    -------
    elements
        Atomic symbols for all atoms in the Molden file.
    frequencies
        Vibrational frequencies in cm^-1.
    eigenvectors
        Array with shape (n_modes, n_atoms, 3).
    """
    lines = input_file.read_text(errors="replace").splitlines()

    atoms_idx = find_section(lines, "[Atoms] AU")
    freq_idx = find_section(lines, "[FREQ]")
    fr_coord_idx = find_section(lines, "[FR-COORD]")
    norm_coord_idx = find_section(lines, "[FR-NORM-COORD]")

    # Atom symbols
    atom_lines = lines[atoms_idx + 1 : freq_idx]
    elements: list[str] = []

    for line in atom_lines:
        fields = line.split()
        if not fields:
            continue
        elements.append(fields[0])

    if not elements:
        raise ValueError("No atoms found in [Atoms] AU section.")

    # Frequencies
    frequencies = np.array(
        [
            float(line.strip())
            for line in lines[freq_idx + 1 : fr_coord_idx]
            if line.strip()
        ],
        dtype=float,
    )

    if frequencies.size == 0:
        raise ValueError("No frequencies found in [FREQ] section.")

    # Normal-mode vectors
    n_atoms = len(elements)
    n_modes = len(frequencies)
    eigenvectors = np.zeros((n_modes, n_atoms, 3), dtype=float)

    i = norm_coord_idx + 1

    for mode in range(n_modes):
        while i < len(lines) and not lines[i].strip():
            i += 1

        if i >= len(lines) or not lines[i].strip().lower().startswith("vibration"):
            raise ValueError(
                f"Expected 'vibration' header for mode {mode + 1}, "
                f"but found: {lines[i] if i < len(lines) else 'EOF'}"
            )

        i += 1

        for atom in range(n_atoms):
            if i >= len(lines):
                raise ValueError(
                    f"Unexpected end of file while reading mode {mode + 1}."
                )

            fields = lines[i].split()
            if len(fields) < 3:
                raise ValueError(
                    f"Invalid eigenvector line at mode {mode + 1}, "
                    f"atom {atom + 1}: {lines[i]}"
                )

            eigenvectors[mode, atom] = [
                float(fields[0]),
                float(fields[1]),
                float(fields[2]),
            ]
            i += 1

    return elements, frequencies, eigenvectors


def default_prefix(input_file: Path) -> str:
    """Create a clean default prefix from the CP2K Molden filename."""
    stem = input_file.stem
    upper = stem.upper()

    marker = "-VIBRATIONS-"
    if marker in upper:
        stem = stem[: upper.index(marker)]

    return stem


def calculate_element_fractions(
    elements: list[str],
    eigenvectors: np.ndarray,
) -> tuple[list[str], dict[str, np.ndarray]]:
    """
    Calculate qualitative per-mode eigenvector fractions by element.

    The contribution for one element is:
        sum(|e|^2 for atoms of that element)
        -------------------------------------
        sum(|e|^2 for all atoms)

    These are qualitative eigenvector participation fractions, not
    IR intensities, Raman intensities, or neutron-weighted PDOS.
    """
    element_array = np.asarray(elements)
    unique_elements = list(dict.fromkeys(elements))

    total_weight = np.sum(eigenvectors**2, axis=(1, 2))

    fractions: dict[str, np.ndarray] = {}

    for element in unique_elements:
        mask = element_array == element
        weight = np.sum(eigenvectors[:, mask, :] ** 2, axis=(1, 2))

        fractions[element] = np.divide(
            weight,
            total_weight,
            out=np.zeros_like(weight),
            where=total_weight > 0,
        )

    return unique_elements, fractions


def gaussian_broaden(
    frequencies: np.ndarray,
    weights: np.ndarray,
    x: np.ndarray,
    sigma: float,
) -> np.ndarray:
    """Gaussian broaden discrete vibrational modes."""
    delta = x[None, :] - frequencies[:, None]
    gaussian = np.exp(-0.5 * (delta / sigma) ** 2)
    return np.sum(weights[:, None] * gaussian, axis=0)


def main() -> int:
    args = parse_arguments()

    if not args.input_file.is_file():
        print(f"Error: input file not found: {args.input_file}", file=sys.stderr)
        return 1

    if args.sigma <= 0:
        print("Error: --sigma must be positive.", file=sys.stderr)
        return 1

    if args.step <= 0:
        print("Error: --step must be positive.", file=sys.stderr)
        return 1

    if args.xmax <= args.xmin:
        print("Error: --xmax must be greater than --xmin.", file=sys.stderr)
        return 1

    prefix = args.prefix or default_prefix(args.input_file)

    try:
        elements, frequencies, eigenvectors = parse_molden_vibrations(
            args.input_file
        )
    except Exception as exc:
        print(f"Error while parsing {args.input_file}: {exc}", file=sys.stderr)
        return 1

    unique_elements, element_fractions = calculate_element_fractions(
        elements,
        eigenvectors,
    )

    n_atoms = len(elements)
    n_modes = len(frequencies)
    n_negative = int(np.sum(frequencies < 0))
    n_zero_like = int(np.sum(np.abs(frequencies) < 10))

    # ------------------------------------------------------------
    # 1. Mode table CSV
    # ------------------------------------------------------------
    mode_data: dict[str, np.ndarray] = {
        "Mode": np.arange(1, n_modes + 1),
        "Frequency_cm-1": frequencies,
    }

    for element in unique_elements:
        mode_data[f"{element}_eigenvector_fraction"] = element_fractions[element]

    mode_df = pd.DataFrame(mode_data)

    mode_csv = Path(f"{prefix}_vibrational_modes.csv")
    mode_df.to_csv(mode_csv, index=False)

    # ------------------------------------------------------------
    # 2. Broadened VDOS
    # ------------------------------------------------------------
    if args.include_negative:
        use_mask = np.ones(n_modes, dtype=bool)
    else:
        use_mask = frequencies > 0

    selected_freq = frequencies[use_mask]

    x = np.arange(
        args.xmin,
        args.xmax + 0.5 * args.step,
        args.step,
    )

    total_raw = gaussian_broaden(
        selected_freq,
        np.ones(selected_freq.size),
        x,
        args.sigma,
    )

    projected_raw: dict[str, np.ndarray] = {}

    for element in unique_elements:
        projected_raw[element] = gaussian_broaden(
            selected_freq,
            element_fractions[element][use_mask],
            x,
            args.sigma,
        )

    # Normalize all curves by the same total-VDOS maximum.
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

    spectrum_data: dict[str, np.ndarray] = {
        "Wavenumber_cm-1": x,
        "Total_VDOS": total_vdos,
    }

    for element in unique_elements:
        spectrum_data[f"{element}_projected_qualitative"] = projected_vdos[
            element
        ]

    spectrum_df = pd.DataFrame(spectrum_data)

    spectrum_csv = Path(f"{prefix}_broadened_vdos.csv")
    spectrum_df.to_csv(spectrum_csv, index=False)

    # ------------------------------------------------------------
    # 3. PNG plot
    # ------------------------------------------------------------
    plt.figure(figsize=(8, 5))
    plt.plot(x, total_vdos, label="Total", linewidth=1.8)

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
    plt.savefig(plot_png, dpi=args.dpi, bbox_inches="tight")
    plt.close()

    # ------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------
    print()
    print("CP2K vibrational analysis post-processing complete")
    print("-------------------------------------------------")
    print(f"Input file          : {args.input_file}")
    print(f"Atoms               : {n_atoms}")
    print(f"Modes               : {n_modes}")
    print(f"Frequency range     : {frequencies.min():.4f} to "
          f"{frequencies.max():.4f} cm^-1")
    print(f"Negative modes      : {n_negative}")
    print(f"|frequency| < 10    : {n_zero_like}")
    print(f"Gaussian sigma      : {args.sigma:.3f} cm^-1")
    print(f"Elements            : {', '.join(unique_elements)}")
    print()
    print("Created:")
    print(f"  {mode_csv}")
    print(f"  {spectrum_csv}")
    print(f"  {plot_png}")
    print()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
