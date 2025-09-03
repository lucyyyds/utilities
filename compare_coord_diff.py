#!/usr/bin/env python3
"""
compare_struct.py

Compare two GROMACS structure files (.gro or .g96) atom by atom, reporting coordinate
and box-vector differences strictly greater than a given tolerance (default 0.001 nm).

Usage:
    python compare_struct.py file1.gro file2.gro [--tol 0.001]
    python compare_struct.py file1.g96 file2.g96 [--tol 0.001]
    python compare_struct.py file1.gro file2.g96 [--tol 0.001]

Outputs a summary of any atoms whose x, y, or z differ by > tol, in the format:

    Atom   12: Δx=0.002 Δy=0.000 Δz=0.005

And box-vector differences if any component differs by > tol.
"""
import sys
import os
import argparse

def parse_gro(path):
    """Parse a .gro file into (title, atoms, box)."""
    with open(path, 'r') as f:
        title = f.readline().rstrip('\n')
        natoms = int(f.readline().split()[0])
        atoms = []
        for _ in range(natoms):
            line = f.readline()
            resnr   = int(line[0:5])
            resname = line[5:10].strip()
            atomname= line[10:15].strip()
            atomnr  = int(line[15:20])
            x = float(line[20:28]); y = float(line[28:36]); z = float(line[36:44])
            vel = None
            if len(line) >= 68:
                vx = float(line[44:52]); vy = float(line[52:60]); vz = float(line[60:68])
                vel = (vx, vy, vz)
            atoms.append({
                'resnr': resnr, 'resname': resname,
                'atomname': atomname, 'atomnr': atomnr,
                'x': x, 'y': y, 'z': z, 'vel': vel
            })
        box = list(map(float, f.readline().split()))
    return title, atoms, box


def parse_g96(path):
    """Parse a .g96 file into (title, atoms, box)."""
    with open(path, 'r') as f:
        lines = f.readlines()

    # Initialize
    title = None
    atoms = []
    box = []
    i = 0
    n = len(lines)

    # Parse TITLE block if present
    if lines[i].strip().upper() == 'TITLE':
        i += 1
        title_lines = []
        while i < n and lines[i].strip().upper() != 'END':
            title_lines.append(lines[i].rstrip('\n'))
            i += 1
        title = '\n'.join(title_lines)
    else:
        title = lines[i].rstrip('\n')

    # Advance to POSITION block
    while i < n and lines[i].strip().upper() != 'POSITION':
        i += 1
    # Skip 'POSITION' line
    if i < n and lines[i].strip().upper() == 'POSITION':
        i += 1

    # Parse atoms until next 'END'
    while i < n and lines[i].strip().upper() != 'END':
        line = lines[i].strip()
        i += 1
        if not line:
            continue
        tokens = line.split()
        try:
            x = float(tokens[-3]); y = float(tokens[-2]); z = float(tokens[-1])
        except (ValueError, IndexError):
            continue
        # Optionally capture residue/atom info
        resnr = None; resname = None; atomname = None; atomnr = None
        if len(tokens) >= 6:
            try:
                resnr = int(tokens[0])
            except ValueError:
                resnr = None
            resname = tokens[1]
            atomname = tokens[2]
            try:
                atomnr = int(tokens[3])
            except ValueError:
                atomnr = None
        atoms.append({
            'resnr': resnr, 'resname': resname,
            'atomname': atomname, 'atomnr': atomnr,
            'x': x, 'y': y, 'z': z, 'vel': None
        })

    # Advance to BOX block
    while i < n and lines[i].strip().upper() != 'BOX':
        i += 1
    # Skip 'BOX' line
    if i < n and lines[i].strip().upper() == 'BOX':
        i += 1
    # Parse box vectors until next 'END'
    while i < n and lines[i].strip().upper() != 'END':
        toks = lines[i].split()
        floats = []
        for tok in toks:
            try:
                floats.append(float(tok))
            except ValueError:
                pass
        if len(floats) >= 3:
            box = floats[:3]
            break
        i += 1

    return title, atoms, box


def parse_structure(path):
    """Dispatch to the correct parser based on file extension."""
    ext = os.path.splitext(path)[1].lower()
    if ext == '.gro':
        return parse_gro(path)
    elif ext == '.g96':
        return parse_g96(path)
    else:
        raise ValueError(f"Unsupported file format '{ext}' for file {path}")


def compare_structures(f1, f2, tol):
    t1, a1, box1 = parse_structure(f1)
    t2, a2, box2 = parse_structure(f2)

    if len(a1) != len(a2):
        print(f"Atom count differs: {len(a1)} vs {len(a2)}")
        sys.exit(1)

    print(f"Comparing {f1} vs {f2} with tolerance {tol:.3f} nm (reporting only Δ > tol)\n")

    diffs = []
    for i, (u, v) in enumerate(zip(a1, a2), start=1):
        dx = abs(u['x'] - v['x'])
        dy = abs(u['y'] - v['y'])
        dz = abs(u['z'] - v['z'])
        if dx > tol or dy > tol or dz > tol:
            diffs.append((i, dx, dy, dz))

    if diffs:
        print("Atoms with coordinate differences > tol:\n")
        for idx, dx, dy, dz in diffs:
            print(f"Atom {idx:5d}: Δx={dx:.3f} Δy={dy:.3f} Δz={dz:.3f}")
    else:
        print("No atom coordinate differences > tol.\n")

    # Compare box vectors
    bz = []
    for comp, b1, b2 in zip('xyz', box1, box2):
        d = abs(b1 - b2)
        if d > tol:
            bz.append((comp, b1, b2, d))

    if bz:
        print("\nBox-vector differences > tol:")
        for comp, b1, b2, d in bz:
            print(f"  box[{comp}]: {b1:.5f} vs {b2:.5f} (Δ={d:.5f})")
    else:
        print("No box-vector differences > tol.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compare two GROMACS structure files (.gro or .g96)."
    )
    parser.add_argument("file1", help="first structure file (.gro or .g96)")
    parser.add_argument("file2", help="second structure file (.gro or .g96)")
    parser.add_argument("--tol", type=float, default=0.001,
                        help="tolerance in nm (default: 0.001); differences ≤ tol are ignored")
    args = parser.parse_args()

    compare_structures(args.file1, args.file2, args.tol)

