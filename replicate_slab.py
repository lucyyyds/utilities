#!/usr/bin/env python3
"""
replicate_and_interleave.py

Read a .gro or .g96 coordinate file, tile it nx × ny times in x–y,
then for As and In separately:
  • group into z‐layers,
  • sort each layer by x ascending,
  • sort layers by z ascending,
and finally interleave: As‐layer1, In‐layer1, As‐layer2, In‐layer2, …

Write back out in the same format (.gro → .gro, .g96 → .g96),
renumbering atoms sequentially and scaling x/y box vectors.

Usage:
    python replicate_and_interleave.py input.{gro,g96} output.{gro,g96} [nx] [ny]

Defaults: nx = ny = 4
"""
import sys, os
from collections import defaultdict

def parse_gro(path):
    with open(path) as f:
        title = f.readline().rstrip("\n")
        n_atoms = int(f.readline().split()[0])
        atom_lines = [f.readline() for _ in range(n_atoms)]
        box = list(map(float, f.readline().split()))
    atoms = []
    for line in atom_lines:
        resnr   = int(line[0:5])
        resname = line[5:10].strip()
        atomnm  = line[10:15].strip()
        x = float(line[20:28])
        y = float(line[28:36])
        z = float(line[36:44])
        vel = None
        if len(line) >= 68:
            vx = float(line[44:52]); vy = float(line[52:60]); vz = float(line[60:68])
            vel = (vx, vy, vz)
        atoms.append({
            'resnr': resnr,
            'resname': resname,
            'atomnm': atomnm,
            'x': x, 'y': y, 'z': z,
            'vel': vel
        })
    return title, box, atoms

def parse_g96(path):
    with open(path) as f:
        lines = [l.rstrip("\n") for l in f]
    i = 0
    title = ""
    # TITLE block
    if lines[i].strip() == "TITLE":
        i += 1
        title = lines[i]
        i += 2  # skip the "END" after title
    # POSITION block
    while lines[i].strip() != "POSITION":
        i += 1
    i += 1
    atoms = []
    while lines[i].strip() != "END":
        parts = lines[i].split()
        resnr   = int(parts[0])
        resname = parts[1]
        atomnm  = parts[2]
        x, y, z = map(float, parts[4:7])
        atoms.append({
            'resnr': resnr,
            'resname': resname,
            'atomnm': atomnm,
            'x': x, 'y': y, 'z': z,
            'vel': None
        })
        i += 1
    # BOX block
    while lines[i].strip() != "BOX":
        i += 1
    i += 1
    box = list(map(float, lines[i].split()))
    return title, box, atoms

def parse_input(path):
    ext = os.path.splitext(path)[1].lower()
    if ext == ".gro":
        return parse_gro(path), ".gro"
    elif ext == ".g96":
        return parse_g96(path), ".g96"
    else:
        raise ValueError(f"Unsupported extension '{ext}'")

def replicate_atoms(atoms, box, nx, ny):
    out = []
    for iy in range(ny):
        for ix in range(nx):
            dx = ix * box[0]
            dy = iy * box[1]
            for a in atoms:
                out.append({ **a,
                             'x': a['x'] + dx,
                             'y': a['y'] + dy })
    return out

def interleave_layers(atoms, z_tol=0.01, x_tol=0.01):
    # Group As and In into z-layers within tolerance z_tol
    layers = {'AS': defaultdict(list), 'IN': defaultdict(list)}
    for a in atoms:
        elem = a['atomnm'].upper()
        if elem in layers:
            zkey = round(a['z'] / z_tol) * z_tol
            layers[elem][zkey].append(a)
    # Within each z-layer, group by x within tolerance x_tol
    for elem in layers:
        for z in list(layers[elem].keys()):
            x_groups = defaultdict(list)
            for a in layers[elem][z]:
                xkey = round(a['x'] / x_tol) * x_tol
                x_groups[xkey].append(a)
            # flatten groups sorted by xkey, then atoms by ascending y
            sorted_atoms = []
            for xkey in sorted(x_groups):
                sorted_atoms.extend(sorted(x_groups[xkey], key=lambda a: a['y']))
            layers[elem][z] = sorted_atoms
    # Sorted z-values
    zs_AS = sorted(layers['AS'])
    zs_IN = sorted(layers['IN'])
    ordered = []
    # Interleave: As-layer1, In-layer1, ...
    for i in range(max(len(zs_AS), len(zs_IN))):
        if i < len(zs_AS): ordered.extend(layers['AS'][zs_AS[i]])
        if i < len(zs_IN): ordered.extend(layers['IN'][zs_IN[i]])
    return ordered

def write_gro(path, title, atoms, box, nx, ny):
    with open(path, 'w') as f:
        f.write(title + "\n")
        f.write(f"{len(atoms)}\n")
        for idx, a in enumerate(atoms, start=1):
            f.write(f"{a['resnr']:5d}{a['resname']:>5s}{a['atomnm']:>5s}"
                    f"{idx:10d}"
                    f"{a['x']:17.9f}{a['y']:17.9f}{a['z']:17.9f}\n")
        f.write(f"{box[0]*nx:16.9f}{box[1]*ny:16.9f}{box[2]:16.9f}\n")

def write_g96(path, title, atoms, box, nx, ny):
    with open(path, 'w') as f:
        f.write("TITLE\n")
        f.write(f"{title}\n")
        f.write("END\n")
        f.write("POSITION\n")
        for idx, a in enumerate(atoms, start=1):
            # exactly matching your example’s field widths & decimals:
            f.write(
                f"    {a['resnr']:d} {a['resname']:5s} {a['atomnm']:3s} {idx:8d} "           
                f"{a['x']:14.9f} {a['y']:14.9f} {a['z']:14.9f}\n"     
            )
        f.write("END\n")
        f.write("BOX\n")
        # three 16-wide, 9-dec fields exactly like your example
        f.write(
            f"{box[0]*nx:15.9f}"
            f"{box[1]*ny:15.9f}"
            f"{box[2]:15.9f}\n"
        )
        f.write("END\n")

def main():
    if len(sys.argv) < 3:
        print("Usage: replicate_and_interleave.py input.{gro,g96} output.{gro,g96} [nx] [ny]")
        sys.exit(1)

    inp, outp = sys.argv[1], sys.argv[2]
    nx = int(sys.argv[3]) if len(sys.argv) > 3 else 4
    ny = int(sys.argv[4]) if len(sys.argv) > 4 else 4

    (title, box, atoms), fmt = parse_input(inp)
    reps    = replicate_atoms(atoms, box, nx, ny)
    ordered = interleave_layers(reps)

    if fmt == ".gro":
        write_gro(outp, title, ordered, box, nx, ny)
    else:
        write_g96(outp, title, ordered, box, nx, ny)

    print(f"Wrote {fmt} → {outp} ({nx}×{ny}, {len(ordered)} atoms)")

if __name__ == "__main__":
    main()

