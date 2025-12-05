import numpy as np
import math
import argparse


def read_xyz(path):
    """Read XYZ file and return coordinates as an Nx3 numpy array."""
    with open(path) as f:
        nat = int(f.readline().strip())
        _comment = f.readline()
        coords = []
        for line in f:
            parts = line.split()
            if len(parts) >= 4:
                coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return np.array(coords)


def compute_layer_avg_z(coords, layer_size):
    """Compute average z-coordinate for each layer."""
    nat = coords.shape[0]
    if nat % layer_size != 0:
        raise ValueError(f"Atom count {nat} not divisible by layer size {layer_size}")

    nlayer = nat // layer_size
    avg_z_list = []

    for L in range(nlayer):
        sl = slice(L * layer_size, (L + 1) * layer_size)
        avg_z = np.mean(coords[sl, 2])
        print(f"{L} {avg_z}")
        avg_z_list.append(avg_z)

    return avg_z_list


def compute_rmsd(coords1, coords2, layer_size):
    """Compute per-layer RMSD and overall RMSD."""
    assert coords1.shape == coords2.shape, "Structures must have same number of atoms"
    nat = coords1.shape[0]

    if nat % layer_size != 0:
        raise ValueError(f"Atom count {nat} not divisible by layer size {layer_size}")

    nlayer = nat // layer_size
    print(f"Total atoms: {nat}, Layer size: {layer_size}, Total layers: {nlayer}\n")

    # Compute average z-coordinates of layers (from coords1)
    avg_z = compute_layer_avg_z(coords1, layer_size)

    for L in range(nlayer):
        sl = slice(L * layer_size, (L + 1) * layer_size)
        diff = coords1[sl] - coords2[sl]
        rmsd = math.sqrt(np.mean(np.sum(diff**2, axis=1)))
        rmsd_xyz = np.sqrt(np.mean(diff**2, axis=0))

        print(
            f"Layer {L+1:2d}: RMSD = {rmsd:.5f} Å  "
            f"(x={rmsd_xyz[0]:.5e}, y={rmsd_xyz[1]:.5e}, z={rmsd_xyz[2]:.5e})  "
            f" | avg(z) = {avg_z[L]:.5f}"
        )

    overall = math.sqrt(np.mean(np.sum((coords1 - coords2)**2, axis=1)))
    print(f"\nOverall RMSD: {overall:.5f} Å")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Layer-based RMSD calculator for XYZ files")
    parser.add_argument("file1", type=str, help="First XYZ file")
    parser.add_argument("file2", type=str, help="Second XYZ file")
    parser.add_argument("-l", "--layer", type=int, default=8, help="Layer size (# atoms per layer)")
    args = parser.parse_args()

    c1 = read_xyz(args.file1)
    c2 = read_xyz(args.file2)

    compute_layer_avg_z(c1, args.layer)
    compute_layer_avg_z(c2, args.layer)


