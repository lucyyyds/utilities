import MDAnalysis as mda

# Input files
gro = "md-5fs.gro"
trr = "md-5fs.trr"

# Output files
coord_out = "coord_InAs.xyz"
vel_out   = "vel_InAs.xyz"

# Load trajectory
u = mda.Universe(gro, trr)

# Select atoms; adjust names if needed
sel = u.select_atoms("name In As")

# If your atom names are uppercase:
# sel = u.select_atoms("name IN AS")

with open(coord_out, "w") as fc, open(vel_out, "w") as fv:
    for ts in u.trajectory:
        n = len(sel)

        # Coordinates in Å from MDAnalysis by default
        pos = sel.positions

        # Velocities are usually read from GROMACS in Å/ps after MDAnalysis unit conversion
        vel = sel.velocities

        fc.write(f"{n}\n")
        fc.write(f"Frame {ts.frame} time_ps {ts.time:.6f}\n")

        fv.write(f"{n}\n")
        fv.write(f"Frame {ts.frame} time_ps {ts.time:.6f}\n")

        for atom, r, v in zip(sel.atoms, pos, vel):
            name = atom.name
            fc.write(f"{name:2s} {r[0]:16.8f} {r[1]:16.8f} {r[2]:16.8f}\n")
            fv.write(f"{name:2s} {v[0]:16.8f} {v[1]:16.8f} {v[2]:16.8f}\n")
