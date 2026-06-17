#!/usr/bin/env python3
import sys

inp, out, skip = sys.argv[1], sys.argv[2], int(sys.argv[3])

with open(inp) as f, open(out, "w") as g:
    frame = 0
    while True:
        line = f.readline()
        if not line:
            break

        natom = int(line.strip())
        comment = f.readline()
        atoms = [f.readline() for _ in range(natom)]

        if frame >= skip:
            g.write(f"{natom}\n")
            g.write(comment)
            g.writelines(atoms)

        frame += 1
