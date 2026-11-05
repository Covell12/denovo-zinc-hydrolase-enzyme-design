#!/usr/bin/env python3
"""Center a motif PDB on its own centroid and append an ORI pseudo-atom."""
import sys

src, dst = sys.argv[1], sys.argv[2]
lines = [l.rstrip("\n") for l in open(src)]

coords = []
for l in lines:
    if l.startswith(("ATOM", "HETATM")):
        coords.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))

cx = sum(c[0] for c in coords) / len(coords)
cy = sum(c[1] for c in coords) / len(coords)
cz = sum(c[2] for c in coords) / len(coords)
print(f"centroid of {len(coords)} atoms: {cx:.3f} {cy:.3f} {cz:.3f}")

out, serial = [], 0
for l in lines:
    if l.startswith(("ATOM", "HETATM")):
        serial = int(l[6:11])
        x = float(l[30:38]) - cx
        y = float(l[38:46]) - cy
        z = float(l[46:54]) - cz
        out.append(f"{l[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{l[54:]}")
    elif l.startswith("END"):
        continue
    else:
        out.append(l)

ori = f"HETATM{serial+1:>5}  ORI ORI B{999:>4}    {0.0:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           X  "
out.append(ori)
out.append("END")

open(dst, "w").write("\n".join(out) + "\n")
print("wrote", dst)
