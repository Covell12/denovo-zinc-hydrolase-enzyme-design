import sys, math
fp = sys.argv[1]
atoms = []
for l in open(fp):
    if l.startswith(("ATOM", "HETATM")):
        atoms.append((l[17:20].strip(), int(l[22:26]), l[12:16].strip(),
                      float(l[30:38]), float(l[38:46]), float(l[46:54])))
res = {(a[0], a[1]) for a in atoms}
print(f"atoms: {len(atoms)}   residues: {len(res)}")
zn = [a for a in atoms if a[0] == "ZN"]
if not zn:
    print("NO ZINC -- the metal was dropped")
    sys.exit()
z = zn[0]
print("zinc present")
near = []
for r, n, nm, x, y, zc in atoms:
    if r == "ZN":
        continue
    d = math.dist((x, y, zc), (z[3], z[4], z[5]))
    if d < 2.8:
        near.append((r, n, nm, d))
print(f"atoms within 2.8 A of zinc: {len(near)}")
for r, n, nm, d in sorted(near, key=lambda t: t[3]):
    print(f"  {r}{n:>4} {nm:<4} {d:.2f} A")
