#!/usr/bin/env python3
"""Build LigandMPNN batch inputs for designs that passed the geometry filter."""
import csv, glob, os, json, math, sys

CSV = sys.argv[1] if len(sys.argv) > 1 else "results/motif_geometry.csv"
GLOB = sys.argv[2] if len(sys.argv) > 2 else os.path.expanduser(
    "~/RFdiffusion2/pipeline_outputs/*/run_znhydrolase*.pdb")

passing = sorted({r["design"] for r in csv.DictReader(open(CSV)) if r["pass"] == "1"})
print(f"{len(passing)} designs passed the geometry filter")

paths, fixed = {}, {}
seen = set()
for fp in sorted(glob.glob(GLOB)):
    name = os.path.basename(fp).replace(".pdb", "")
    if name not in passing or name in seen:
        continue
    seen.add(name)
    atoms = []
    for l in open(fp):
        if l.startswith(("ATOM", "HETATM")):
            try:
                atoms.append((l[17:20].strip(), l[21], int(l[22:26]), l[12:16].strip(),
                              float(l[30:38]), float(l[38:46]), float(l[46:54])))
            except ValueError:
                pass
    zn = [a for a in atoms if a[0] == "ZN"]
    if not zn:
        continue
    z = zn[0]
    zp = (z[4], z[5], z[6])
    cat = []
    for a in atoms:
        if a[0] == "HIS" and a[3] in ("NE2", "ND1") and 0 < math.dist((a[4], a[5], a[6]), zp) < 2.6:
            cat.append((a[1], a[2]))
    thrs = [(a[1], a[2], math.dist((a[4], a[5], a[6]), zp))
            for a in atoms if a[0] == "THR" and a[3] == "OG1"]
    if thrs:
        best = min(thrs, key=lambda t: t[2])
        cat.append((best[0], best[1]))
    cat = sorted(set(cat))
    paths[fp] = ""
    fixed[fp] = " ".join(f"{ch}{num}" for ch, num in cat)

json.dump(paths, open("passing.json", "w"), indent=2)
json.dump(fixed, open("fixed.json", "w"), indent=2)
print(f"wrote passing.json and fixed.json for {len(paths)} files")
print()
for fp, res in list(fixed.items())[:5]:
    print(f"  {os.path.basename(fp)[:46]:<46} lock: {res}")
if len(fixed) > 5:
    print(f"  ... and {len(fixed)-5} more")
nres = {len(v.split()) for v in fixed.values()}
print()
print(f"residues locked per design: {sorted(nres)} (expect 4)")
