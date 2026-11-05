#!/usr/bin/env python3
"""Geometry filter for RFdiffusion2 zinc-hydrolase designs."""
import sys, math, csv, os, itertools

REF_D = [2.03, 2.04, 2.02]
REF_A = [105.1, 113.6, 99.5]
REF_THR = 3.83
D_TOL, A_TOL = 0.25, 20.0


def load(fp):
    out = []
    for l in open(fp):
        if l.startswith(("ATOM", "HETATM")):
            try:
                out.append((l[17:20].strip(), int(l[22:26]), l[12:16].strip(),
                            float(l[30:38]), float(l[38:46]), float(l[46:54])))
            except ValueError:
                pass
    return out


def angle(p, c, q):
    v1 = [p[i] - c[i] for i in range(3)]
    v2 = [q[i] - c[i] for i in range(3)]
    d = sum(v1[i] * v2[i] for i in range(3))
    n1 = math.sqrt(sum(x * x for x in v1))
    n2 = math.sqrt(sum(x * x for x in v2))
    if n1 == 0 or n2 == 0:
        return float("nan")
    return math.degrees(math.acos(max(-1.0, min(1.0, d / (n1 * n2)))))


def analyse(fp):
    atoms = load(fp)
    r = {"design": os.path.basename(fp).replace(".pdb", ""),
         "residues": len({(a[0], a[1]) for a in atoms}), "zinc": 0,
         "n_coord": 0, "d1": "", "d2": "", "d3": "", "max_d_dev": "",
         "max_ang_dev": "", "thr_og1": "", "blocking_atom": "", "blocked": "",
         "pass": ""}
    zn = [a for a in atoms if a[0] == "ZN"]
    if not zn:
        return r
    r["zinc"] = 1
    z = zn[0]
    zp = (z[3], z[4], z[5])
    coord = [a for a in atoms
             if a[0] == "HIS" and a[2] in ("NE2", "ND1")
             and 0 < math.dist((a[3], a[4], a[5]), zp) < 2.6]
    coord.sort(key=lambda a: math.dist((a[3], a[4], a[5]), zp))
    r["n_coord"] = len(coord)
    ds = [math.dist((a[3], a[4], a[5]), zp) for a in coord]
    for i, d in enumerate(ds[:3]):
        r[f"d{i+1}"] = round(d, 2)
    if ds:
        r["max_d_dev"] = round(max(abs(d - 2.03) for d in ds), 2)
    if len(coord) >= 3:
        angs = [angle((a[3], a[4], a[5]), zp, (b[3], b[4], b[5]))
                for a, b in itertools.combinations(coord[:3], 2)]
        r["max_ang_dev"] = round(max(abs(a - b) for a, b in zip(sorted(angs), sorted(REF_A))), 1)
    thr = [math.dist((a[3], a[4], a[5]), zp)
           for a in atoms if a[0] == "THR" and a[2] == "OG1"]
    if thr:
        r["thr_og1"] = round(min(thr), 2)
    coord_ids = {(a[0], a[1], a[2]) for a in coord}
    intruders = []
    for a in atoms:
        if a[0] == "ZN" or (a[0], a[1], a[2]) in coord_ids:
            continue
        d = math.dist((a[3], a[4], a[5]), zp)
        if d < 2.8:
            intruders.append((f"{a[0]}{a[1]}:{a[2]}", round(d, 2)))
    intruders.sort(key=lambda t: t[1])
    r["blocked"] = 1 if intruders else 0
    r["blocking_atom"] = f"{intruders[0][0]}@{intruders[0][1]}" if intruders else ""
    ok = (r["n_coord"] == 3
          and r["max_d_dev"] != "" and r["max_d_dev"] <= D_TOL
          and r["max_ang_dev"] != "" and r["max_ang_dev"] <= A_TOL
          and not intruders)
    r["pass"] = 1 if ok else 0
    return r


rows = [analyse(fp) for fp in sorted(sys.argv[1:])]
cols = ["design", "residues", "zinc", "n_coord", "d1", "d2", "d3",
        "max_d_dev", "max_ang_dev", "thr_og1", "blocked", "blocking_atom", "pass"]
with open("results/motif_geometry.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols)
    w.writeheader()
    w.writerows(rows)

hdr = f"{'design':<42} {'res':>4} {'d1':>5} {'d2':>5} {'d3':>5} {'dDev':>5} {'aDev':>5} {'Thr':>5} {'blk':>4} {'pass':>4}"
print(hdr)
print("-" * len(hdr))
for r in rows:
    print(f"{r['design'][:42]:<42} {r['residues']:>4} {str(r['d1']):>5} {str(r['d2']):>5} "
          f"{str(r['d3']):>5} {str(r['max_d_dev']):>5} {str(r['max_ang_dev']):>5} "
          f"{str(r['thr_og1']):>5} {str(r['blocked']):>4} {str(r['pass']):>4}")

n = len(rows)
if n:
    print()
    print(f"designs analysed      : {n}")
    print(f"zinc retained         : {sum(r['zinc'] for r in rows)}/{n}")
    print(f"three His coordinating: {sum(1 for r in rows if r['n_coord']==3)}/{n}")
    print(f"4th position blocked  : {sum(1 for r in rows if r['blocked']==1)}/{n}")
    print(f"PASS all geometry     : {sum(1 for r in rows if r['pass']==1)}/{n}")
    print()
    print("reference (3KS3): d = 2.03 2.04 2.02 A, angles 105.1 113.6 99.5 deg, Thr OG1 3.83 A")
    print("wrote results/motif_geometry.csv")
