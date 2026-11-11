#!/usr/bin/env python3
"""Score Boltz-2 predictions against their designs."""
import glob, json, os, math, csv, sys, itertools

MANIFEST = "fasta_manifest.json"
FOLDED = sys.argv[1] if len(sys.argv) > 1 else "folded"
DESIGN_GLOB = os.path.expanduser("~/RFdiffusion2/pipeline_outputs/*/run_znhydrolase*.pdb")

PLDDT_MIN, RMSD_MAX, COORD_TOL = 80.0, 2.0, 0.30


def load_atoms(fp):
    out = []
    for l in open(fp):
        if l.startswith(("ATOM", "HETATM")):
            try:
                out.append({"res": l[17:20].strip(), "chain": l[21],
                            "num": int(l[22:26]), "name": l[12:16].strip(),
                            "xyz": (float(l[30:38]), float(l[38:46]), float(l[46:54]))})
            except ValueError:
                pass
    return out


def kabsch_rmsd(P, Q):
    import numpy as np
    n = len(P)
    cp = [sum(p[i] for p in P) / n for i in range(3)]
    cq = [sum(q[i] for q in Q) / n for i in range(3)]
    P = [[p[i] - cp[i] for i in range(3)] for p in P]
    Q = [[q[i] - cq[i] for i in range(3)] for q in Q]
    H = [[sum(P[k][i] * Q[k][j] for k in range(n)) for j in range(3)] for i in range(3)]
    U, S, Vt = np.linalg.svd(np.array(H))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    diff = np.array(P) @ R.T - np.array(Q)
    return float(np.sqrt((diff ** 2).sum() / n))


def zn_geometry(atoms):
    zn = [a for a in atoms if a["res"] in ("ZN", " ZN") or a["name"] == "ZN"]
    if not zn:
        return None
    zp = zn[0]["xyz"]
    nn = [a for a in atoms
          if a["res"] == "HIS" and a["name"] in ("NE2", "ND1")
          and 0 < math.dist(a["xyz"], zp) < 2.8]
    nn.sort(key=lambda a: math.dist(a["xyz"], zp))
    ds = [math.dist(a["xyz"], zp) for a in nn]
    angs = []
    if len(nn) >= 3:
        for a, b in itertools.combinations(nn[:3], 2):
            v1 = [a["xyz"][i] - zp[i] for i in range(3)]
            v2 = [b["xyz"][i] - zp[i] for i in range(3)]
            dot = sum(v1[i] * v2[i] for i in range(3))
            n1 = math.sqrt(sum(x * x for x in v1)); n2 = math.sqrt(sum(x * x for x in v2))
            angs.append(math.degrees(math.acos(max(-1, min(1, dot / (n1 * n2))))))
    return {"n": len(nn), "d": ds[:3], "ang": angs}


manifest = {m["tag"]: m for m in json.load(open(MANIFEST))}
designs = {}
for fp in glob.glob(DESIGN_GLOB):
    designs[os.path.basename(fp).split("-atomized")[0]] = fp

rows = []
for pdb in sorted(glob.glob(os.path.join(FOLDED, "**", "*_model_0.pdb"), recursive=True)):
    tag = os.path.basename(pdb).replace("_model_0.pdb", "")
    m = manifest.get(tag)
    if m is None:
        continue
    conf_fp = glob.glob(os.path.join(os.path.dirname(pdb), f"confidence_{tag}_model_0.json")) \
        or glob.glob(os.path.join(os.path.dirname(pdb), "confidence*.json"))
    plddt = ""
    if conf_fp:
        c = json.load(open(conf_fp[0]))
        for k in ("complex_plddt", "plddt", "confidence_score"):
            if k in c:
                v = c[k]
                plddt = round(v * 100 if v <= 1 else v, 1)
                break
    pred = load_atoms(pdb)
    dfp = designs.get(m["design"])
    rmsd = ""
    if dfp:
        des = load_atoms(dfp)
        pc = [a["xyz"] for a in pred if a["name"] == "CA"]
        dc = [a["xyz"] for a in des if a["name"] == "CA"]
        k = min(len(pc), len(dc))
        if k >= 50:
            rmsd = round(kabsch_rmsd(pc[:k], dc[:k]), 2)
    g = zn_geometry(pred)
    row = {"tag": tag, "set": m["set"], "design": m["design"],
           "mpnn_conf": m["conf"], "ala_pct": m["ala_pct"],
           "plddt": plddt, "rmsd": rmsd,
           "n_coord": g["n"] if g else 0,
           "zn_d1": round(g["d"][0], 2) if g and len(g["d"]) > 0 else "",
           "zn_d2": round(g["d"][1], 2) if g and len(g["d"]) > 1 else "",
           "zn_d3": round(g["d"][2], 2) if g and len(g["d"]) > 2 else "",
           "max_d_dev": round(max(abs(d - 2.03) for d in g["d"]), 2) if g and g["d"] else "",
           "seq": m["seq"]}
    ok = (plddt != "" and plddt >= PLDDT_MIN
          and rmsd != "" and rmsd <= RMSD_MAX
          and row["n_coord"] == 3
          and row["max_d_dev"] != "" and row["max_d_dev"] <= COORD_TOL)
    row["pass"] = 1 if ok else 0
    rows.append(row)

if not rows:
    print("no predictions found yet")
    sys.exit()

os.makedirs("results", exist_ok=True)
cols = ["tag", "set", "design", "mpnn_conf", "ala_pct", "plddt", "rmsd",
        "n_coord", "zn_d1", "zn_d2", "zn_d3", "max_d_dev", "pass", "seq"]
with open("results/folded_metrics.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols)
    w.writeheader()
    w.writerows(rows)

rows.sort(key=lambda r: (-(r["pass"]), -(r["plddt"] or 0), r["rmsd"] if r["rmsd"] != "" else 99))
hdr = f"{'tag':<30} {'set':<7} {'pLDDT':>6} {'RMSD':>5} {'nZn':>3} {'dDev':>5} {'Ala%':>5} {'pass':>4}"
print(hdr); print("-" * len(hdr))
for r in rows[:25]:
    print(f"{r['tag'][:30]:<30} {r['set']:<7} {str(r['plddt']):>6} {str(r['rmsd']):>5} "
          f"{r['n_coord']:>3} {str(r['max_d_dev']):>5} {r['ala_pct']:>5} {r['pass']:>4}")
if len(rows) > 25:
    print(f"... and {len(rows)-25} more (full table in results/folded_metrics.csv)")

print()
import statistics as st
for s in ("nobias", "bias"):
    sub = [r for r in rows if r["set"] == s]
    if not sub:
        continue
    pl = [r["plddt"] for r in sub if r["plddt"] != ""]
    rm = [r["rmsd"] for r in sub if r["rmsd"] != ""]
    print(f"{s:<8} n={len(sub):>3}  pLDDT median {st.median(pl) if pl else 0:>5.1f}  "
          f"RMSD median {st.median(rm) if rm else 0:>5.2f}  "
          f"3-coord {sum(1 for r in sub if r['n_coord']==3):>3}/{len(sub)}  "
          f"PASS {sum(r['pass'] for r in sub):>3}/{len(sub)}")
print()
print(f"TOTAL PASS: {sum(r['pass'] for r in rows)}/{len(rows)}")
print(f"criteria: pLDDT>={PLDDT_MIN}, RMSD<={RMSD_MAX} A, 3 His coordinating, max Zn-N deviation <={COORD_TOL} A")
