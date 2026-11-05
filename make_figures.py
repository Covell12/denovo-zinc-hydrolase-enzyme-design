#!/usr/bin/env python3
"""Figures 1 and 2 for the de novo zinc hydrolase design."""
import sys, math, os, itertools

argv = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else sys.argv[1:]
PRED, REF = argv[0], argv[1]
OUT = argv[2] if len(argv) > 2 else "figures"
os.makedirs(OUT, exist_ok=True)

from pymol import cmd
import numpy as np


def atoms_of(fp):
    out = []
    for l in open(fp):
        if l.startswith(("ATOM", "HETATM")):
            try:
                out.append({"res": l[17:20].strip(), "ch": l[21], "num": int(l[22:26]),
                            "name": l[12:16].strip(),
                            "xyz": (float(l[30:38]), float(l[38:46]), float(l[46:54]))})
            except ValueError:
                pass
    return out


def site(fp):
    A = atoms_of(fp)
    zn = [a for a in A if a["res"] == "ZN" or a["name"] == "ZN"]
    if not zn:
        return [], None, None
    z = zn[0]
    his = [a for a in A if a["res"] == "HIS" and a["name"] in ("NE2", "ND1")
           and 0 < math.dist(a["xyz"], z["xyz"]) < 2.8]
    his.sort(key=lambda a: math.dist(a["xyz"], z["xyz"]))
    thr = [a for a in A if a["res"] == "THR" and a["name"] == "OG1"]
    t = min(thr, key=lambda a: math.dist(a["xyz"], z["xyz"])) if thr else None
    return his[:3], t, z


def kabsch(P, Q):
    P, Q = np.array(P, float), np.array(Q, float)
    P = P - P.mean(0); Q = Q - Q.mean(0)
    U, S, Vt = np.linalg.svd(P.T @ Q)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return float(np.sqrt((((P @ R.T) - Q) ** 2).sum() / len(P)))


p_his, p_thr, p_zn = site(PRED)
r_his, r_thr, r_zn = site(REF)
print(f"predicted: {len(p_his)} coordinating His, Thr={'yes' if p_thr else 'no'}")
print(f"reference: {len(r_his)} coordinating His, Thr={'yes' if r_thr else 'no'}")

best, best_rms = None, 1e9
if len(p_his) == 3 and len(r_his) == 3:
    for perm in itertools.permutations(range(3)):
        P = [p_his[perm[i]]["xyz"] for i in range(3)] + [p_zn["xyz"]]
        Q = [r_his[i]["xyz"] for i in range(3)] + [r_zn["xyz"]]
        if p_thr and r_thr:
            P.append(p_thr["xyz"]); Q.append(r_thr["xyz"])
        rms = kabsch(P, Q)
        if rms < best_rms:
            best_rms, best = rms, perm
    print(f"best motif pairing {best}: RMSD {best_rms:.3f} A")
    with open(os.path.join(OUT, "motif_rmsd.txt"), "w") as fh:
        fh.write(f"{best_rms:.3f}\n")

cmd.set("ray_opaque_background", 0)
cmd.set("antialias", 2)
cmd.set("stick_radius", 0.17)
cmd.set("ray_shadows", 0)
cmd.bg_color("white")

PSEL = " or ".join(f"(chain {a['ch']} and resi {a['num']})" for a in p_his)
if p_thr:
    PSEL += f" or (chain {p_thr['ch']} and resi {p_thr['num']})"
RSEL = " or ".join(f"(chain {a['ch']} and resi {a['num']})" for a in r_his)
if r_thr:
    RSEL += f" or (chain {r_thr['ch']} and resi {r_thr['num']})"

cmd.load(PRED, "design")
cmd.hide("everything")
cmd.show("cartoon", "design and polymer")
cmd.spectrum("count", "deepblue_white_firebrick", "design and polymer and name CA")
cmd.set("cartoon_transparency", 0.12)
cmd.show("spheres", "design and (resn ZN or name ZN)")
cmd.color("purple", "design and (resn ZN or name ZN)")
cmd.set("sphere_scale", 0.42, "design and (resn ZN or name ZN)")
cmd.show("sticks", f"design and ({PSEL}) and not hydro")
cmd.color("yellow", f"design and ({PSEL}) and elem C")
cmd.util.cnc(f"design and ({PSEL})")
cmd.orient("design and polymer")
cmd.zoom("design and polymer", 2.5)
cmd.png(os.path.join(OUT, "fig1_design.png"), width=2400, height=1900, dpi=300, ray=1)
print("wrote fig1_design.png")

cmd.set("cartoon_transparency", 0.55)
cmd.orient(f"design and ({PSEL})")
cmd.zoom(f"design and ({PSEL})", 3.5)
cmd.png(os.path.join(OUT, "fig1b_active_site.png"), width=2400, height=1900, dpi=300, ray=1)
print("wrote fig1b_active_site.png")

cmd.delete("all")
cmd.load(PRED, "design")
cmd.load(REF, "native")
if best is not None:
    pa = [f"design and chain {p_his[best[i]]['ch']} and resi {p_his[best[i]]['num']} and name {p_his[best[i]]['name']}" for i in range(3)]
    qa = [f"native and chain {r_his[i]['ch']} and resi {r_his[i]['num']} and name {r_his[i]['name']}" for i in range(3)]
    pa.append("design and (resn ZN or name ZN)"); qa.append("native and (resn ZN or name ZN)")
    if p_thr and r_thr:
        pa.append(f"design and chain {p_thr['ch']} and resi {p_thr['num']} and name OG1")
        qa.append(f"native and chain {r_thr['ch']} and resi {r_thr['num']} and name OG1")
    cmd.pair_fit(*[x for pair in zip(pa, qa) for x in pair])

cmd.hide("everything")
cmd.show("sticks", f"design and ({PSEL}) and not hydro")
cmd.show("sticks", f"native and ({RSEL}) and not hydro")
cmd.color("marine", f"design and ({PSEL})")
cmd.color("orange", f"native and ({RSEL})")
cmd.show("spheres", "(resn ZN or name ZN)")
cmd.color("marine", "design and (resn ZN or name ZN)")
cmd.color("orange", "native and (resn ZN or name ZN)")
cmd.set("sphere_scale", 0.28)
cmd.orient(f"design and ({PSEL})")
cmd.zoom(f"design and ({PSEL})", 2.0)
cmd.png(os.path.join(OUT, "fig2_motif_overlay.png"), width=2400, height=1900, dpi=300, ray=1)
print("wrote fig2_motif_overlay.png")
print()
print(f"motif RMSD = {best_rms:.3f} A" if best is not None else "overlay skipped")
