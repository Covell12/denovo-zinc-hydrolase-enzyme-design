import sys, math, os, itertools
argv = sys.argv[sys.argv.index("--")+1:] if "--" in sys.argv else sys.argv[1:]
PRED, REF, OUT = argv[0], argv[1], argv[2]
os.makedirs(OUT, exist_ok=True)
from pymol import cmd
import numpy as np

def atoms_of(fp):
    out=[]
    for l in open(fp):
        if l.startswith(("ATOM","HETATM")):
            try: out.append({"res":l[17:20].strip(),"ch":l[21],"num":int(l[22:26]),
                             "name":l[12:16].strip(),
                             "xyz":(float(l[30:38]),float(l[38:46]),float(l[46:54]))})
            except ValueError: pass
    return out

def site(fp):
    A=atoms_of(fp)
    z=[a for a in A if a["res"]=="ZN" or a["name"]=="ZN"][0]
    his=sorted([a for a in A if a["res"]=="HIS" and a["name"] in ("NE2","ND1")
                and 0<math.dist(a["xyz"],z["xyz"])<2.8],
               key=lambda a: math.dist(a["xyz"],z["xyz"]))[:3]
    thr=min([a for a in A if a["res"]=="THR" and a["name"]=="OG1"],
            key=lambda a: math.dist(a["xyz"],z["xyz"]))
    return his,thr,z

def kabsch(P,Q):
    P,Q=np.array(P,float),np.array(Q,float); P=P-P.mean(0); Q=Q-Q.mean(0)
    U,S,Vt=np.linalg.svd(P.T@Q); d=np.sign(np.linalg.det(Vt.T@U.T))
    return float(np.sqrt((((P@(Vt.T@np.diag([1,1,d])@U.T).T)-Q)**2).sum()/len(P)))

ph,pt,pz=site(PRED); rh,rt,rz=site(REF)
best,brms=None,1e9
for pm in itertools.permutations(range(3)):
    r=kabsch([ph[pm[i]]["xyz"] for i in range(3)]+[pz["xyz"]],
             [rh[i]["xyz"] for i in range(3)]+[rz["xyz"]])
    if r<brms: brms,best=r,pm
full=min(kabsch([ph[pm[i]]["xyz"] for i in range(3)]+[pz["xyz"],pt["xyz"]],
                [rh[i]["xyz"] for i in range(3)]+[rz["xyz"],rt["xyz"]])
         for pm in itertools.permutations(range(3)))
open(os.path.join(OUT,"stats.txt"),"w").write(
    f"core {brms:.3f}\nfull {full:.3f}\n"
    f"thr_design {math.dist(pt['xyz'],pz['xyz']):.2f}\nthr_native {math.dist(rt['xyz'],rz['xyz']):.2f}\n")
print(f"core {brms:.3f}  full {full:.3f}")

cmd.set("ray_opaque_background",1); cmd.set("antialias",2); cmd.set("ray_shadows",0)
cmd.set("stick_radius",0.17); cmd.set("dash_radius",0.05); cmd.set("dash_gap",0.3)
cmd.set("orthoscopic",1); cmd.set("ray_trace_mode",0)
cmd.bg_color("white")
cmd.set_color("dsgn",[0.1647,0.4706,0.8392])
cmd.set_color("natv",[0.9216,0.4078,0.2039])
cmd.load(PRED,"design"); cmd.load(REF,"native")

pa,qa=[],[]
for i in range(3):
    a=ph[best[i]]; b=rh[i]
    pa.append(f"design and chain {a['ch']} and resi {a['num']} and name {a['name']}")
    qa.append(f"native and chain {b['ch']} and resi {b['num']} and name {b['name']}")
pa.append("design and name ZN"); qa.append("native and name ZN")
cmd.pair_fit(*[x for p in zip(pa,qa) for x in p])

PH=" or ".join(f"(chain {a['ch']} and resi {a['num']})" for a in ph)
RH=" or ".join(f"(chain {a['ch']} and resi {a['num']})" for a in rh)
PT=f"(chain {pt['ch']} and resi {pt['num']})"; RT=f"(chain {rt['ch']} and resi {rt['num']})"

def style(include_thr):
    cmd.hide("everything")
    psel = f"{PH} or {PT}" if include_thr else PH
    rsel = f"{RH} or {RT}" if include_thr else RH
    cmd.show("sticks", f"(design and ({psel})) and not hydro")
    cmd.show("sticks", f"(native and ({rsel})) and not hydro")
    cmd.color("dsgn", f"design and ({psel})")
    cmd.color("natv",  f"native and ({rsel})")
    cmd.set("stick_transparency",0.5,"native")
    cmd.show("spheres","name ZN")
    cmd.set("sphere_scale",0.46,"name ZN")
    cmd.color("dsgn","design and name ZN"); cmd.color("natv","native and name ZN")
    cmd.set("sphere_transparency",0.45,"native and name ZN")
    for i,a in enumerate(ph):
        cmd.distance(f"cd{i}", f"design and chain {a['ch']} and resi {a['num']} and name {a['name']}",
                     "design and name ZN")
        cmd.color("grey40", f"cd{i}"); cmd.hide("labels", f"cd{i}")

style(False)
cmd.orient(f"design and ({PH}) or design and name ZN")
cmd.zoom(f"design and ({PH} or {PT}) or design and name ZN", 1.0)
cmd.turn("x",-10)
VIEW = cmd.get_view()
cmd.png(os.path.join(OUT,"panelA.png"), width=1700, height=1700, dpi=300, ray=1)

style(True)
cmd.distance("td", f"design and {PT} and name OG1", "design and name ZN")
cmd.distance("tn", f"native and {RT} and name OG1", "native and name ZN")
cmd.color("dsgn","td"); cmd.color("natv","tn")
cmd.hide("labels","td"); cmd.hide("labels","tn")
cmd.set("dash_radius",0.07,"td"); cmd.set("dash_radius",0.07,"tn")
cmd.set_view(VIEW)
cmd.png(os.path.join(OUT,"panelB.png"), width=1700, height=1700, dpi=300, ray=1)
print("panels rendered")
