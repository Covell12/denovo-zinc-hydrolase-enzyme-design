#!/usr/bin/env python3
"""Compare LigandMPNN sequence sets on composition and confidence."""
import glob, os, sys, re, statistics as st
from collections import Counter


def parse(folder):
    recs = []
    for fa in sorted(glob.glob(os.path.join(folder, "seqs", "*.fa"))):
        hdr = None
        for line in open(fa):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                hdr = line
                continue
            if hdr is None or "overall_confidence" not in hdr:
                hdr = None
                continue
            def g(k):
                m = re.search(rf"{k}=([0-9.]+)", hdr)
                return float(m.group(1)) if m else None
            recs.append({"file": os.path.basename(fa), "seq": line,
                         "conf": g("overall_confidence"),
                         "lig": g("ligand_confidence"),
                         "rec": g("seq_rec")})
            hdr = None
    return recs


def summarise(name, recs):
    if not recs:
        print(f"{name}: no sequences found")
        return None
    ala = [100 * s["seq"].count("A") / len(s["seq"]) for s in recs]
    comp = Counter("".join(s["seq"] for s in recs))
    tot = sum(comp.values())
    conf = [s["conf"] for s in recs if s["conf"] is not None]
    lig = [s["lig"] for s in recs if s["lig"] is not None]
    rec = [s["rec"] for s in recs if s["rec"] is not None]
    print(f"--- {name} ---")
    print(f"  sequences          : {len(recs)}")
    print(f"  length             : {st.median([len(s['seq']) for s in recs]):.0f}")
    print(f"  alanine            : median {st.median(ala):.1f}%  range {min(ala):.1f}-{max(ala):.1f}%")
    if conf: print(f"  overall_confidence : median {st.median(conf):.4f}")
    if lig:  print(f"  ligand_confidence  : median {st.median(lig):.4f}")
    if rec:  print(f"  seq_rec            : median {st.median(rec):.4f}")
    top = "  ".join(f"{aa}:{100*c/tot:.0f}%" for aa, c in comp.most_common(6))
    print(f"  composition        : {top}")
    return {"ala": st.median(ala), "conf": st.median(conf) if conf else 0,
            "lig": st.median(lig) if lig else 0, "n": len(recs)}


sets = sys.argv[1:] or ["sequences", "sequences_bias"]
out = {}
for folder in sets:
    label = os.path.basename(folder.rstrip("/"))
    out[label] = summarise(label, parse(folder))
    print()

print("natural reference: Ala ~8%, Leu ~10%, Glu ~6%, Gly ~7%")
print()
valid = {k: v for k, v in out.items() if v}
if len(valid) == 2:
    pick = min(valid.items(), key=lambda kv: abs(kv[1]["ala"] - 8))
    print(f"closer to natural composition: {pick[0]} ({pick[1]['ala']:.1f}% Ala)")
    for k, v in valid.items():
        print(f"  {k:<18} Ala {v['ala']:>5.1f}%   conf {v['conf']:.4f}   lig {v['lig']:.4f}")
