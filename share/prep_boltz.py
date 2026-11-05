#!/usr/bin/env python3
"""Write Boltz-2 FASTA inputs for the top-confidence designs from each set."""
import glob, os, re, sys, json

N_PER_SET = int(os.environ.get("N_PER_SET", "30"))
OUT = "fasta"
SETS = sys.argv[1:] or ["sequences", "sequences_bias"]


def parse(folder):
    recs = []
    for fa in sorted(glob.glob(os.path.join(folder, "seqs", "*.fa"))):
        design = os.path.basename(fa).split("-atomized")[0]
        hdr, idx = None, 0
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
                return float(m.group(1)) if m else 0.0
            idx += 1
            recs.append({"design": design, "idx": idx, "seq": line,
                         "conf": g("overall_confidence"), "lig": g("ligand_confidence")})
            hdr = None
    return recs


os.makedirs(OUT, exist_ok=True)
manifest = []
for folder in SETS:
    label = "bias" if "bias" in folder else "nobias"
    recs = parse(folder)
    if not recs:
        print(f"{folder}: no sequences")
        continue
    recs.sort(key=lambda r: -r["conf"])
    chosen = recs[:N_PER_SET]
    for r in chosen:
        tag = f"{label}_{r['design'].replace('run_znhydrolase_','')}_s{r['idx']}"
        with open(os.path.join(OUT, tag + ".fasta"), "w") as fh:
            fh.write(f">A|protein|empty\n{r['seq']}\n>Z|ccd\nZN\n")
        manifest.append({"tag": tag, "set": label, "design": r["design"],
                         "idx": r["idx"], "conf": r["conf"], "lig": r["lig"],
                         "length": len(r["seq"]),
                         "ala_pct": round(100 * r["seq"].count("A") / len(r["seq"]), 1),
                         "seq": r["seq"]})
    print(f"{label:<8} {len(recs):>4} sequences available, wrote top {len(chosen)} "
          f"(conf {chosen[-1]['conf']:.4f}-{chosen[0]['conf']:.4f})")

json.dump(manifest, open("fasta_manifest.json", "w"), indent=2)
print()
print(f"wrote {len(manifest)} fasta files to {OUT}/ and fasta_manifest.json")
