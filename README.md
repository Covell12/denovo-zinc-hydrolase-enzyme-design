# De novo zinc hydrolase design

A direct test of whether the current AI protein-design pipeline can
build a new protein around the catalytic zinc site of human carbonic
anhydrase II without starting from carbonic anhydrase / any natural protein.

Note that nothing was synthesized, only predicted.

## Target

The metal-binding catalytic core of carbonic anhydrase II, taken from PDB entry
[3KS3](https://www.rcsb.org/structure/3KS3) at 0.9 Å: three histidine
nitrogens and one threonine oxygen arranged around a zinc ion.

## Pipeline

| Stage | Tool | In → Out |
|---|---|---|
| Backbone generation | RFdiffusion2 | motif → 32 backbones |
| Geometry filter | — | 32 → 20 |
| Sequence design | LigandMPNN | 20 → 320 sequences |
| Structure prediction | Boltz-2 | 60 → 5 passing |

A design passed if Boltz-2 gave it pLDDT ≥ 80, its predicted structure fell
within 2.0 Å of the proposed backbone, exactly three histidines coordinated
the zinc, and no histidine nitrogen sat more than 0.30 Å from its target
position.

## Results

Five of the sixty sequences passed, a success rate of 8.3%.

**Top candidate — `bias_cond0_14_s1`**, 127 residues, pLDDT 90.3, 1.92 Å from
its design.

| Atom | Native | This design |
|---|---|---|
| His94 N → His6 | 2.03 Å | 1.94 Å |
| His96 N → His8 | 2.04 Å | 2.00 Å |
| His119 N → His16 | 2.02 Å | 2.02 Å |
| Thr199 O → Thr98 | 3.83 Å | 5.13 Å |

The zinc and the three histidines reproduce the native geometry to 0.153 Å.
Including the threonine, the match is 0.817 Å, however, the threonine oxygen landed
1.3 Å too far out to.

## Contents

- `inputs/` — the four-atom motif extracted from it (the catalytic core; this was the only input to RFdiffusion2)
- `results/` — `motif_geometry.csv` (32 backbones scored), `folded_metrics.csv` (60 predictions scored)
- `designs/` — the five passing designs
- `figures/` — figures of some designs
- `models/` — models of the top designs

## Hardware

NVIDIA GeForce RTX 4080 SUPER; it ran about 11 hours total (majority RFdiffusion2 runtime, which lasted for 9 h 36 m).

Made by Covell Bernstein as a supplement to an 8th grade history research paper. Thanks for your interest!
