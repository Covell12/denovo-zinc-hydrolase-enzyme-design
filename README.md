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

In order to pass, a design had to satisfy four conditions: it required a
confidence score of at least 80 from Boltz-2, Boltz-2's predicted structure
had to be within 2.0 Å of the backbone that RFdiffusion2 had proposed,
the three histidines had to surround the zinc ion, and the histidine nitrogens had
to deviate no farther than 0.30 Å from their positions in the target motif.

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

### Design labels

The paper names designs `ZH⟨backbone⟩-⟨sequence⟩` in table 2; these proteins correspond to the following files: ZH14-1,
ZH14-4, ZH14-5 and ZH14-6 are `bias_cond0_14_s1`, `_s4`, `_s5` and `_s6`;
ZH28-5 is `nobias_cond0_28_s5`.

Made by Covell Bernstein as a supplement to an 8th grade history research paper. Thanks for your interest!
