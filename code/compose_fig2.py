import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.lines import Line2D
import os, sys

D = sys.argv[1] if len(sys.argv) > 1 else "figures/panels"
OUT = sys.argv[2] if len(sys.argv) > 2 else "figures/fig2_motif_overlay.png"
stats = dict(l.split() for l in open(os.path.join(D, "stats.txt")))

SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#8a8982"
BLUE, ORANGE = "#2a78d6", "#eb6834"

fig, axes = plt.subplots(1, 2, figsize=(10.6, 5.9), facecolor=SURFACE)
fig.subplots_adjust(left=0.02, right=0.98, top=0.80, bottom=0.10, wspace=0.03)

for ax, img, title, sub in (
    (axes[0], "panelA.png", "A   Catalytic core",
     f"Zn + three His nitrogens\nRMSD {stats['core']} \u00c5"),
    (axes[1], "panelB.png", "B   With the second-shell threonine",
     f"all five motif atoms\nRMSD {stats['full']} \u00c5"),
):
    ax.imshow(mpimg.imread(os.path.join(D, img)))
    ax.axis("off")
    ax.set_title(title, fontsize=11.6, color=INK, loc="left", pad=6)
    ax.text(0.015, 0.03, sub, transform=ax.transAxes, fontsize=9.6,
            color=INK2, va="bottom", ha="left", linespacing=1.5)

axes[1].text(0.985, 0.975,
             f"Thr OG1 \u2192 Zn    design {stats['thr_design']} \u00c5    native {stats['thr_native']} \u00c5",
             transform=axes[1].transAxes, fontsize=9.6, color=INK2,
             va="top", ha="right")

handles = [
    Line2D([], [], marker="o", color="none", markerfacecolor=BLUE,
           markeredgecolor=BLUE, markersize=10, label="This design (de novo, 128 aa)"),
    Line2D([], [], marker="o", color="none", markerfacecolor=ORANGE,
           markeredgecolor=ORANGE, markersize=10, label="Human carbonic anhydrase II (PDB 3KS3)"),
]
fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.02, 0.915),
           ncol=2, frameon=False, fontsize=9.4, labelcolor=INK2,
           handletextpad=0.5, columnspacing=2.2)

fig.suptitle("A protein with no evolutionary relationship to carbonic anhydrase "
             "reproduces its zinc site to 0.15 \u00c5",
             x=0.02, y=0.965, ha="left", fontsize=12.6, color=INK)
fig.text(0.02, 0.022,
         "Superimposed on the four core atoms. Grey dashes: Zn\u2013N coordination. "
         "The metal centre transfers within crystallographic error; the threonine that "
         "orients the\ncatalytic hydroxide sits 1.3 \u00c5 further out than in the native enzyme, "
         "beyond hydrogen-bonding range.",
         fontsize=7.8, color=MUTED, ha="left", va="bottom")

fig.savefig(OUT, dpi=300, facecolor=SURFACE)
print("wrote", OUT)
