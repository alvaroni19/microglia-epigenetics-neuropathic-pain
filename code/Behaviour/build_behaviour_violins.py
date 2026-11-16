"""Behaviour change-score violins drawn with the shared group palette
(light tint = control, saturated = CCI-IoN).
Inputs are the saved outputs of code/Behaviour/behaviour_reanalysis.R (no statistics recomputed here)."""
import sys, glob, numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib import font_manager as fm
for f in glob.glob(sys.prefix + "/fonts/arial*.ttf") + glob.glob("C:/Windows/Fonts/arial*.ttf"): fm.fontManager.addfont(f)
from matplotlib.lines import Line2D
from pathlib import Path
ROOT = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(ROOT / "code" / "Figures_v2"))
import original_style as ost  # untrimmed violins shared with the other panels
B = str(ROOT / "generated_outputs" / "Behaviour_reanalysis") + "/"
OUT = B
mpl.rcParams.update({"font.family": "Arial", "font.size": 7, "axes.linewidth": 0.8,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8, "pdf.fonttype": 42, "svg.fonttype": "none"})
GROUPS = ["Male", "Proestrus", "Estrus", "Diestrus"]
SAT = {"Male": "#FFA367", "Proestrus": "#78BCE5", "Estrus": "#80CA76", "Diestrus": "#FFDB3E"}
LIGHT = {"Male": "#FFE6D7", "Proestrus": "#DAECF8", "Estrus": "#DEF1D9", "Diestrus": "#FFFFCB"}
d = pd.read_csv(B + "plot_data.csv")
emm = pd.read_csv(B + "emm_HC3.csv")
st = pd.read_excel(B + "behaviour_reanalysis_stats.xlsx", sheet_name=None)
inj = st["injury_effect_by_group"].set_index("Group")
cci = st["CCI_phase_comparisons"].set_index("contrast")
pint = st["ANOVA_females_HC3"].set_index("term").loc["Condition:Group", "Pr(>F)"]
def star(p): return "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < .01 else "*" if p < .05 else "n.s."

fig = plt.figure(figsize=(125/25.4, 72/25.4))
axb = fig.add_axes([0.14, 0.24, 0.84, 0.62])
rng = np.random.default_rng(1); pos = {}; ext = {}
for i, g in enumerate(GROUPS):
    for j, (cond, col) in enumerate([("Control", LIGHT[g]), ("CCI", SAT[g])]):
        x = i * 2.6 + j; pos[(g, cond)] = x
        y = d.query("Group==@g and Condition==@cond")["Delta"].values
        lo_, hi_ = ost.violin(axb, y, x, col, edge="black", lw=0.7, width=0.8, point_size=0, median_w=0)
        ext[(g, cond)] = (lo_, hi_)
        axb.scatter(x + rng.uniform(-.12, .12, len(y)), y, s=4, c="black", lw=0, zorder=3)
        e = emm.query("Group==@g and Condition==@cond").iloc[0]
        axb.plot([x-.22, x+.22], [e.emmean]*2, c="black", lw=1.6, zorder=4)
        axb.plot([x, x], [e["lower.CL"], e["upper.CL"]], c="black", lw=0.9, zorder=4)
def bracket(ax, x1, x2, y, txt, h=.08):
    ax.plot([x1, x1, x2, x2], [y-h, y, y, y-h], c="black", lw=0.7)
    ax.text((x1+x2)/2, y+.03, txt, ha="center", va="bottom", fontsize=7 if txt == "n.s." else 8)
for g in GROUPS:
    top = max(ext[(g, "Control")][1], ext[(g, "CCI")][1])
    bracket(axb, pos[(g, "Control")], pos[(g, "CCI")], top + .25, star(inj.loc[g, "p_holm"]))
ytop = max(v[1] for v in ext.values()) + .75
bracket(axb, pos[("Proestrus", "CCI")], pos[("Diestrus", "CCI")], ytop, star(cci.loc["Proestrus - Diestrus", "p_holm_females"]))
axb.axhline(0, ls=":", c="grey", lw=0.7, zorder=0)
axb.set_xticks([pos[k] for k in pos]); axb.set_xticklabels(["Control", "CCI-IoN"] * 4, rotation=45, ha="right", rotation_mode="anchor")
axb.set_ylabel("Δ mechanical response (0.16 g)\n(test − baseline)", fontweight="bold", fontsize=8)
axb.set_ylim(min(v[0] for v in ext.values()) - .15, ytop + 1.0); axb.set_xlim(-.7, pos[("Diestrus", "CCI")] + .7)
axb.text(.01, .995, f"Condition × phase (females): P = {pint:.3f}", transform=axb.transAxes, ha="left", va="top", fontsize=6.5)
for s in ["top", "right"]: axb.spines[s].set_visible(False)
hd = [Line2D([], [], marker="o", ls="", mfc=SAT[g], mec="black", mew=.7, ms=6, label=g) for g in GROUPS]
fig.legend(handles=hd, loc="upper center", bbox_to_anchor=(0.58, 1.0), ncol=4, frameon=False, fontsize=8,
           handletextpad=.1, columnspacing=.6, title="Biological phase", title_fontsize=8)
fig.legends[0]._legend_box.align = "left"
# title placed inline, as in the in vivo RT-qPCR panels
lg = fig.legends[0]; lg.set_title(None)
fig.text(0.58 - 0.30, 0.945, "Biological phase", ha="right", va="center", fontsize=8)
for ext in ["png", "pdf", "svg"]:
    fig.savefig(f"{OUT}/behaviour_change_score_violins.{ext}", dpi=600 if ext == "png" else None)
print("saved")
