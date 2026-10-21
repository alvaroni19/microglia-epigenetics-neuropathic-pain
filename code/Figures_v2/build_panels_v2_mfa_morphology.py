"""Figure panels 4c, 5c, 5e and 6c-f.

Values are read from the release tables; no test is recomputed except the
display of already-defined quantities:
  marker-pain Spearman       generated_outputs/PCR_MFA/invivo_marker_pain_spearman.csv
  in vivo MFA                generated_outputs/PCR_MFA/invivo_mfa_coords.csv, invivo_mfa_eig.csv
  in vitro estrogenic score  generated_outputs/PCR_260514/invitro_estrogenic_score_per_replicate.csv (directed
           delta score = treatment value - median of its reference group, as in
           the release builder), invitro_estrogen_genes_pca_scores.csv / _variance.csv
  ramification index         generated_outputs/MorphOMICs/cam_sholl_per_cell.csv (ipsilateral,
           lamina I-II) and the mixed-model table
           cam_ramification_index_layer12_animal_and_mixedmodel.csv (stars)
  morphology UMAP            data/MorphOMICs/saved_instances_cam/...cam_reduced.csv (ipsilateral,
           lamina I-II; zoom on the main cloud, all bags in the inset)
  persistence difference     persistence-image centroids and representative SWC traces, as in
           code/MorphOMICs/plot_persistence_diff_recolored.py (difference map only)
  ramification-pain          partial Spearman of ranks adjusted for sex/estrous phase, as in
           code/MorphOMICs/plot_ramification_pain_correlation_layer12_english.py
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Ellipse, Patch
from scipy.ndimage import gaussian_filter
from scipy.stats import chi2

ROOT = Path(__file__).absolute().parents[2]
CODE, GEN = ROOT / "code", ROOT / "generated_outputs"
OUT = GEN / "Figures_v2"

QPCR_VITRO = ROOT / "data" / "PCR_260514" / "invitro_qpcr_relative_expression.xlsx"
QPCR_VIVO = ROOT / "data" / "PCR_260514" / "invivo_qpcr_delta_ct.xlsx"
BEHAVIOUR = ROOT / "data" / "Allodynia_0_16" / "allodynia_change_scores.xlsx"
VON_FREY = ROOT / "data" / "Allodynia_0_16" / "allodynia_von_frey_scores.xlsx"
LINKAGE = ROOT / "data" / "MorphOMICs" / "morphology_behaviour_linkage.xlsx"
RESTRICTED = ("behaviour and RT-qPCR inputs are not distributed "
              "(available from the corresponding author on request)")


def _gate(label, *paths):
    """True when every input exists; otherwise report the panel as skipped."""
    missing = [p for p in paths if not Path(p).exists()]
    if missing:
        print(f"[SKIP] {label}: {missing[0].name} not found; {RESTRICTED}")
    return not missing

for sub in ("Figures_v2", "_helpers", "MorphOMICs"):
    sys.path.insert(0, str(CODE / sub))
import original_style as st  # noqa: E402

COND = {"Control": "#999999", "CCI": "#E57AAA"}     # morphology palette
SEX_FILL = {"Male": "#FFA367", "Female": "#6666FF"}
PHASE_FILL = {"Male": st.SAT["Male"], "Proestrus": st.SAT["Proestrus"],
              "Estrus": st.SAT["Estrus"], "Diestrus": st.SAT["Diestrus"]}
PHASE_EDGE = {"Male": "#E0661E", "Proestrus": "#3E95C8", "Estrus": "#2E9E57", "Diestrus": "#E2BC00"}
CYCLE_EN = {"Male": "Male", "Proestrus": "Proestrus", "Estrus": "Estrus", "Diestrus": "Diestrus"}
GREEK = {"Gadd45a": "Gadd45α", "Gadd45b": "Gadd45β", "Gadd45g": "Gadd45γ", "Tnf": "Tnf-α"}


def _cov_ellipse(ax, x, y, face, edge, z, level=0.95):
    if len(x) < 3:
        return
    vals, vecs = np.linalg.eigh(np.cov(x, y))
    o = vals.argsort()[::-1]
    vals, vecs = vals[o], vecs[:, o]
    ang = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    s = np.sqrt(chi2.ppf(level, 2))
    w, h = 2 * s * np.sqrt(np.maximum(vals, 1e-9))
    ax.add_patch(Ellipse((np.mean(x), np.mean(y)), w, h, angle=ang, facecolor=face, alpha=0.22,
                         edgecolor="none", zorder=z))
    ax.add_patch(Ellipse((np.mean(x), np.mean(y)), w, h, angle=ang, facecolor="none", edgecolor=edge,
                         lw=0.9, zorder=z + 0.1))


# ----------------------------------------------------------------------------- marker-pain Spearman
def panel_invivo_marker_pain_spearman():
    df = pd.read_csv(GEN / "PCR_MFA" / "invivo_marker_pain_spearman.csv").sort_values("rho")
    df["lab"] = df["gene_label"].map(lambda g: GREEK.get(g, g))
    y = np.arange(len(df))
    fig, ax = plt.subplots(figsize=(62 * st.MM, (10 + 4.2 * len(df)) * st.MM))
    ax.hlines(y, 0, df["rho"], color="#BBBBBB", lw=0.7, zorder=1)
    for yi, (_, r) in zip(y, df.iterrows()):
        col = "#D7301F" if r["rho"] > 0 else "#2C7FB8"
        sig = r["q_BH"] < 0.05
        ax.scatter(r["rho"], yi, s=70 if sig else 14, color=col, edgecolor="#555555" if sig else "none",
                   lw=0.8, zorder=3)
        if sig:
            ax.text(r["rho"], yi, "*", ha="center", va="center", fontsize=9, color="white",
                    fontweight="bold", zorder=4)
    ax.axvline(0, color="#777777", lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(df["lab"], fontsize=7, fontstyle="italic")
    ax.set_xlim(-1, 1)
    ax.set_xticks(np.arange(-1, 1.01, 0.25))
    ax.set_xticklabels([f"{v:.2f}".replace("-", "−") for v in np.arange(-1, 1.01, 0.25)], fontsize=6)
    ax.set_xlabel("Spearman correlation ρ", fontsize=7, labelpad=1)
    ax.set_ylim(-0.6, len(df) - 0.4)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="y", length=0, pad=1)
    st.save(fig, OUT / "panels" / "invivo_marker_pain_spearman")


# ----------------------------------------------------------------------------- in vivo MFA
def panel_invivo_mfa_phase_condition():
    df = pd.read_csv(GEN / "PCR_MFA" / "invivo_mfa_coords.csv")
    eig = pd.read_csv(GEN / "PCR_MFA" / "invivo_mfa_eig.csv")
    fig, ax = plt.subplots(figsize=(80 * st.MM, 58 * st.MM))
    ax.axvline(0, color="#BBBBBB", lw=0.5, zorder=1)
    ax.axhline(0, color="#BBBBBB", lw=0.5, zorder=1)
    for i, ph in enumerate(["Male", "Diestrus", "Estrus", "Proestrus"]):
        sub = df[df["phase"] == ph]
        _cov_ellipse(ax, sub["Dim1"].to_numpy(), sub["Dim2"].to_numpy(), PHASE_FILL[ph], PHASE_EDGE[ph], 2 + 2 * i,
                     level=0.68)  # 68% normal data ellipse (~1 SD) per phase
        for cond, mk in (("Control", "o"), ("CCI", "^")):
            sub_c = sub[sub["condition"] == cond]
            ax.scatter(sub_c["Dim1"], sub_c["Dim2"], marker=mk, s=20, facecolor=PHASE_FILL[ph], edgecolor="#333333",
                       lw=0.4, zorder=3 + 2 * i)
    h1 = [Line2D([0], [0], marker="o", ls="", mfc=PHASE_FILL[p], mec="#333333", ms=4.5, label=p)
          for p in ["Male", "Proestrus", "Estrus", "Diestrus"]]
    h2 = [Line2D([0], [0], marker=m, ls="", mfc="#999999", mec="#333333", ms=4.5, label=l)
          for m, l in (("o", "Control"), ("^", "CCI-IoN"))]
    l1 = ax.legend(handles=h1, title="Biological phase", loc="upper left", bbox_to_anchor=(1.0, 1.02),
                   frameon=False, fontsize=6.5, title_fontsize=7, handletextpad=0.1, labelspacing=0.15,
                   borderaxespad=0.2)
    l1._legend_box.align = "left"
    ax.add_artist(l1)
    l2 = ax.legend(handles=h2, title="Condition", loc="upper left", bbox_to_anchor=(1.0, 0.50), frameon=False,
                   fontsize=6.5, title_fontsize=7, handletextpad=0.1, labelspacing=0.15, borderaxespad=0.2)
    l2._legend_box.align = "left"
    ax.set_xlabel(f"Dimension 1 ({float(eig['var1'][0]):.1f}%)", fontsize=7, fontweight="bold", labelpad=1)
    ax.set_ylabel(f"Dimension 2 ({float(eig['var2'][0]):.1f}%)", fontsize=7, fontweight="bold", labelpad=1)
    st.finish(ax)
    st.save(fig, OUT / "panels" / "invivo_mfa_phase_condition")


# ----------------------------------------------------------------------------- in vitro estrogenic score and PCA
def panel_invitro_estrogenic_score_and_pca():
    scr = GEN / "PCR_260514"
    score = pd.read_csv(scr / "invitro_estrogenic_score_per_replicate.csv")
    sc = pd.read_csv(scr / "invitro_estrogen_genes_pca_scores.csv")
    var = pd.read_csv(scr / "invitro_estrogen_genes_pca_variance.csv").set_index("PC")["variance_explained_pct"]
    contrasts = [("Male", "CSF1", "Control", "M-CSF1\nvs M-C"), ("Male", "CSF1+E2", "CSF1", "M-CSF1+E2\nvs M-CSF1"),
                 ("Female", "CSF1", "Control", "F-CSF1\nvs F-C"), ("Female", "CSF1+E2", "CSF1", "F-CSF1+E2\nvs F-CSF1")]
    fig, ax = plt.subplots(figsize=(82 * st.MM, 46 * st.MM))
    rng = np.random.default_rng(20260608)
    rows = []
    for x, (sex, num, den, lab) in enumerate(contrasts):
        a = score.loc[(score.sex == sex) & (score.treatment == num), "ESI_full"].dropna().to_numpy(float)
        b = score.loc[(score.sex == sex) & (score.treatment == den), "ESI_full"].dropna().to_numpy(float)
        d = a - np.median(b)
        rows += [dict(contrast=lab.replace("\n", " "), replicate=i + 1, directed_delta=v) for i, v in enumerate(d)]
        ax.bar(x, np.median(d), width=0.62, color=SEX_FILL[sex], alpha=0.35, edgecolor=SEX_FILL[sex], lw=0.6, zorder=2)
        ax.scatter(x + rng.normal(0, 0.05, len(d)), d, s=16, facecolor=SEX_FILL[sex], edgecolor="#333333",
                   lw=0.4, zorder=3)
    pd.DataFrame(rows).to_csv(OUT / "panels" / "invitro_estrogenic_score_values.csv", index=False)
    ax.axhline(0, color="#555555", lw=0.6, ls=(0, (1, 1)))
    ax.axvline(1.5, color="#BBBBBB", lw=0.5)
    ax.set_xticks(range(4))
    ax.set_xticklabels([c[3] for c in contrasts], fontsize=6)
    ax.set_ylabel("Δ estrogenic score", fontsize=7, fontweight="bold", labelpad=1)
    ax.grid(axis="y", color="#E6E6E6", lw=0.4)
    ax.set_axisbelow(True)
    st.finish(ax)
    ax.legend(handles=[Patch(fc=SEX_FILL[s], alpha=0.5, ec=SEX_FILL[s], label=s) for s in ("Male", "Female")],
              title="Biological group", loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, frameon=False,
              fontsize=6.5, title_fontsize=7, handlelength=1.2, columnspacing=1.0)
    st.save(fig, OUT / "panels" / "invitro_estrogenic_score")

    fig, ax = plt.subplots(figsize=(78 * st.MM, 52 * st.MM))
    for sex in ("Male", "Female"):
        s = sc[sc.sex == sex]
        _cov_ellipse(ax, s["PC1"].to_numpy(), s["PC2"].to_numpy(), SEX_FILL[sex], SEX_FILL[sex], 1)
    shapes = {"Control": "o", "CSF1": "s", "CSF1+E2": "^"}
    for tr, mk in shapes.items():
        for sex in ("Male", "Female"):
            g = sc[(sc.sex == sex) & (sc.treatment == tr)]
            ax.scatter(g["PC1"], g["PC2"], marker=mk, s=20, facecolor=SEX_FILL[sex], edgecolor="#333333",
                       lw=0.4, zorder=3)
    ax.axhline(0, color="#777777", lw=0.4, ls=(0, (1, 1)))
    ax.axvline(0, color="#777777", lw=0.4, ls=(0, (1, 1)))
    ax.set_xlabel(f"Dimension 1 ({var['PC1']:.0f}%)", fontsize=7, fontweight="bold", labelpad=1)
    ax.set_ylabel(f"Dimension 2 ({var['PC2']:.0f}%)", fontsize=7, fontweight="bold", labelpad=1)
    st.finish(ax)
    lg = ax.legend(handles=[Line2D([0], [0], marker=m, ls="", mfc="white", mec="black", ms=4.5, label=t)
                            for t, m in shapes.items()], title="Condition", loc="lower right", frameon=False,
                   fontsize=6.5, title_fontsize=7, handletextpad=0.1, labelspacing=0.15)
    lg._legend_box.align = "left"
    st.save(fig, OUT / "panels" / "invitro_estrogen_genes_pca")


# ----------------------------------------------------------------------------- ramification index
def _stars(p):
    return "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else ""


def panel_morphology_ramification_index():
    pc = GEN / "MorphOMICs" / "cam_sholl_per_cell.csv"
    if not pc.exists():
        pc = ROOT / "data" / "MorphOMICs" / "reports_cam" / "cam_sholl_per_cell.csv"
    d = pd.read_csv(pc)
    d = d[(d["Layer"] == 12) & (d["Side"] == "IPSI")].dropna(subset=["ShollRamificationIndex"])
    lmm = pd.read_csv(GEN / "MorphOMICs" / "cam_ramification_index_layer12_animal_and_mixedmodel.csv").set_index("cycle")
    cycles = ["Male", "Proestrus", "Estrus", "Diestrus"]
    fig, axes = plt.subplots(1, 4, figsize=(104 * st.MM, 42 * st.MM), sharey=True,
                             gridspec_kw={"wspace": 0.05})
    rng = np.random.default_rng(7)
    ri = "ShollRamificationIndex"
    lo = d[ri].min()
    for k, (ax, cy) in enumerate(zip(axes, cycles)):
        sub = d[d["Phase"].replace({"Metestrus": "Diestrus"}) == cy]
        tops = []
        for x, g in ((1, "Control"), (2, "CCI")):
            v = sub.loc[sub["Condition"] == g, ri].to_numpy(float)
            _, hi = st.violin(ax, v, x, COND[g], width=0.8, lw=0.5, point_size=0, median_w=0)
            ax.scatter(x + rng.normal(0, 0.06, len(v)), v, s=4, color="#333333", alpha=0.35, lw=0, zorder=3)
            ax.scatter([x], [v.mean()], s=60, color="black", edgecolor="white", lw=1.2, zorder=5)
            tops.append(hi)
        s = _stars(float(lmm.loc[CYCLE_EN[cy], "p_value"]))
        if s:
            yb = max(tops) + 0.15
            st.bracket(ax, 1, 2, yb, 0.25, s, fontsize=8)
        ax.set_title(CYCLE_EN[cy].upper(), fontsize=7.5, fontweight="bold", pad=2)
        ax.set_xticks([1, 2])
        ax.set_xticklabels(["Control", "CCI-IoN"], fontsize=6.5)
        ax.set_xlim(0.4, 2.6)
        st.finish(ax)
        ax.tick_params(axis="x", length=0)
        if k:
            ax.spines["left"].set_visible(True)
            ax.tick_params(axis="y", left=False)
    axes[0].set_ylabel("Ramification index", fontsize=7, fontweight="bold", labelpad=1)
    axes[0].set_ylim(min(0.8, lo - 0.2), axes[0].get_ylim()[1] + 0.4)
    axes[0].set_yticks(range(1, int(axes[0].get_ylim()[1]) + 1))
    st.save(fig, OUT / "panels" / "morphology_ramification_index")


# ----------------------------------------------------------------------------- morphology UMAP
def panel_morphology_umap():
    """UMAP panels zoomed on the main cloud (robust distance <= 6 MAD units from the
    median, as in code/MorphOMICs/plot_umap_layer12_full_and_zoom.py); every
    ipsilateral lamina I-II bag stays visible in the 'Full view' inset."""
    from matplotlib.patches import Rectangle
    u = pd.read_csv(ROOT / "data" / "MorphOMICs" / "saved_instances_cam" /
                    "Morphomics.PID_CAM_CCI_260224.cam_reduced.csv", index_col=0)
    u["Phase"] = u["Phase"].astype(str).replace({"Metestrus": "Diestrus"})
    u = u[(u["Layer"].astype(str) == "12") & (u["Side"] == "IPSI")].copy()
    xy = u[["umap___1", "umap___2"]]
    med = xy.median()
    mad = (xy - med).abs().median().replace(0, 1)
    u["in_zoom"] = np.sqrt((((xy - med) / mad) ** 2).sum(axis=1)) <= 6.0

    def lim(v, f=0.04):
        lo, hi = float(v.min()), float(v.max())
        pad = max((hi - lo) * f, 0.25)
        return lo - pad, hi + pad
    zx, zy = lim(u.loc[u["in_zoom"], "umap___1"]), lim(u.loc[u["in_zoom"], "umap___2"])
    # headroom so that the 'Full view' inset (upper-right 38% of each panel) never covers bags
    IN = (0.61, 0.60, 0.37, 0.38)
    x_in = zx[0] + IN[0] * (zx[1] - zx[0])
    y_need = u.loc[u["in_zoom"] & (u["umap___1"] >= x_in - 0.3), "umap___2"].max() + 0.3
    zy = (zy[0], max(zy[1], zy[0] + (y_need - zy[0]) / (IN[1] - 0.02)))
    fx, fy = lim(u["umap___1"]), lim(u["umap___2"])
    rep = (u.groupby(["Phase", "Condition"])["in_zoom"].agg(total_bags="size", bags_in_zoom="sum").reset_index())
    rep["bags_outside_zoom"] = rep["total_bags"] - rep["bags_in_zoom"]
    rep.to_csv(OUT / "panels" / "morphology_umap_zoom_report.csv", index=False)
    fig, axes = plt.subplots(1, 4, figsize=(118 * st.MM, 40 * st.MM), sharex=True, sharey=True,
                             gridspec_kw={"wspace": 0.10})
    for k, (ax, cy) in enumerate(zip(axes, ["Male", "Proestrus", "Estrus", "Diestrus"])):
        sub = u[u["Phase"] == cy]
        ins = ax.inset_axes(list(IN))
        for g in ("Control", "CCI"):
            s_ = sub[sub["Condition"] == g]
            z = 3 if g == "CCI" else 2
            zs = s_[s_["in_zoom"]]
            ax.scatter(zs["umap___1"], zs["umap___2"], s=4, color=COND[g], lw=0, alpha=0.9, zorder=z)
            ins.scatter(s_["umap___1"], s_["umap___2"], s=1.5, color=COND[g], lw=0, alpha=0.9, zorder=z)
        ins.set_xlim(fx); ins.set_ylim(fy); ins.set_xticks([]); ins.set_yticks([])
        for sp in ins.spines.values():
            sp.set_visible(True); sp.set_linewidth(0.5)
        ins.add_patch(Rectangle((zx[0], zy[0]), zx[1] - zx[0], zy[1] - zy[0], fill=False, ec="#333333", lw=0.5))
        ins.text(0.05, 0.04, "Full view", transform=ins.transAxes, fontsize=5.5, ha="left", va="bottom")
        ax.set_title(CYCLE_EN[cy], fontsize=7.5, fontweight="bold", pad=2)
        ax.set_xlim(zx); ax.set_ylim(zy)
        ax.set_xlabel("UMAP 1", fontsize=7, fontweight="bold", labelpad=1)
        st.finish(ax, labelsize=5.5)
        if k:
            ax.tick_params(axis="y", left=False)
    axes[0].set_ylabel("UMAP 2", fontsize=7, fontweight="bold", labelpad=1)
    st.save(fig, OUT / "panels" / "morphology_umap")


# ----------------------------------------------------------------------------- persistence difference
def panel_morphology_persistence_difference():
    from morphomics_plot_utils import compute_persistence_centroids, pick_representative_swc, read_swc_dataframe
    data = ROOT / "data" / "MorphOMICs"
    with np.load(data / "saved_instances_cam" / "cam_persistence_vectors.npz", allow_pickle=False) as z:
        cols = ["Condition", "Phase", "Sex", "Side", "Layer"]
        dv = pd.DataFrame({c: z[c] for c in cols})
        dv["pi"] = list(z["persistence_images"])
    dv["Phase"] = dv["Phase"].replace({"Metestrus": "Diestrus"})
    cen = compute_persistence_centroids(dv, group_columns=["Condition", "Phase", "Side", "Layer"], vector_column="pi")
    man = pd.read_csv(data / "reports_cam" / "cam_manifest_selected.csv")
    man["Phase"] = man["Phase"].astype(str).replace({"Metestrus": "Diestrus"})
    m12 = man[(man["side"].astype(str) == "IPSI") & (man["layer"].astype(str) == "12")]
    diff_cmap = LinearSegmentedColormap.from_list("ctrl_cci", ["#3C3C3C", "#F4F4F4", COND["CCI"]])
    keyed = {}
    for k, v in cen.items():
        c = v["combo"]
        if str(c.get("Side")) == "IPSI" and str(c.get("Layer")) == "12":
            keyed[(c["Condition"], c["Phase"])] = v
    for cy in ["Male", "Proestrus", "Estrus", "Diestrus"]:
        a, b = keyed[("CCI", cy)], keyed[("Control", cy)]
        ai, bi = gaussian_filter(a["centroid"], 1.2), gaussian_filter(b["centroid"], 1.2)
        diff = gaussian_filter(ai / (ai.max() + 1e-12) - bi / (bi.max() + 1e-12), 1.0)
        vmax = float(np.abs(diff).max()) or 1.0
        fig = plt.figure(figsize=(86 * st.MM, 40 * st.MM))
        ax = fig.add_axes([0.10, 0.20, 0.36, 0.68])
        cax = fig.add_axes([0.475, 0.20, 0.018, 0.68])
        im = ax.imshow(diff, origin="lower", cmap=diff_cmap, vmin=-vmax, vmax=vmax)
        ax.set_xlabel("Birth", fontsize=7, fontweight="bold", labelpad=1)
        ax.set_ylabel("Death", fontsize=7, fontweight="bold", labelpad=1)
        ax.tick_params(labelsize=6, length=2, width=0.5, pad=1)
        n = diff.shape[0]
        ticks = np.linspace(0, n - 1, 5)
        ax.set_xticks(ticks); ax.set_yticks(ticks)
        ax.set_xticklabels([f"{t:.0f}" for t in np.linspace(0, 80, 5)])
        ax.set_yticklabels([f"{t:.0f}" for t in np.linspace(0, 80, 5)])
        for sp in ax.spines.values():
            sp.set_linewidth(0.5)
        cb = fig.colorbar(im, cax=cax)
        cb.ax.tick_params(labelsize=5.5, length=1.5, width=0.4, pad=1)
        cb.outline.set_linewidth(0.4)
        ph = CYCLE_EN[cy]
        ax.set_title(f"{ph} CCI-IoN / {ph} C", fontsize=7, fontweight="bold", pad=3, loc="left")
        tax = fig.add_axes([0.60, 0.20, 0.38, 0.62])
        for cond, col, lab in (("CCI", COND["CCI"], "CCI-IoN"), ("Control", "#8C8C8C", "Control")):
            sub = m12[(m12["Condition"].astype(str) == cond) & (m12["Phase"] == cy)]
            paths = [Path(p) if Path(str(p)).is_absolute() else ROOT / str(p) for p in sub["target_path"].dropna()]
            rep = pick_representative_swc([p for p in paths if p.exists()])
            sw = read_swc_dataframe(rep)
            ctr = sw.loc[sw["parent"] == -1, ["x", "y"]].iloc[0].to_numpy(float) if (sw["parent"] == -1).any() \
                else sw[["x", "y"]].mean().to_numpy(float)
            xy = {int(r.n): np.array([r.x, r.y]) - ctr for r in sw.itertuples(index=False)}
            for r in sw.itertuples(index=False):
                if int(r.parent) in xy:
                    p1, p2 = xy[int(r.n)], xy[int(r.parent)]
                    tax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=col, lw=0.7)
            tax.plot([], [], color=col, lw=2.0, label=lab)
        tax.set_aspect("equal"); tax.axis("off")
        tax.set_title("Representative\nSWC traces", fontsize=7, pad=2)
        tax.legend(loc="upper left", bbox_to_anchor=(0.55, 0.08), frameon=False, fontsize=6.5, handlelength=1.0,
                   handletextpad=0.3, labelspacing=0.1, borderaxespad=0)
        st.save(fig, OUT / "panels" / f"morphology_persistence_difference_{ph}")


# ----------------------------------------------------------------------------- ramification-pain partial Spearman
def panel_morphology_ramification_pain_partial_spearman():
    spec = importlib.util.spec_from_file_location(
        "corr", CODE / "MorphOMICs" / "plot_ramification_pain_correlation_layer12_english.py")
    corr = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(corr)
    src = GEN / "MorphOMICs" / "cam_ramification_pain_correlation_layer12_IPSI_source.csv"
    if not src.exists():
        src = ROOT / "data" / "MorphOMICs" / "cam_ramification_pain_correlation_layer12_IPSI_source.csv"
    d = pd.read_csv(src)
    d = d[(d["side"] == "IPSI") & (d["lamina"] == 12)][["animal_id", "condition", "sex_or_phase",
                                                         "ShollRamificationIndex", "test_0.16g"]].dropna()
    d["sex_or_phase"] = d["sex_or_phase"].replace({"Metestrus": "Diestrus"})
    d = d.drop_duplicates("animal_id").reset_index(drop=True)
    rho, p, dof, fx, py = corr.partial_spearman(d, "ShollRamificationIndex", "test_0.16g", ["sex_or_phase"])
    pd.DataFrame([dict(n_animals=len(d), partial_rho=rho, p_value=p, df=dof)]).to_csv(
        OUT / "panels" / "morphology_ramification_pain_partial_spearman_values.csv", index=False)
    d["fx"], d["py"] = fx, py
    marks = {"Male": "s", "Proestrus": "o", "Estrus": "^", "Diestrus": "D"}
    fig, ax = plt.subplots(figsize=(64 * st.MM, 48 * st.MM))
    ax.grid(True, color="#E3E3E3", lw=0.6, zorder=0)
    b, a = np.polyfit(fx, py, 1)
    xs = np.linspace(fx.min(), fx.max(), 50)
    ax.plot(xs, a + b * xs, color="#202020", lw=1.0, zorder=2)
    for _, r in d.iterrows():
        ax.scatter(r["fx"], r["py"], marker=marks[r["sex_or_phase"]], s=24,
                   facecolor=COND["CCI"] if r["condition"] == "CCI" else COND["Control"],
                   edgecolor="#333333", lw=0.5, zorder=3)
    ax.set_xlabel("Ramification residual rank", fontsize=7, fontweight="bold", labelpad=1)
    ax.set_ylabel("Pain residual rank", fontsize=7, fontweight="bold", labelpad=1)
    st.finish(ax)
    h1 = [Line2D([0], [0], marker="o", ls="", mfc=COND[c], mec="#333333", ms=4.5, label=l)
          for c, l in (("Control", "Control"), ("CCI", "CCI-IoN"))]
    h2 = [Line2D([0], [0], marker=m, ls="", mfc="white", mec="black", ms=4.5, label=k) for k, m in marks.items()]
    l1 = ax.legend(handles=h1, title="Condition", loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False,
                   fontsize=6.5, title_fontsize=7, handletextpad=0.1, labelspacing=0.15)
    l1._legend_box.align = "left"
    ax.add_artist(l1)
    l2 = ax.legend(handles=h2, title="Sex/estrous phase", loc="upper left", bbox_to_anchor=(1.02, 0.62),
                   frameon=False, fontsize=6.5, title_fontsize=7, handletextpad=0.1, labelspacing=0.15)
    l2._legend_box.align = "left"
    st.save(fig, OUT / "panels" / "morphology_ramification_pain_partial_spearman")
    return rho, p, len(d)


def main():
    plt.rcParams.update({"font.family": "Arial", "svg.fonttype": "none", "pdf.fonttype": 42})
    panel_morphology_ramification_index(); panel_morphology_umap(); panel_morphology_persistence_difference()
    if _gate("in vitro estrogenic score panel", QPCR_VITRO):
        panel_invitro_estrogenic_score_and_pca()
    if _gate("in vivo correlation and MFA panels", QPCR_VIVO, VON_FREY):
        panel_invivo_marker_pain_spearman(); panel_invivo_mfa_phase_condition()
    if _gate("ramification-pain panel", VON_FREY, LINKAGE):
        rho, p, n = panel_morphology_ramification_pain_partial_spearman()
        print(f"[OK] ramification-pain partial rho={rho:.3f} P={p:.4f} n={n}")
    print("[OK] correlation, MFA, estrogenic-score and morphology panels")


if __name__ == "__main__":
    main()
