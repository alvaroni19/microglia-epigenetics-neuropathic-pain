"""Drawing helpers that reproduce the design of the manuscript figures.

Conventions
-----------
MeDIP panels (directionality, GO, heatmap, Csf1r locus, UpSet, window map) write CSF1 contrasts
treatment first (F-CSF1 vs F-C, M-CSF1 vs M-C), as in every other figure,
and sex contrasts male first (M-C vs F-C, M-CSF1 vs F-CSF1).  The release
tables store F_CSF_vs_F_C, M_CSF_vs_M_C (kept as they are) and F_C_vs_M_C,
F_CSF_vs_M_CSF (log2FC multiplied by -1 and Hyper/Hypo exchanged); the ERE
scan table stores the CSF1 contrasts control first (F_C_vs_F_CSF,
M_C_vs_M_CSF), which are inverted as well.  "Hyper" always means higher
methylation in the first group of the label.

RT-qPCR panels keep treatment-versus-reference contrasts (M-CSF1 vs M-C,
M-CCI-IoN vs M-C), as in the manuscript text.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle

MM = 1 / 25.4

# ----------------------------------------------------------------------------- palettes
# log2FC scale sampled from the colour bar of the manuscript bubble matrices.
LOGFC_STOPS = ["#002E66", "#0069B1", "#0097C7", "#78C8E0", "#CBE7F1", "#F8F7F7",
               "#FFDCC8", "#FFA680", "#F5604A", "#D21128", "#7C0019"]
LOGFC_CMAP = LinearSegmentedColormap.from_list("manuscript_logfc", LOGFC_STOPS, N=256)

CSF1_HYPO, CSF1_HYPER = "#0073B6", "#F45E00"         # DMR directionality, CSF1 effect
SEX_FEMALE_HIGH, SEX_MALE_HIGH = "#00A071", "#E57AAA"  # DMR directionality, sex effect
HYPER, HYPO = "#C0392B", "#2471A3"                   # DMR direction

SAT = {"Male": "#FFA367", "Proestrus": "#78BCE5", "Estrus": "#80CA76", "Diestrus": "#FFDB3E"}
LIGHT = {"Male": "#FFE6D7", "Proestrus": "#DAECF8", "Estrus": "#DEF1D9", "Diestrus": "#FFFFCB"}
SEX_EDGE = {"male": "#FFA367", "female": "#6666FF"}
TREAT_FILL = {"control": "#FFFFFF", "CSF1": "#DADADA", "CSF1E2": "#808080"}

# ----------------------------------------------------------------------------- MeDIP contrasts
# release key -> (display label, sign applied to log2FC)
MEDIP_CONTRAST = {
    "F_CSF_vs_F_C": ("F-CSF1 vs F-C", 1),
    "M_CSF_vs_M_C": ("M-CSF1 vs M-C", 1),
    "F_C_vs_M_C": ("M-C vs F-C", -1),
    "F_CSF_vs_M_CSF": ("M-CSF1 vs F-CSF1", -1),
    # keys written control first in the ERE scan table
    "F_C_vs_F_CSF": ("F-CSF1 vs F-C", -1),
    "M_C_vs_M_CSF": ("M-CSF1 vs M-C", -1),
}
SWAP_DIRECTION = {"HYPER": "HYPO", "HYPO": "HYPER", "Hyper": "Hypo", "Hypo": "Hyper"}


def medip_label(key: str) -> str:
    return MEDIP_CONTRAST[key][0]


def medip_direction(key: str, direction: str) -> str:
    """Direction relative to the first group of the displayed label."""
    if MEDIP_CONTRAST[key][1] == -1:
        return SWAP_DIRECTION.get(direction, direction)
    return direction


# ----------------------------------------------------------------------------- generic
def save(fig, stem, dpi=600):
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    from matplotlib.legend import Legend
    extra = fig.findobj(Legend)   # legends added with add_artist are otherwise ignored
    for ext in ("svg", "pdf", "png"):
        fig.savefig(stem.with_suffix(f".{ext}"), dpi=dpi, bbox_inches="tight",
                    bbox_extra_artists=extra, pad_inches=0.02, facecolor="white")
    plt.close(fig)


def finish(ax, labelsize=6):
    ax.spines[["top", "right"]].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(0.6)
    ax.tick_params(labelsize=labelsize, width=0.6, length=2, pad=1)


def bracket(ax, x0, x1, y, h, text, fontsize=7):
    ax.plot([x0, x0, x1, x1], [y, y + h, y + h, y], color="black", lw=0.6, zorder=7)
    ax.text((x0 + x1) / 2, y + 0.2 * h, text, ha="center", va="bottom",
            fontsize=fontsize, zorder=7)


# ----------------------------------------------------------------------------- violins
def _bw_nrd0(x):
    """R's default bandwidth (bw.nrd0), used by ggplot2::geom_violin."""
    x = np.asarray(x, float)
    sd = np.std(x, ddof=1)
    iqr = np.subtract(*np.percentile(x, [75, 25]))
    lo = min(sd, iqr / 1.34) if iqr > 0 else sd
    if lo <= 0:
        lo = abs(x[0]) if x[0] != 0 else 1.0
    return 0.9 * lo * len(x) ** (-0.2)


def violin(ax, vals, x, face, edge="black", lw=0.6, width=0.8, point_size=3.0,
           median_w=0.36, seed=None):
    """Violin with untrimmed Gaussian KDE tails (ggplot2 trim = FALSE).

    The density is evaluated from min - 3*bw to max + 3*bw, so both ends taper
    to a point.  Returns the lower and upper y extent of the violin.
    """
    vals = np.asarray(vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        return np.nan, np.nan
    if len(vals) >= 2 and np.ptp(vals) > 0:
        bw = _bw_nrd0(vals)
        ys = np.linspace(vals.min() - 3 * bw, vals.max() + 3 * bw, 512)
        dens = np.exp(-0.5 * ((ys[:, None] - vals[None, :]) / bw) ** 2).sum(axis=1)
        half = dens / dens.max() * width / 2
        ax.fill_betweenx(ys, x - half, x + half, facecolor=face, edgecolor=edge,
                         lw=lw, zorder=2, joinstyle="round")
        lo, hi = ys.min(), ys.max()
    else:
        lo, hi = vals.min(), vals.max()
    med = np.median(vals)
    ax.plot([x - median_w / 2, x + median_w / 2], [med, med], color="black", lw=1.1,
            zorder=5, solid_capstyle="butt")
    rng = np.random.default_rng(seed if seed is not None else int(abs(x) * 100) + len(vals))
    ax.scatter(x + rng.uniform(-0.06, 0.06, len(vals)), vals, s=point_size, color="black",
               lw=0, zorder=6)
    return lo, hi


# ----------------------------------------------------------------------------- bubble matrix
def bubble_matrix(df, comps, genes, col_labels, cell_mm=8.0, lim=None,
                  gene_fs=7, tick_fs=6.5, cbar_label="log$_2$ fold change"):
    """Compact bubble matrix: adjacent grey cells, gene names and colour bar
    placed next to the grid.

    df needs columns comparison, gene_label, log2fc, stars.
    """
    df = df[df["comparison"].isin(comps)].copy()
    genes = [g for g in genes if g in set(df["gene_label"])]
    nr, nc = len(genes), len(comps)
    if lim is None:
        lim = float(np.nanmax(np.abs(df["log2fc"].to_numpy())))
        lim = (lim if np.isfinite(lim) and lim > 0 else 1.0) * 1.05
    grid_w, grid_h = nc * cell_mm, nr * cell_mm
    left, bottom, gap, cb_w = 14.0, 4.0, 2.0, 2.6
    fig_w, fig_h = left + grid_w + gap + cb_w + 12.0, bottom + grid_h + 2.0
    fig = plt.figure(figsize=(fig_w * MM, fig_h * MM))
    ax = fig.add_axes([left / fig_w, bottom / fig_h, grid_w / fig_w, grid_h / fig_h])
    ymap = {g: nr - 1 - i for i, g in enumerate(genes)}
    xmap = {c: i for i, c in enumerate(comps)}
    for x in range(nc):
        for y in range(nr):
            ax.add_patch(Rectangle((x - 0.5, y - 0.5), 1, 1, facecolor="#F7F7F7",
                                   edgecolor="#E2E2E2", lw=0.6, zorder=0))
    norm = Normalize(-lim, lim)
    diam_pt = 0.82 * cell_mm * MM * 72
    for _, r in df.iterrows():
        if r["gene_label"] not in ymap:
            continue
        x, y = xmap[r["comparison"]], ymap[r["gene_label"]]
        col = LOGFC_CMAP(norm(r["log2fc"])) if np.isfinite(r["log2fc"]) else (1, 1, 1, 1)
        ax.scatter(x, y, s=diam_pt ** 2, facecolor=col, edgecolor="#333333", lw=0.5, zorder=3)
        if isinstance(r["stars"], str) and r["stars"]:
            t = ax.text(x, y, r["stars"], ha="center", va="center", fontsize=7.5,
                        fontweight="bold", color="white", zorder=4)
            t.set_path_effects([pe.withStroke(linewidth=1.3, foreground="black")])
    ax.set_xlim(-0.5, nc - 0.5)
    ax.set_ylim(-0.5, nr - 0.5)
    ax.set_xticks(range(nc))
    ax.set_xticklabels([col_labels.get(c, c) for c in comps], rotation=30, ha="right",
                       rotation_mode="anchor", fontsize=tick_fs)
    ax.set_yticks([ymap[g] for g in genes])
    ax.set_yticklabels(genes, fontsize=gene_fs, fontstyle="italic")
    ax.tick_params(length=0, pad=1.5)
    for sp in ax.spines.values():
        sp.set_visible(False)
    cb_h = min(grid_h, max(24.0, 0.8 * grid_h))
    cax = fig.add_axes([(left + grid_w + gap) / fig_w,
                        (bottom + (grid_h - cb_h) / 2) / fig_h, cb_w / fig_w, cb_h / fig_h])
    sm = plt.cm.ScalarMappable(cmap=LOGFC_CMAP, norm=norm)
    cb = fig.colorbar(sm, cax=cax)
    cb.solids.set_rasterized(False)
    cb.set_label(cbar_label, fontsize=7, fontweight="bold", labelpad=2)
    cb.ax.tick_params(labelsize=6, width=0.4, length=1.5, pad=1)
    cb.outline.set_linewidth(0.4)
    return fig
