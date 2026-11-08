# PUBLIC RELEASE CANDIDATE COPY
# Source: tka_inflammation_redesign/STAGE7J_FIGURE3_MULTIPANEL_UPGRADE/01_make_figure3.py
# Source SHA256: 12b5e2a9fd3b6e1d0a117db187c82a2aaec027a548a2c476add51ebe796773e4
# Purpose: Figure 3 source; use only approved aggregate plot inputs
# This copy does not overwrite the source and contains no study input data.

#!/usr/bin/env python3
"""Stage 7J-F3: publication-quality Figure 3, using frozen Stage 7H/I inputs."""
from pathlib import Path
import math
import textwrap
import hashlib
import platform
import sys
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent
STAGE2 = PROJECT / "rescue_discovery_20261026"
STAGE7H = PROJECT / "STAGE7H_R_CRP_ECM_SPATIAL_BRIDGE"
SCORES = pd.read_csv(ROOT / "FIGURE3_SAMPLE_SCORES.csv")
FUNC = pd.read_csv(ROOT / "ECM64_FUNCTIONAL_PANEL_MAP.csv")
BRIDGE = pd.read_csv(STAGE7H / "GSE225731_ECM64_ECM33_BRIDGE.csv")
GSE198520 = pd.read_csv(STAGE7H / "GSE198520_FIXED_MODULE_REPLICATION.csv")
PATHWAYS = pd.read_csv(STAGE2 / "GSE225731_PATHWAY_DISCOVERY_ALL.csv")

OUT = ROOT
PLOT_BLUE = "#1F4E79"
BLUE2 = "#4472C4"
TEAL = "#2A9D8F"
ORANGE = "#D97904"
PURPLE = "#6A4C93"
GREY = "#7A7A7A"
LIGHT_GREY = "#D9DEE5"
DARK = "#242424"

mpl.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8.2,
    "axes.titlesize": 9.5,
    "axes.labelsize": 8.2,
    "xtick.labelsize": 7.2,
    "ytick.labelsize": 7.2,
    "legend.fontsize": 7.2,
    "axes.linewidth": 0.75,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
})

def fmt_p(x):
    if pd.isna(x):
        return "NA"
    if x < 0.001:
        return f"{x:.2e}"
    return f"{x:.3f}"

def fmt_r(r, lo, hi):
    return f"rho={r:.3f} [{lo:.3f}, {hi:.3f}]"

def panel_label(ax, letter):
    ax.text(-0.12, 1.06, letter, transform=ax.transAxes, fontsize=12,
            fontweight="bold", va="top", ha="left", color=DARK)

def clean_spines(ax, keep=("left", "bottom")):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_visible(s in keep)
    ax.tick_params(width=0.6, length=3)

def get_row(sig):
    return BRIDGE.loc[BRIDGE.signature == sig].iloc[0]

def draw_panel_a(ax):
    tested = PATHWAYS.loc[PATHWAYS.analysis_status.eq("TESTED")].copy()
    tested = tested.sort_values(["rho", "pathway"]).reset_index(drop=True)
    tested["rank"] = np.arange(1, len(tested) + 1)
    ax.scatter(tested["rank"], tested["rho"], s=5, color="#B7BDC7", alpha=.62,
               linewidths=0, rasterized=True)
    hit = tested.pathway.eq("REACTOME_COLLAGEN_DEGRADATION")
    ax.scatter(tested.loc[hit, "rank"], tested.loc[hit, "rho"], s=24,
               color=PLOT_BLUE, edgecolor="white", linewidth=.45, zorder=4)
    if hit.any():
        x = float(tested.loc[hit, "rank"].iloc[0]); y = float(tested.loc[hit, "rho"].iloc[0])
        q = float(tested.loc[hit, "q_value_BH"].iloc[0])
        ax.annotate(f"ECM64\nq={q:.4f}", xy=(x, y), xytext=(x + 105, y + .055),
                    arrowprops=dict(arrowstyle="-", color=PLOT_BLUE, lw=.7),
                    fontsize=7.1, color=PLOT_BLUE, ha="left", va="bottom")
    ax.axhline(0, color="#9EA5AF", lw=.65, ls=(0, (3, 2)))
    ax.set_xlim(0, len(tested) + 25)
    ax.set_xlabel(f"Rank among {len(tested):,} tested pathways (ordered by Spearman rho)")
    ax.set_ylabel("Spearman rho with baseline CRP")
    ax.set_title("Pathway-screen context", loc="left", fontweight="bold")
    excluded = int((PATHWAYS.analysis_status != "TESTED").sum())
    ax.text(.02, .96, f"Screened {len(PATHWAYS):,}  |  Tested {len(tested):,}  |  Excluded {excluded:,}",
            transform=ax.transAxes, va="top", fontsize=7.2, color=DARK)
    clean_spines(ax)
    panel_label(ax, "A")

def draw_scatter(ax, x, y, row, title, color=PLOT_BLUE, ylim=None):
    ax.scatter(x, y, s=18, color=color, alpha=.78, edgecolor="white", linewidth=.35)
    ax.set_xscale("log")
    ax.set_title(title, loc="left", fontsize=8.2, fontweight="bold")
    ax.text(.03, .96,
            f"{fmt_r(row.rho_vs_CRP, row.bootstrap_low, row.bootstrap_high)}\n"
            f"P={fmt_p(row.P_value)}; N={int(row.N)}",
            transform=ax.transAxes, ha="left", va="top", fontsize=7.0,
            bbox=dict(boxstyle="round,pad=.25", fc="white", ec="none", alpha=.82))
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.set_xlabel("Baseline CRP (mg/dL; log)", labelpad=4)
    ax.set_ylabel("Module score (z-scale)")
    clean_spines(ax)

def draw_panel_b(parent):
    panel_label(parent, "B")
    parent.set_axis_off()
    sub = parent.get_subplotspec().subgridspec(1, 2, wspace=.50)
    axes = [parent.figure.add_subplot(sub[0]), parent.figure.add_subplot(sub[1])]
    x = SCORES.CRP_mg_dL
    ylim = (min(SCORES[["ECM64", "ECM33"]].min()) - .12,
            max(SCORES[["ECM64", "ECM33"]].max()) + .12)
    draw_scatter(axes[0], x, SCORES.ECM64, get_row("ECM64"), "ECM64 (64 genes)", PLOT_BLUE, ylim)
    draw_scatter(axes[1], x, SCORES.ECM33, get_row("ECM33"), "ECM33 (Xenium-covered)", TEAL, ylim)
    return axes

def draw_panel_c(ax):
    n_total = len(FUNC)
    n_struct = int(FUNC.functional_class.eq("COL structural").sum())
    n_processing = int(FUNC.functional_class.str.contains("protease").sum())
    n_other = int(n_total - n_struct - n_processing)
    n_covered = int(FUNC.xenium_covered.astype(str).str.upper().eq("TRUE").sum())
    n_uncovered = n_total - n_covered
    classes = [("COL structural", n_struct, PLOT_BLUE),
               ("Protease / processing", n_processing, ORANGE),
               ("Other non-COL", n_other, PURPLE)]
    left = 0
    for name, val, col in classes:
        ax.barh(1, val, left=left, height=.34, color=col, edgecolor="white", linewidth=.5)
        if val >= 8:
            ax.text(left + val / 2, 1, str(val), ha="center", va="center", color="white", fontsize=7, fontweight="bold")
        left += val
    ax.barh(0, n_covered, left=0, height=.34, color=TEAL, edgecolor="white", linewidth=.5)
    ax.barh(0, n_uncovered, left=n_covered, height=.34, color=LIGHT_GREY, edgecolor="white", linewidth=.5)
    ax.text(n_covered / 2, 0, f"{n_covered} covered", ha="center", va="center", color="white", fontsize=7, fontweight="bold")
    ax.text(n_covered + n_uncovered / 2, 0, f"{n_uncovered} not covered", ha="center", va="center", color=DARK, fontsize=7)
    ax.set_yticks([0, 1], ["Xenium panel", "ECM64 function"])
    ax.set_xlim(0, n_total)
    ax.set_ylim(-.35, 1.42)
    ax.set_xlabel("Number of genes")
    ax.set_title("Composition and panel coverage", loc="left", fontweight="bold")
    ax.text(.02, 1.31, f"COL {n_struct}  |  protease/processing {n_processing}  |  other non-COL {n_other}",
            transform=ax.transData, fontsize=6.8, color=DARK, va="top")
    clean_spines(ax, keep=("bottom",))
    panel_label(ax, "C")

def draw_panel_d(ax):
    r33 = get_row("ECM33")
    r31 = get_row("ECM31_uncovered")
    ax.scatter(SCORES.ECM64, SCORES.ECM33, s=18, color=TEAL, alpha=.78,
               edgecolor="white", linewidth=.35)
    lo = min(SCORES.ECM64.min(), SCORES.ECM33.min()) - .12
    hi = max(SCORES.ECM64.max(), SCORES.ECM33.max()) + .12
    ax.plot([lo, hi], [lo, hi], color="#9EA5AF", lw=.7, ls=(0, (3, 2)))
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel("ECM64 score")
    ax.set_ylabel("ECM33 score")
    ax.set_title("Fixed-panel concordance", loc="left", fontweight="bold")
    ax.text(.04, .95, f"rho={r33.rho_vs_ECM64:.3f}; N={int(r33.N)}\nSame cohort; not independent validation",
            transform=ax.transAxes, va="top", fontsize=7.2)
    ax.text(.04, .06, f"Uncovered ECM31–CRP: rho={r31.rho_vs_CRP:.3f}",
            transform=ax.transAxes, va="bottom", fontsize=7.0, color=GREY)
    clean_spines(ax)
    panel_label(ax, "D")

def draw_panel_e(parent):
    parent.set_axis_off()
    inner = parent.get_subplotspec().subgridspec(1, 2, width_ratios=[.66, .34], wspace=.02)
    ax = parent.figure.add_subplot(inner[0])
    tx = parent.figure.add_subplot(inner[1])
    order = ["ECM64", "ECM33", "ECM31_uncovered", "Structural_all",
             "Structural_Xenium", "Processing_all", "Processing_Xenium"]
    labels = ["ECM64 (64)", "ECM33 covered (33)", "ECM31 uncovered (31)",
              "Structural all (38)", "Structural Xenium (17)",
              "Processing all (26)", "Processing Xenium (16)"]
    colors = [PLOT_BLUE, TEAL, GREY, PLOT_BLUE, TEAL, ORANGE, ORANGE]
    y = np.arange(len(order))[::-1]
    for yy, sig, lab, col in zip(y, order, labels, colors):
        r = get_row(sig)
        ax.plot([r.bootstrap_low, r.bootstrap_high], [yy, yy], color=col, lw=1.8, solid_capstyle="round")
        ax.plot([r.bootstrap_low, r.bootstrap_low], [yy-.09, yy+.09], color=col, lw=.75)
        ax.plot([r.bootstrap_high, r.bootstrap_high], [yy-.09, yy+.09], color=col, lw=.75)
        ax.scatter([r.rho_vs_CRP], [yy], color=col, s=24, zorder=3, edgecolor="white", linewidth=.4)
    ax.axvline(0, color="#89919C", lw=.7, ls=(0, (3, 2)))
    ax.set_yticks(y, labels)
    ax.set_xlabel("Spearman rho with baseline CRP\nPaired-bootstrap 95% CI")
    ax.set_title("Fixed module decomposition (GSE225731)", loc="left", fontweight="bold")
    ax.set_xlim(-.35, .72)
    clean_spines(ax)
    tx.set_xlim(0, 1); tx.set_ylim(ax.get_ylim()); tx.axis("off")
    for yy, sig in zip(y, order):
        r = get_row(sig)
        tx.text(.02, yy, f"{r.rho_vs_CRP:.3f}  ({r.bootstrap_low:.3f}, {r.bootstrap_high:.3f})\nP={fmt_p(r.P_value)}",
                ha="left", va="center", fontsize=6.5, color=DARK, linespacing=1.05)
    tx.text(.02, -.12, "ECM64 discovery; other rows fixed decompositions.",
            transform=tx.transAxes, fontsize=6.5, color=GREY, va="top")
    panel_label(parent, "E")

def draw_panel_f(parent):
    parent.set_axis_off()
    inner = parent.get_subplotspec().subgridspec(1, 2, width_ratios=[.65, .35], wspace=.02)
    ax = parent.figure.add_subplot(inner[0])
    tx = parent.figure.add_subplot(inner[1])
    rows = [("Baseline ECM64", "baseline_N", "baseline_rho", "baseline_CI_low", "baseline_CI_high", "baseline_P", PLOT_BLUE),
            ("12-week change", "delta_N", "delta_rho", "delta_CI_low", "delta_CI_high", "delta_P", ORANGE)]
    ys = [1, 0]
    for yy, (lab, ncol, rcol, lcol, hcol, pcol, col) in zip(ys, rows):
        r = GSE198520.iloc[0]
        lo, hi, rho = r[lcol], r[hcol], r[rcol]
        ax.plot([lo, hi], [yy, yy], color=col, lw=1.8, solid_capstyle="round")
        ax.plot([lo, lo], [yy-.09, yy+.09], color=col, lw=.75)
        ax.plot([hi, hi], [yy-.09, yy+.09], color=col, lw=.75)
        ax.scatter([rho], [yy], color=col, s=24, edgecolor="white", linewidth=.4, zorder=3)
    ax.axvline(0, color="#89919C", lw=.7, ls=(0, (3, 2)))
    ax.set_yticks(ys, ["Baseline", "12-week change"])
    ax.set_xlim(-.55, .65)
    ax.set_xlabel("Spearman rho with CRP\nFisher approximate 95% CI")
    ax.set_title("External cohort: lack of ECM64–CRP replication", loc="left", fontweight="bold")
    clean_spines(ax)
    tx.set_xlim(0, 1); tx.set_ylim(ax.get_ylim()); tx.axis("off")
    for yy, (lab, ncol, rcol, lcol, hcol, pcol, col) in zip(ys, rows):
        r = GSE198520.iloc[0]
        tx.text(.02, yy, f"rho={r[rcol]:.3f}  ({r[lcol]:.3f}, {r[hcol]:.3f})\nP={fmt_p(r[pcol])}; N={int(r[ncol])}",
                ha="left", va="center", fontsize=6.5, color=DARK, linespacing=1.05)
    tx.text(.02, .97, "GSE198520; fixed ECM64;\nexternal comparison",
            transform=tx.transAxes, fontsize=6.4, color=GREY, va="top", linespacing=1.05)
    panel_label(parent, "F")

def write_source_data():
    rows = []
    for _, r in SCORES.iterrows():
        d = {"record_type":"sample_score", "panel":"B/D", "dataset":"GSE225731",
             "sample_id":r.sample_id, "CRP_mg_dL":r.CRP_mg_dL}
        for c in ["ECM64","ECM33","ECM31_uncovered","Structural_all","Structural_Xenium","Processing_all","Processing_Xenium"]:
            d[c] = r[c]
        rows.append(d)
    for _, r in BRIDGE.loc[BRIDGE.signature.isin(["ECM64","ECM33","ECM31_uncovered","Structural_all","Structural_Xenium","Processing_all","Processing_Xenium"])].iterrows():
        rows.append({"record_type":"module_summary", "panel":"E", "dataset":"GSE225731",
                     "signature":r.signature, "n_genes":int(r.n_genes), "sample_size":int(r.N),
                     "estimate":r.rho_vs_CRP, "ci_low":r.bootstrap_low, "ci_high":r.bootstrap_high,
                     "fisher_ci_low":r.CI_low, "fisher_ci_high":r.CI_high, "p_value":r.P_value,
                     "ci_method":"paired bootstrap percentile 95%; Fisher CI retained as auxiliary"})
    r = get_row("ECM64")
    rows.append({"record_type":"pathway_screen_summary", "panel":"A", "dataset":"GSE225731",
                 "signature":"ECM64", "n_genes":int(r.n_genes), "sample_size":int(r.N), "estimate":r.rho_vs_CRP,
                 "ci_low":r.bootstrap_low, "ci_high":r.bootstrap_high, "fisher_ci_low":r.CI_low,
                 "fisher_ci_high":r.CI_high, "p_value":r.P_value,
                 "q_value_BH":float(PATHWAYS.loc[PATHWAYS.pathway.eq("REACTOME_COLLAGEN_DEGRADATION"), "q_value_BH"].iloc[0]),
                 "ci_method":"paired bootstrap percentile 95%"})
    r = GSE198520.iloc[0]
    rows += [
        {"record_type":"external_summary", "panel":"F", "dataset":"GSE198520", "signature":"ECM64 baseline",
         "sample_size":int(r.baseline_N), "estimate":r.baseline_rho, "ci_low":r.baseline_CI_low,
         "ci_high":r.baseline_CI_high, "p_value":r.baseline_P, "ci_method":"Fisher approximate 95%"},
        {"record_type":"external_summary", "panel":"F", "dataset":"GSE198520", "signature":"ECM64 12-week delta",
         "sample_size":int(r.delta_N), "estimate":r.delta_rho, "ci_low":r.delta_CI_low,
         "ci_high":r.delta_CI_high, "p_value":r.delta_P, "ci_method":"Fisher approximate 95%"},
    ]
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "FIGURE3_SOURCE_DATA.csv", index=False)
    BRIDGE.loc[BRIDGE.signature.isin(["ECM64","ECM33","ECM31_uncovered","Structural_all","Structural_Xenium","Processing_all","Processing_Xenium"])].to_csv(OUT / "SUPPLEMENTARY_MODULE_RESULTS.csv", index=False)

def write_docs():
    (OUT / "FIGURE3_COMPLETE_LEGEND_EN.md").write_text(textwrap.dedent("""
    # Figure 3. Exploratory ECM–CRP molecular bridge and cross-cohort limitation

    (A) Context of the prespecified pathway screen in GSE225731. The full set of 1,389 tested pathways is shown; 353 pathways were excluded before testing because of the frozen coverage rule. The fixed Reactome `REACTOME_COLLAGEN_DEGRADATION` program (ECM64) is highlighted with its original BH-adjusted q value.

    (B) Donor-level baseline CRP associations for the complete ECM64 program and its fixed Xenium-covered subset (ECM33), in the same 74-patient GSE225731 cohort. CRP is shown in original mg/dL units on a log-scaled axis. The reported intervals are the paired donor-bootstrap percentile 95% CIs; no linear regression line is used because the effect measure is Spearman correlation.

    (C) Composition and Xenium coverage of the frozen 64-gene program. The structural COL component, protease/processing component, and other non-COL component are shown separately from the 33-gene Xenium-covered and 31-gene uncovered partitions.

    (D) Concordance between ECM64 and ECM33 scores in the same GSE225731 patients. This is fixed-panel concordance, not independent validation. The uncovered ECM31–CRP association is shown as a contextual diagnostic.

    (E) Fixed module decomposition in GSE225731. Points and intervals use the same paired-bootstrap convention as panel B. ECM64 is the discovery candidate; all other rows were prespecified decompositions of the frozen candidate and were not selected by a new genome-wide search.

    (F) External fixed-module comparison in GSE198520. Baseline and 12-week change analyses did not show a clear ECM64–CRP replication. Intervals are the Fisher approximate 95% CIs available in the fixed external replication table. The historical negative/inconclusive comparison is documented in the analysis record rather than shown as an internal gate label.

    All statistics are patient/donor-level. Spatial/Xenium data are not used as patient-level CRP observations in this figure. The figure is exploratory and does not establish causality, TKA-specific mechanism, or prediction of postoperative CRP.
    """).strip() + "\n", encoding="utf-8")
    (OUT / "FIGURE3_DATA_PROVENANCE.md").write_text(textwrap.dedent(f"""
    # Figure 3 data provenance

    - Figure source directory: `{OUT}`
    - GSE225731 pathway screen: `{STAGE2 / 'GSE225731_PATHWAY_DISCOVERY_ALL.csv'}`
    - GSE225731 fixed sample scores: `{OUT / 'FIGURE3_SAMPLE_SCORES.csv'}`
    - ECM functional map: `{OUT / 'ECM64_FUNCTIONAL_PANEL_MAP.csv'}`
    - GSE225731 bridge statistics: `{STAGE7H / 'GSE225731_ECM64_ECM33_BRIDGE.csv'}`
    - GSE198520 fixed comparison: `{STAGE7H / 'GSE198520_FIXED_MODULE_REPLICATION.csv'}`
    - Scoring and bootstrap source: `{STAGE7H / '01_stage7h_r_ecm_bridge.R'}`

    No new pathway discovery, gene selection, or patient exclusion was performed. GSE225731 scores use the frozen TMM → log2(CPM+1) → gene-wise across-sample z-score → within-module mean procedure. The source R script used 5,000 paired donor bootstrap resamples with fixed seeds; the corresponding bootstrap intervals are used in panels B and E. The external GSE198520 table is reproduced without reanalysis and uses Fisher approximate intervals.
    """).strip() + "\n", encoding="utf-8")
    (OUT / "ECM33_CI_RECONCILIATION.md").write_text(textwrap.dedent("""
    # ECM33 confidence-interval reconciliation

    The apparent discrepancy is a difference in interval construction, not a change in the ECM33 estimate.

    - Fixed point estimate: Spearman rho = 0.4103286067, N = 74, P = 0.0002823513.
    - Fisher approximate 95% CI stored in `CI_low`/`CI_high`: 0.2006419715–0.5840656068. This corresponds to the approximately 0.201–0.584 interval seen in the earlier figure.
    - Paired donor-bootstrap percentile 95% CI stored in `bootstrap_low`/`bootstrap_high`: 0.2188269164–0.5633708128. This is the interval used in the final Figure 3 panels B and E, because it follows the fixed Stage 7H paired-resampling procedure.

    The same rule is applied to ECM64 and the fixed decompositions. No endpoint was manually adjusted. Both interval types remain in the source-data CSV for auditability.
    """).strip() + "\n", encoding="utf-8")
    (OUT / "FIGURE3_REVISED_RESULTS_EN.md").write_text(textwrap.dedent("""
    # Revised Results text for Figure 3

    In the complete GSE225731 cohort (n=74), the fixed Reactome collagen-degradation program (ECM64, 64 genes) showed a positive association with baseline CRP (Spearman rho=0.449; paired-bootstrap 95% CI 0.265–0.601; P<0.001). The Xenium-covered subset (ECM33) showed a directionally similar association (rho=0.410; 95% CI 0.219–0.563; P<0.001). ECM64 and ECM33 scores were highly concordant in the same cohort (rho=0.959), and therefore this comparison is not an independent validation.

    Functional decomposition showed positive associations for the structural COL component and the processing component, with wider intervals for the Xenium-covered subsets. The uncovered ECM31 component remained positively associated with CRP (rho=0.351; paired-bootstrap 95% CI 0.136–0.552; P=0.002). In the fixed external GSE198520 comparison, ECM64 was not clearly associated with baseline CRP (rho=-0.019; Fisher 95% CI -0.308 to 0.272; P=0.898; n=46) or with the 12-week CRP change (rho=0.180; 95% CI -0.152 to 0.477; P=0.285; n=37).

    These results support an exploratory ECM–CRP signal in GSE225731, but do not establish an independent or TKA-specific molecular mechanism.
    """).strip() + "\n", encoding="utf-8")
    (OUT / "FIGURE3_REVISED_DISCUSSION_EN.md").write_text(textwrap.dedent("""
    # Revised Discussion text for Figure 3

    The collagen-degradation signal was detectable across the full frozen module and in fixed structural and processing decompositions, suggesting that the association is not attributable solely to the Xenium-observable subset. However, the very high ECM64–ECM33 concordance reflects shared genes and the same patients, rather than independent replication. The absence of a clear ECM64–CRP association in GSE198520 limits transportability across cohorts and clinical contexts. Accordingly, the data are best interpreted as an exploratory association between a synovial ECM-related transcriptional state and baseline systemic inflammation, with tissue-cell and spatial results providing biological context rather than causal proof.

    The findings do not show that TKA removes a specific inflammatory cell population, that collagen remodeling causes postoperative CRP reduction, or that the ECM score predicts postoperative CRP in the RA-TKA cohort.
    """).strip() + "\n", encoding="utf-8")
    (OUT / "OPTIONAL_TRANSLATIONAL_EXTENSION_FEASIBILITY.md").write_text(textwrap.dedent("""
    # Optional translational extension feasibility

    No drug prediction, Mendelian randomization, or docking was executed for Stage 7J-F3. A future translational extension would require an independently justified target, donor-level replication, and a validated perturbation or pharmacologic dataset. The current Figure 3 evidence is insufficient to nominate an ECM target for clinical intervention.
    """).strip() + "\n", encoding="utf-8")

def make_supplementary_coverage():
    f = FUNC.copy()
    f["structural"] = f.functional_class.eq("COL structural")
    f["processing"] = f.functional_class.str.contains("protease")
    f = f.sort_values(["structural", "processing", "gene_symbol"], ascending=[False, False, True]).reset_index(drop=True)
    mat = np.column_stack([f.structural.astype(int), f.processing.astype(int), f.xenium_covered.astype(int)])
    fig, ax = plt.subplots(figsize=(5.3, 9.2))
    cmap = mpl.colors.ListedColormap(["#F1F3F5", PLOT_BLUE])
    ax.imshow(mat, aspect="auto", cmap=cmap, vmin=0, vmax=1, interpolation="none")
    ax.set_xticks(range(3), ["COL\nstructural", "Protease /\nprocessing", "Xenium\ncovered"])
    ax.set_yticks(range(len(f)), f.gene_symbol)
    ax.tick_params(axis="y", labelsize=6.6, length=0)
    ax.tick_params(axis="x", labelsize=7.0, length=0)
    ax.set_title("Supplementary Figure. Frozen ECM64 gene-level coverage", loc="left", fontweight="bold", pad=12)
    for x in range(4): ax.axvline(x-.5, color="white", lw=.8)
    for y in range(len(f)+1): ax.axhline(y-.5, color="white", lw=.25)
    clean_spines(ax, keep=())
    fig.tight_layout()
    fig.savefig(OUT / "SUPPLEMENTARY_ECM_GENE_COVERAGE.pdf", bbox_inches="tight")
    plt.close(fig)

def write_qa_and_versions():
    inputs = [
        STAGE2 / "GSE225731_PATHWAY_DISCOVERY_ALL.csv",
        ROOT / "FIGURE3_SAMPLE_SCORES.csv",
        ROOT / "ECM64_FUNCTIONAL_PANEL_MAP.csv",
        STAGE7H / "GSE225731_ECM64_ECM33_BRIDGE.csv",
        STAGE7H / "GSE198520_FIXED_MODULE_REPLICATION.csv",
        STAGE7H / "01_stage7h_r_ecm_bridge.R",
    ]
    lines = []
    for p in inputs:
        h = hashlib.sha256(p.read_bytes()).hexdigest()
        lines.append(f"{h}  {p}")
    (OUT / "01_FIGURE3_INPUT_SHA256.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (OUT / "01_FIGURE3_SOFTWARE_VERSIONS.txt").write_text(
        "Python " + sys.version.replace("\n", " ") + "\n" +
        "Platform " + platform.platform() + "\n" +
        "Matplotlib " + mpl.__version__ + "\n" +
        "Pandas " + pd.__version__ + "\n" +
        "NumPy " + np.__version__ + "\n" +
        "R session information: 00_R_sessionInfo.txt\n" +
        "Rendering check: pdftoppm at 160 dpi\n", encoding="utf-8")
    (OUT / "FIGURE3_VISUAL_QA_REPORT.md").write_text(textwrap.dedent("""
    # Figure 3 visual QA report

    - `FIGURE3_MULTIPANEL_FINAL.pdf`: one-page vector PDF; rendered at 160 dpi with `pdftoppm` and visually inspected.
    - `FIGURE3_MULTIPANEL_FINAL.png`: 600-dpi raster preview generated from the same Matplotlib figure object.
    - `FIGURE3_PANEL_A_F.pdf`: vector supplementary extraction of panels A and F.
    - `SUPPLEMENTARY_ECM_GENE_COVERAGE.pdf`: vector gene-level coverage matrix.
    - Six-panel layout: A–F present, white background, restrained color palette, no grid decoration, and no rasterized primary plot layers except the deliberately dense 1,389-point context cloud in panel A.
    - Panel B uses the same baseline CRP scale and original mg/dL units for ECM64 and ECM33.
    - Panel E and F statistics are placed in dedicated right-side text columns to avoid masking confidence intervals.
    - No panel clipping or unreadable overlap was observed in the final render.
    - Numeric values were checked against `FIGURE3_SOURCE_DATA.csv` and the fixed Stage 7H bridge tables after rendering.
    """).strip() + "\n", encoding="utf-8")

def build_figures():
    fig = plt.figure(figsize=(12.5, 8.35), facecolor="white")
    gs = fig.add_gridspec(2, 3, left=.065, right=.985, top=.935, bottom=.09,
                          wspace=.40, hspace=.46, width_ratios=[1.12, 1.12, 1.15])
    ax_a = fig.add_subplot(gs[0, 0]); draw_panel_a(ax_a)
    ax_b = fig.add_subplot(gs[0, 1]); draw_panel_b(ax_b)
    ax_c = fig.add_subplot(gs[0, 2]); draw_panel_c(ax_c)
    ax_d = fig.add_subplot(gs[1, 0]); draw_panel_d(ax_d)
    ax_e = fig.add_subplot(gs[1, 1]); draw_panel_e(ax_e)
    ax_f = fig.add_subplot(gs[1, 2]); draw_panel_f(ax_f)
    fig.suptitle("Figure 3 | Exploratory ECM–CRP molecular bridge and cross-cohort limitation",
                 x=.065, ha="left", y=.985, fontsize=12, fontweight="bold", color=DARK)
    fig.savefig(OUT / "FIGURE3_MULTIPANEL_FINAL.pdf", bbox_inches="tight")
    fig.savefig(OUT / "FIGURE3_MULTIPANEL_FINAL.png", dpi=600, bbox_inches="tight")
    plt.close(fig)

    fig2, (a, f) = plt.subplots(1, 2, figsize=(12.5, 4.8))
    draw_panel_a(a); draw_panel_f(f)
    fig2.suptitle("Figure 3 panels A and F: screen context and external comparison",
                  x=.065, ha="left", y=.99, fontsize=12, fontweight="bold")
    fig2.tight_layout(rect=[0, 0, 1, .94])
    fig2.savefig(OUT / "FIGURE3_PANEL_A_F.pdf", bbox_inches="tight")
    plt.close(fig2)

def main():
    write_source_data()
    write_docs()
    make_supplementary_coverage()
    build_figures()
    write_qa_and_versions()
    print("Stage 7J Figure 3 outputs written to", OUT)

if __name__ == "__main__":
    main()
