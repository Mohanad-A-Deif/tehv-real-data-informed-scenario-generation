"""Robustness summary figure (fig_22) from the stored R2 tables.
(a) Feasibility rate of each scenario group under joint +/-50% perturbation of all endpoint coefficients.
(b) Overlap with the reference top-50 feasible set under random perturbation of the objective weights."""
import pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.family": "serif", "font.size": 13})
G = ["baseline", "hydrodynamic_only", "biology_only", "full_data_informed"]
LAB = ["Baseline", "Intermediate-\nprior", "Evidence-centered\nwide-prior", "Evidence-centered\nnarrow-prior"]
COL = ["#7d8187", "#9a4a54", "#4e735f", "#4a5468"]

gl = pd.read_csv("tables/R2_T6_coef_global.csv"); gl = gl[gl.level.astype(str) == "0.5"]
wr = pd.read_csv("tables/R2_T8_weight_random.csv")

fig, ax = plt.subplots(1, 2, figsize=(20, 7.65), dpi=100)
box = dict(widths=0.45, whis=1.5, showfliers=False, patch_artist=True,
           medianprops=dict(color="black", lw=1.5), whiskerprops=dict(lw=1.5), capprops=dict(lw=1.5))
b = ax[0].boxplot([100 * gl["feas_" + g] for g in G], **box)
for p, c in zip(b["boxes"], COL): p.set(facecolor=c, edgecolor="black", lw=1.5)
ax[0].set_xticks(range(1, 5)); ax[0].set_xticklabels(LAB, fontsize=12)
ax[0].set_ylabel("Feasibility rate (%)", fontsize=15)
ax[0].text(0.03, 0.95, "(a)", transform=ax[0].transAxes, fontweight="bold", fontsize=15, va="top")

b = ax[1].boxplot([100 * wr[wr.level.astype(str) == k].top50_overlap for k in ["0.2", "0.5", "dirichlet"]], **box)
for p in b["boxes"]: p.set(facecolor=COL[3], edgecolor="black", lw=1.5)
ax[1].set_xticks(range(1, 4)); ax[1].set_xticklabels(["Weights\n±20%", "Weights\n±50%", "Random\n(Dirichlet)"], fontsize=15)
ax[1].set_ylabel("Overlap with reference top-50 set (%)", fontsize=15)
ax[1].text(0.03, 0.05, "(b)", transform=ax[1].transAxes, fontweight="bold", fontsize=15)

for a in ax:
    a.set_ylim(0, 100); a.grid(axis="y", ls=":", color="0.8"); a.tick_params(direction="in", labelsize=15)
    for s in a.spines.values(): s.set_linewidth(1.5)
fig.tight_layout()
fig.savefig("fig_22_robustness_summary.png", dpi=345)
print("saved fig_22_robustness_summary.png")
