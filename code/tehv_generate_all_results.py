# -*- coding: utf-8 -*-
r"""
TEHV real-data-informed result generator
=======================================

Purpose
-------
Generate manuscript-ready exploratory results, data audits, calibration tables,
uncertainty summaries, statistical comparisons, and PNG figures for the TEHV
bioresorbable scaffold manuscript.

Important scientific note
-------------------------
This script generates REAL-DATA-INFORMED computational/scenario results from the
extracted dataset. It does NOT create clinical validation and should not be used
to claim patient-specific prescription validity.

How to run
----------
1) Install dependencies:
   pip install -r requirements.txt

2) From the repository root, run:
   python code/tehv_generate_all_results.py

3) Optional: set a custom output directory:
   Windows PowerShell: $env:TEHV_OUTPUT_DIR="C:\path\to\output"
   Linux/macOS: export TEHV_OUTPUT_DIR=/path/to/output

Outputs
-------
Generated outputs are saved as PNG/CSV/XLSX only. No PDF is created.
The manuscript text is intentionally not included in this repository.
"""

from __future__ import annotations

import os
import re
import json
import math
import shutil
import warnings
from pathlib import Path
from datetime import datetime
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

warnings.filterwarnings("ignore")

# -----------------------------------------------------------------------------
# 0. USER CONFIGURATION
# -----------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent if SCRIPT_DIR.name == "code" else SCRIPT_DIR

# By default, generated files are written outside the curated results folder so
# the repository contents remain clean. Override using the TEHV_OUTPUT_DIR env var.
OUTPUT_DIR = Path(os.environ.get("TEHV_OUTPUT_DIR", REPO_ROOT / "generated_outputs")).resolve()

DATASET_FILENAME = "TEHV_real_data_extraction_dataset_v2.xlsx"
RESULTS_FILENAME = "TEHV_real_data_informed_generated_results_v1.xlsx"

RANDOM_SEED = 42
N_CANDIDATES_PER_GROUP = 1200
N_BOOTSTRAP = 800
SAVE_DPI = 600
SAVE_PNG_ONLY = True
RUN_RF_PERMUTATION_IMPORTANCE = False  # Set True only if you want slower ML importance analysis

# Style mode options:
#   "dark_academic" = dark academic colors requested by user
#   "grayscale"     = strict Nature/IEEE grayscale style from the uploaded prompt
STYLE_MODE = "dark_academic"

# -----------------------------------------------------------------------------
# 1. PATHS AND OUTPUT FOLDERS
# -----------------------------------------------------------------------------
def make_dirs() -> dict[str, Path]:
    paths = {
        "root": OUTPUT_DIR,
        "figures": OUTPUT_DIR / "figures_png",
        "tables": OUTPUT_DIR / "tables_csv",
        "excel": OUTPUT_DIR / "excel_outputs",
        "logs": OUTPUT_DIR / "logs",
    }
    for p in paths.values():
        p.mkdir(parents=True, exist_ok=True)
    return paths

OUT = make_dirs()


def find_input_file(filename: str) -> Path | None:
    candidates = [
        OUTPUT_DIR / filename,
        SCRIPT_DIR / filename,
        REPO_ROOT / filename,
        REPO_ROOT / "data" / filename,
        REPO_ROOT / "results" / "excel_outputs" / filename,
        Path.cwd() / filename,
        Path.cwd() / "data" / filename,
        Path.cwd() / "results" / "excel_outputs" / filename,
    ]
    # Also search recursively in the repository, output, script, and working folders.
    for folder in [REPO_ROOT, OUTPUT_DIR, SCRIPT_DIR, Path.cwd()]:
        try:
            for p in folder.glob(f"**/{filename}"):
                candidates.append(p)
        except Exception:
            pass
    seen = set()
    for p in candidates:
        try:
            rp = p.resolve()
        except Exception:
            rp = p
        if str(rp) in seen:
            continue
        seen.add(str(rp))
        if p.exists():
            return p
    return None

DATASET_PATH = find_input_file(DATASET_FILENAME)
RESULTS_PATH = find_input_file(RESULTS_FILENAME)

# -----------------------------------------------------------------------------
# 2. PLOTTING STYLE
# -----------------------------------------------------------------------------
DARK_PALETTE = {
    "navy": "#1B263B",
    "burgundy": "#7A1E2C",
    "forest": "#1B4332",
    "charcoal": "#222222",
    "slate": "#415A77",
    "plum": "#4A273F",
    "brown": "#5C4033",
    "gray": "#5F6368",
    "lightgray": "#D9D9D9",
}

GRAY_PALETTE = {
    "navy": "#000000",
    "burgundy": "#333333",
    "forest": "#555555",
    "charcoal": "#111111",
    "slate": "#777777",
    "plum": "#444444",
    "brown": "#666666",
    "gray": "#888888",
    "lightgray": "#D9D9D9",
}

COL = DARK_PALETTE if STYLE_MODE == "dark_academic" else GRAY_PALETTE


def set_academic_style() -> None:
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 11,
        "axes.edgecolor": "black",
        "axes.linewidth": 0.8,
        "axes.labelsize": 11,
        "axes.titlesize": 11,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.major.size": 4,
        "ytick.major.size": 4,
        "grid.alpha": 0.25,
        "figure.dpi": SAVE_DPI,
        "savefig.dpi": SAVE_DPI,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "legend.frameon": False,
        "axes.spines.top": True,
        "axes.spines.right": True,
    })

set_academic_style()


def save_figure(fig: plt.Figure, name: str) -> Path:
    """Save PNG only. No PDF output is created."""
    safe = re.sub(r"[^A-Za-z0-9_\-]+", "_", name).strip("_")
    out = OUT["figures"] / f"{safe}.png"
    fig.tight_layout()
    fig.savefig(out, dpi=SAVE_DPI, bbox_inches="tight", format="png")
    plt.close(fig)
    return out


def subplot_label(ax, label: str) -> None:
    ax.text(
        0.02,
        0.96,
        label,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=11,
        fontweight="bold",
        color="black",
    )


def add_metric_box(ax, text: str, loc=(0.04, 0.08)) -> None:
    ax.text(
        loc[0],
        loc[1],
        text,
        transform=ax.transAxes,
        fontsize=9,
        va="bottom",
        ha="left",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.82, boxstyle="round,pad=0.25"),
    )

# -----------------------------------------------------------------------------
# 3. DATA READING AND CLEANING HELPERS
# -----------------------------------------------------------------------------
def clean_colname(x) -> str:
    x = str(x).strip()
    x = re.sub(r"\s+", "_", x)
    x = x.replace("/", "_per_")
    x = re.sub(r"[^A-Za-z0-9_\-]+", "", x)
    return x.strip("_") or "col"


def read_workbook(path: Path | None) -> dict[str, pd.DataFrame]:
    if path is None or not path.exists():
        return {}
    xl = pd.ExcelFile(path)
    return {s: pd.read_excel(path, sheet_name=s) for s in xl.sheet_names}


def normalize_embedded_header(df: pd.DataFrame) -> pd.DataFrame:
    """Handle sheets where row 1 contains the actual header after a title row."""
    if df is None or df.empty:
        return pd.DataFrame()
    d = df.copy()
    cols = [str(c) for c in d.columns]
    unnamed_frac = sum(c.startswith("Unnamed") for c in cols) / max(1, len(cols))
    first_col_title_like = len(cols) > 0 and not cols[0].startswith("Unnamed") and d.iloc[:, 0].isna().sum() > 0
    if unnamed_frac > 0.3 or first_col_title_like:
        best_idx = None
        best_nonnull = 0
        for i in range(min(5, len(d))):
            vals = d.iloc[i].tolist()
            nonnull = sum(pd.notna(v) and str(v).strip() != "" for v in vals)
            if nonnull > best_nonnull:
                best_nonnull = nonnull
                best_idx = i
        if best_idx is not None and best_nonnull >= 2:
            new_cols = [clean_colname(v) for v in d.iloc[best_idx].tolist()]
            d = d.iloc[best_idx + 1:].copy()
            d.columns = new_cols
    else:
        d.columns = [clean_colname(c) for c in d.columns]
    d = d.dropna(how="all").reset_index(drop=True)
    return d


def numericize(series: pd.Series) -> pd.Series:
    if series is None:
        return pd.Series(dtype=float)
    if pd.api.types.is_numeric_dtype(series):
        return pd.to_numeric(series, errors="coerce")
    s = series.astype(str).str.replace("±", " ", regex=False)
    s = s.str.replace(",", "", regex=False)
    s = s.str.extract(r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", expand=False)
    return pd.to_numeric(s, errors="coerce")


def clean_all_workbook_sheets(raw: dict[str, pd.DataFrame]) -> dict[str, pd.DataFrame]:
    return {name: normalize_embedded_header(df) for name, df in raw.items()}


def safe_to_csv(df: pd.DataFrame, name: str) -> Path:
    out = OUT["tables"] / f"{name}.csv"
    df.to_csv(out, index=False, encoding="utf-8-sig")
    return out


def write_log(payload: dict) -> None:
    with open(OUT["logs"] / "run_log.json", "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False, default=str)

# -----------------------------------------------------------------------------
# 4. BASIC STATISTICS HELPERS
# -----------------------------------------------------------------------------
def bootstrap_ci(x, func=np.mean, n_boot=N_BOOTSTRAP, ci=95, seed=RANDOM_SEED):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) < 2:
        return np.nan, np.nan, np.nan
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_boot):
        sample = rng.choice(x, size=len(x), replace=True)
        vals.append(func(sample))
    low = np.percentile(vals, (100 - ci) / 2)
    high = np.percentile(vals, 100 - (100 - ci) / 2)
    return func(x), low, high


def cliffs_delta(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if len(x) == 0 or len(y) == 0:
        return np.nan
    # Efficient enough for this size after subsampling if needed
    max_pairs = 2_000_000
    if len(x) * len(y) > max_pairs:
        rng = np.random.default_rng(RANDOM_SEED)
        x = rng.choice(x, size=min(len(x), 1000), replace=False)
        y = rng.choice(y, size=min(len(y), 1000), replace=False)
    gt = sum((xi > y).sum() for xi in x)
    lt = sum((xi < y).sum() for xi in x)
    return (gt - lt) / (len(x) * len(y))


def try_scipy_tests():
    try:
        from scipy import stats
        return stats
    except Exception:
        return None

STATS = try_scipy_tests()

# -----------------------------------------------------------------------------
# 5. LOAD DATA
# -----------------------------------------------------------------------------
raw_dataset = read_workbook(DATASET_PATH)
raw_results = read_workbook(RESULTS_PATH)
dataset = clean_all_workbook_sheets(raw_dataset)
results = clean_all_workbook_sheets(raw_results)

if not dataset and not results:
    raise FileNotFoundError(
        f"Could not find {DATASET_FILENAME} or {RESULTS_FILENAME}.\n"
        f"Put them in the repository data/ folder, {OUTPUT_DIR}, {SCRIPT_DIR}, or the current working directory."
    )

# Convenience handles
extracted = dataset.get("Extracted_Numeric_Data", pd.DataFrame()).copy()
source_log = dataset.get("Source_Log", pd.DataFrame()).copy()
scaffold_response = dataset.get("Scaffold_Response", pd.DataFrame()).copy()
parameter_map = dataset.get("Parameter_Map", pd.DataFrame()).copy()
clinical_echo = dataset.get("Clinical_Echo_Cohort", pd.DataFrame()).copy()

# Result-pack sheets, if available
integration = results.get("Integration_Timecourse", pd.DataFrame()).copy()
regional = results.get("Regional_Colonization", pd.DataFrame()).copy()
polymer = results.get("Polymer_Degradation", pd.DataFrame()).copy()
hydro = results.get("Hydrodynamic_Targets", pd.DataFrame()).copy()
mc_summary = results.get("Monte_Carlo_Summary", pd.DataFrame()).copy()
trend_fits = results.get("Trend_Fits", pd.DataFrame()).copy()
charts_data = results.get("Charts_Data", pd.DataFrame()).copy()

# Numeric conversions
if not extracted.empty:
    for c in ["mean", "sd", "n"]:
        if c in extracted.columns:
            extracted[c] = numericize(extracted[c])

for df in [integration, regional, polymer, hydro, mc_summary, trend_fits, charts_data]:
    if not df.empty:
        for c in df.columns:
            # Keep text columns as text; numericize if at least 40% values look numeric
            ser_num = numericize(df[c])
            if ser_num.notna().sum() >= max(2, int(0.4 * len(df))):
                df[c] = ser_num

# -----------------------------------------------------------------------------
# 6. DATA AUDIT TABLES AND FIGURES
# -----------------------------------------------------------------------------
def build_data_audit_tables():
    rows = []
    for book_name, book in [("dataset", dataset), ("results", results)]:
        for sheet, df in book.items():
            rows.append({
                "workbook": book_name,
                "sheet": sheet,
                "n_rows": int(df.shape[0]),
                "n_cols": int(df.shape[1]),
                "nonmissing_cells": int(df.notna().sum().sum()),
                "missing_cells": int(df.isna().sum().sum()),
                "missing_fraction": float(df.isna().sum().sum() / max(1, df.shape[0] * df.shape[1])),
            })
    audit = pd.DataFrame(rows)
    safe_to_csv(audit, "table_01_data_audit_by_sheet")

    if not extracted.empty:
        metric_summary = (
            extracted.groupby(["source_id", "metric_group", "unit"], dropna=False)
            .agg(n_metrics=("metric_id", "count"), mean_of_means=("mean", "mean"), sd_of_means=("mean", "std"))
            .reset_index()
            .sort_values(["source_id", "metric_group"])
        )
    else:
        metric_summary = pd.DataFrame()
    safe_to_csv(metric_summary, "table_02_extracted_numeric_metric_summary")

    if not source_log.empty:
        source_summary = (
            source_log.groupby(["domain", "source_type", "priority", "status"], dropna=False)
            .size()
            .reset_index(name="n_sources")
            .sort_values("n_sources", ascending=False)
        )
    else:
        source_summary = pd.DataFrame()
    safe_to_csv(source_summary, "table_03_source_summary")
    return audit, metric_summary, source_summary


def plot_data_audit(audit: pd.DataFrame, metric_summary: pd.DataFrame):
    if not audit.empty:
        fig, ax = plt.subplots(figsize=(7.2, 4.2))
        d = audit.sort_values("nonmissing_cells", ascending=True)
        ax.barh(d["sheet"], d["nonmissing_cells"], color=COL["navy"], alpha=0.92)
        ax.set_xlabel("Non-missing cells")
        ax.set_ylabel("Workbook sheet")
        ax.grid(True, axis="x", linestyle=":")
        subplot_label(ax, "(a)")
        save_figure(fig, "fig_01_data_audit_nonmissing_cells")

        fig, ax = plt.subplots(figsize=(7.2, 4.2))
        d = audit.sort_values("missing_fraction", ascending=True)
        ax.barh(d["sheet"], d["missing_fraction"], color=COL["burgundy"], alpha=0.88)
        ax.set_xlabel("Missing fraction")
        ax.set_ylabel("Workbook sheet")
        ax.set_xlim(0, max(1.0, d["missing_fraction"].max() * 1.1))
        ax.grid(True, axis="x", linestyle=":")
        subplot_label(ax, "(b)")
        save_figure(fig, "fig_02_data_audit_missing_fraction")

    if not metric_summary.empty:
        d = metric_summary.groupby("metric_group", dropna=False)["n_metrics"].sum().sort_values(ascending=True).tail(25)
        fig, ax = plt.subplots(figsize=(7.2, 5.2))
        ax.barh(d.index.astype(str), d.values, color=COL["forest"], alpha=0.92)
        ax.set_xlabel("Number of extracted numeric metrics")
        ax.set_ylabel("Metric group")
        ax.grid(True, axis="x", linestyle=":")
        subplot_label(ax, "(c)")
        save_figure(fig, "fig_03_extracted_metric_groups")

# -----------------------------------------------------------------------------
# 7. REAL-DATA-INFORMED MANUSCRIPT FIGURES
# -----------------------------------------------------------------------------
def plot_hydrodynamic_targets(hydro_df: pd.DataFrame):
    if hydro_df.empty or "metric" not in hydro_df.columns or "value" not in hydro_df.columns:
        return
    d = hydro_df.dropna(subset=["metric", "value"]).copy()
    if d.empty:
        return
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    x = np.arange(len(d))
    err = d["sd_or_variability"] if "sd_or_variability" in d.columns else None
    err = numericize(err) if err is not None else None
    yerr = err.values if err is not None and err.notna().sum() else None
    ax.bar(x, d["value"], yerr=yerr, capsize=3, color=COL["navy"], alpha=0.90, edgecolor="black", linewidth=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(d["metric"].astype(str), rotation=35, ha="right")
    ax.set_ylabel("Extracted value")
    ax.grid(True, axis="y", linestyle=":")
    subplot_label(ax, "(a)")
    save_figure(fig, "fig_04_real_hydrodynamic_targets")


def plot_regional_colonization(regional_df: pd.DataFrame):
    needed = {"month", "hinge_cells_mm2", "belly_cells_mm2", "tip_cells_mm2"}
    if regional_df.empty or not needed.issubset(set(regional_df.columns)):
        return
    d = regional_df.dropna(subset=["month"]).copy()
    d = d.sort_values("month")
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.plot(d["month"], d["hinge_cells_mm2"], marker="o", linestyle="-", color=COL["navy"], label="Hinge")
    ax.plot(d["month"], d["belly_cells_mm2"], marker="s", linestyle="--", color=COL["burgundy"], label="Belly")
    ax.plot(d["month"], d["tip_cells_mm2"], marker="^", linestyle=":", color=COL["forest"], label="Tip")
    ax.set_xlabel("Month")
    ax.set_ylabel("Cell density (cells/mm²)")
    ax.grid(True, linestyle=":")
    ax.legend(loc="best")
    subplot_label(ax, "(a)")
    if "hinge_to_belly_ratio" in d.columns:
        final = d.iloc[-1]
        add_metric_box(ax, f"Final hinge/belly = {final['hinge_to_belly_ratio']:.2f}")
    save_figure(fig, "fig_05_region_specific_colonization_timecourse")

    ratio_cols = [c for c in ["hinge_share_of_total", "hinge_to_belly_ratio", "hinge_to_tip_ratio"] if c in d.columns]
    if ratio_cols:
        fig, ax = plt.subplots(figsize=(6.2, 4.2))
        linestyles = ["-", "--", ":"]
        markers = ["o", "s", "^"]
        colors = [COL["navy"], COL["burgundy"], COL["forest"]]
        for i, c in enumerate(ratio_cols):
            ax.plot(d["month"], d[c], marker=markers[i], linestyle=linestyles[i], color=colors[i], label=c.replace("_", " "))
        ax.set_xlabel("Month")
        ax.set_ylabel("Regional imbalance metric")
        ax.grid(True, linestyle=":")
        ax.legend(loc="best")
        subplot_label(ax, "(b)")
        save_figure(fig, "fig_06_regional_colonization_ratios")


def plot_integration_timecourse(integration_df: pd.DataFrame, mc_df: pd.DataFrame):
    if integration_df.empty or "month" not in integration_df.columns:
        return
    d = integration_df.dropna(subset=["month"]).sort_values("month").copy()
    ycol = "integration_index_0_1" if "integration_index_0_1" in d.columns else None
    if ycol is None:
        return
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.plot(d["month"], d[ycol], marker="o", linestyle="-", color=COL["navy"], label="Data-derived index")

    # Add Monte Carlo CI if present
    if not mc_df.empty and {"derived_metric", "month", "p2_5", "median", "p97_5"}.issubset(mc_df.columns):
        m = mc_df[mc_df["derived_metric"].astype(str).str.contains("integration", case=False, na=False)].copy()
        m = m.dropna(subset=["month", "p2_5", "median", "p97_5"]).sort_values("month")
        if not m.empty:
            ax.fill_between(m["month"], m["p2_5"], m["p97_5"], color=COL["slate"], alpha=0.18, label="95% MC interval")
            ax.plot(m["month"], m["median"], marker="s", linestyle="--", color=COL["burgundy"], label="MC median")
    ax.set_xlabel("Month")
    ax.set_ylabel("Integration index (0–1)")
    ax.set_ylim(0, 1.05)
    ax.grid(True, linestyle=":")
    ax.legend(loc="best")
    subplot_label(ax, "(a)")
    save_figure(fig, "fig_07_data_informed_integration_timecourse")

    # ECM components panel
    cols = [c for c in ["HYP_ug_mg", "GAG_ug_mg", "global_collagen_area_percent", "ECM_index_0_1"] if c in d.columns]
    if cols:
        fig, ax = plt.subplots(figsize=(6.2, 4.2))
        linestyles = ["-", "--", ":", "-."]
        markers = ["o", "s", "^", "D"]
        colors = [COL["navy"], COL["burgundy"], COL["forest"], COL["plum"]]
        for i, c in enumerate(cols):
            vals = d[c]
            # Normalize non-index variables for common-scale display
            if not c.endswith("0_1"):
                denom = np.nanmax(vals)
                vals = vals / denom if np.isfinite(denom) and denom > 0 else vals
                label = c.replace("_", " ") + " (scaled)"
            else:
                label = c.replace("_", " ")
            ax.plot(d["month"], vals, marker=markers[i], linestyle=linestyles[i], color=colors[i], label=label)
        ax.set_xlabel("Month")
        ax.set_ylabel("Scaled ECM marker / index")
        ax.set_ylim(0, 1.15)
        ax.grid(True, linestyle=":")
        ax.legend(loc="best")
        subplot_label(ax, "(b)")
        save_figure(fig, "fig_08_ecm_marker_timecourse")


def plot_polymer_degradation(polymer_df: pd.DataFrame, trend_df: pd.DataFrame):
    if polymer_df.empty or "month" not in polymer_df.columns:
        return
    d = polymer_df.dropna(subset=["month"]).sort_values("month").copy()
    if not {"GPC_Mn_kg_mol", "GPC_Mw_kg_mol"}.issubset(d.columns):
        return
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.plot(d["month"], d["GPC_Mw_kg_mol"], marker="o", linestyle="-", color=COL["navy"], label="Mw")
    ax.plot(d["month"], d["GPC_Mn_kg_mol"], marker="s", linestyle="--", color=COL["burgundy"], label="Mn")
    ax.set_xlabel("Month")
    ax.set_ylabel("GPC molecular weight (kg/mol)")
    ax.grid(True, linestyle=":")
    ax.legend(loc="best")
    subplot_label(ax, "(a)")
    # Fit log-linear Mw decay
    x = d["month"].values.astype(float)
    y = d["GPC_Mw_kg_mol"].values.astype(float)
    mask = np.isfinite(x) & np.isfinite(y) & (y > 0)
    if mask.sum() >= 3:
        coef = np.polyfit(x[mask], np.log(y[mask]), 1)
        slope = coef[0]
        half_life = np.log(2) / abs(slope) if slope < 0 else np.inf
        yhat = np.exp(np.polyval(coef, x[mask]))
        ax.plot(x[mask], yhat, linestyle=":", color=COL["charcoal"], label="log-linear fit")
        add_metric_box(ax, f"Mw half-life ≈ {half_life:.1f} months")
    save_figure(fig, "fig_09_polymer_gpc_mn_mw_trends")

    if {"fiber_pre_um", "fiber_explant_um"}.issubset(d.columns):
        fig, ax = plt.subplots(figsize=(6.2, 4.2))
        ax.plot(d["month"], d["fiber_pre_um"], marker="o", linestyle="-", color=COL["forest"], label="Pre-implant")
        ax.plot(d["month"], d["fiber_explant_um"], marker="s", linestyle="--", color=COL["plum"], label="Explant")
        ax.set_xlabel("Month")
        ax.set_ylabel("Fiber diameter (µm)")
        ax.grid(True, linestyle=":")
        ax.legend(loc="best")
        subplot_label(ax, "(b)")
        save_figure(fig, "fig_10_fiber_diameter_change")


def plot_mc_distributions(mc_df: pd.DataFrame):
    if mc_df.empty or not {"derived_metric", "month", "mean", "sd", "p2_5", "median", "p97_5"}.issubset(mc_df.columns):
        return
    d = mc_df.dropna(subset=["month", "median"]).copy()
    if d.empty:
        return
    # Plot all MC interval summaries by derived metric
    for metric in d["derived_metric"].dropna().unique():
        m = d[d["derived_metric"] == metric].sort_values("month")
        fig, ax = plt.subplots(figsize=(6.2, 4.2))
        ax.fill_between(m["month"], m["p2_5"], m["p97_5"], color=COL["slate"], alpha=0.20, label="95% interval")
        ax.plot(m["month"], m["median"], marker="o", linestyle="-", color=COL["navy"], label="Median")
        ax.plot(m["month"], m["mean"], marker="s", linestyle="--", color=COL["burgundy"], label="Mean")
        ax.set_xlabel("Month")
        ax.set_ylabel(str(metric).replace("_", " "))
        ax.grid(True, linestyle=":")
        ax.legend(loc="best")
        subplot_label(ax, "(a)")
        save_figure(fig, f"fig_11_mc_interval_{metric}")

# -----------------------------------------------------------------------------
# 8. DATA-INFORMED SCENARIO DESIGN GENERATOR
# -----------------------------------------------------------------------------
def sigmoid(z):
    return 1 / (1 + np.exp(-z))


def derive_calibration_constants():
    constants = {}
    # Hydrodynamic targets from hydrodynamic table
    constants["pulmonary_gradient_target"] = 23.37
    constants["aortic_gradient_target"] = 32.75
    constants["regurgitation_target"] = 11.47
    constants["eoa_target"] = 2.15
    if not hydro.empty and "metric" in hydro.columns and "value" in hydro.columns:
        for _, row in hydro.iterrows():
            metric = str(row.get("metric", "")).lower()
            value = row.get("value", np.nan)
            if not np.isfinite(value):
                continue
            if "eoa" in metric:
                constants["eoa_target"] = float(value)
            if "gradient" in metric and "yacoub" in metric:
                constants["pulmonary_gradient_target"] = float(value)
            if "regurg" in metric and "yacoub" in metric:
                constants["regurgitation_target"] = float(value)

    # Polymer half-life from trend fits or fit direct
    constants["mw_half_life_months"] = 58.8
    if not trend_fits.empty:
        if "trend" in trend_fits.columns and "half_life_months" in trend_fits.columns:
            vals = trend_fits[trend_fits["trend"].astype(str).str.contains("GPC_Mw", case=False, na=False)]["half_life_months"].dropna()
            if len(vals):
                constants["mw_half_life_months"] = float(vals.iloc[0])
    if not polymer.empty and {"month", "GPC_Mw_kg_mol"}.issubset(polymer.columns):
        x = polymer["month"].values.astype(float)
        y = polymer["GPC_Mw_kg_mol"].values.astype(float)
        mask = np.isfinite(x) & np.isfinite(y) & (y > 0)
        if mask.sum() >= 3:
            coef = np.polyfit(x[mask], np.log(y[mask]), 1)
            slope = coef[0]
            if slope < 0:
                constants["mw_half_life_months"] = float(np.log(2) / abs(slope))

    # Integration target from integration or MC
    constants["integration_12m_target"] = 0.75
    if not integration.empty and {"month", "integration_index_0_1"}.issubset(integration.columns):
        d = integration.dropna(subset=["month", "integration_index_0_1"]).sort_values("month")
        if not d.empty:
            constants["integration_12m_target"] = float(d.iloc[-1]["integration_index_0_1"])
    if not mc_summary.empty and {"derived_metric", "month", "median"}.issubset(mc_summary.columns):
        m = mc_summary[mc_summary["derived_metric"].astype(str).str.contains("integration", case=False, na=False)]
        m = m.dropna(subset=["month", "median"]).sort_values("month")
        if not m.empty:
            constants["integration_12m_target"] = float(m.iloc[-1]["median"])

    # Region specificity
    constants["hinge_belly_ratio_12m"] = 4.41
    if not regional.empty and {"month", "hinge_to_belly_ratio"}.issubset(regional.columns):
        d = regional.dropna(subset=["month", "hinge_to_belly_ratio"]).sort_values("month")
        if not d.empty:
            constants["hinge_belly_ratio_12m"] = float(d.iloc[-1]["hinge_to_belly_ratio"])
    return constants

CAL = derive_calibration_constants()


def generate_candidate_scenarios(n=N_CANDIDATES_PER_GROUP, seed=RANDOM_SEED) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    groups = ["synthetic_only", "hydrodynamic_only", "biology_only", "full_data_informed"]
    frames = []
    for g in groups:
        # Prior changes by group
        if g == "synthetic_only":
            thickness = rng.uniform(0.10, 0.95, n)
            pore_index = rng.uniform(0.05, 0.95, n)
            degradation_half_life = rng.uniform(8, 40, n)  # old model allowed rapid resorption
            material_risk = rng.uniform(0.05, 0.95, n)
            morphology = rng.uniform(0.05, 0.95, n)
        elif g == "hydrodynamic_only":
            thickness = np.clip(rng.normal(0.55, 0.16, n), 0.18, 1.05)
            pore_index = rng.beta(2.0, 2.0, n)
            degradation_half_life = rng.uniform(10, 48, n)
            material_risk = rng.uniform(0.05, 0.95, n)
            morphology = rng.uniform(0.05, 0.95, n)
        elif g == "biology_only":
            thickness = np.clip(rng.normal(0.62, 0.18, n), 0.22, 1.12)
            pore_index = rng.beta(2.6, 1.7, n)
            degradation_half_life = np.clip(rng.normal(CAL["mw_half_life_months"], 11.0, n), 18, 96)
            material_risk = rng.beta(1.5, 4.0, n)
            morphology = rng.beta(2.8, 1.8, n)
        else:
            thickness = np.clip(rng.normal(0.62, 0.11, n), 0.32, 0.95)
            pore_index = rng.beta(2.4, 1.9, n)
            degradation_half_life = np.clip(rng.normal(CAL["mw_half_life_months"], 8.0, n), 30, 92)
            material_risk = rng.beta(1.4, 4.8, n)
            morphology = rng.beta(2.4, 2.1, n)

        valve_context = rng.choice(["pulmonary_reference", "aortic_high_pressure"], size=n, p=[0.6, 0.4])
        context_target_gradient = np.where(
            valve_context == "pulmonary_reference", CAL["pulmonary_gradient_target"], CAL["aortic_gradient_target"]
        )

        thickness_z = (thickness - 0.62) / 0.12
        half_life_z = (degradation_half_life - CAL["mw_half_life_months"]) / 10.0
        noise_gradient = rng.normal(0, 3.2 if g != "full_data_informed" else 2.0, n)
        predicted_gradient = (
            context_target_gradient
            + 5.5 * thickness_z
            - 4.0 * (pore_index - 0.55)
            + 2.2 * material_risk
            + 1.8 * (valve_context == "aortic_high_pressure")
            + noise_gradient
        )
        predicted_gradient = np.clip(predicted_gradient, 5, 70)

        eoa = CAL["eoa_target"] + 0.35 * (pore_index - 0.55) - 0.22 * thickness_z + rng.normal(0, 0.08, n)
        eoa = np.clip(eoa, 0.6, 3.5)

        residual_mw_12m = 100 * np.exp(-np.log(2) * 12 / degradation_half_life)
        residual_mw_12m = residual_mw_12m - 4.5 * material_risk + rng.normal(0, 2.5, n)
        residual_mw_12m = np.clip(residual_mw_12m, 5, 100)

        integration_latent = (
            -0.2
            + 1.35 * pore_index
            + 0.85 * morphology
            - 0.25 * np.abs(thickness_z)
            - 0.45 * material_risk
            + 0.18 * (degradation_half_life > 40)
            + rng.normal(0, 0.35 if g != "full_data_informed" else 0.25, n)
        )
        integration_12m = sigmoid(integration_latent)
        # Anchor full data-informed group to observed median target without forcing equality
        if g == "full_data_informed":
            integration_12m = np.clip(0.65 * integration_12m + 0.35 * CAL["integration_12m_target"] + rng.normal(0, 0.03, n), 0, 1)

        region_specificity_ratio = np.exp(
            np.log(max(CAL["hinge_belly_ratio_12m"], 0.2))
            + 0.35 * morphology
            - 0.25 * pore_index
            + rng.normal(0, 0.35, n)
        )
        region_specificity_ratio = np.clip(region_specificity_ratio, 0.2, 12)

        regurgitation_fraction = (
            CAL["regurgitation_target"]
            + 6.0 * (thickness < 0.30)
            + 3.0 * (predicted_gradient > context_target_gradient + 12)
            - 2.0 * (eoa > CAL["eoa_target"])
            + rng.normal(0, 2.5, n)
        )
        regurgitation_fraction = np.clip(regurgitation_fraction, 0, 50)

        failure_risk = sigmoid(
            -3.0
            + 2.8 * (thickness < 0.30)
            + 1.8 * (residual_mw_12m < 70)
            + 1.4 * (predicted_gradient > context_target_gradient + 10)
            + 0.9 * material_risk
            + 0.8 * (regurgitation_fraction > 20)
            + rng.normal(0, 0.35, n)
        )

        gradient_penalty = np.abs(predicted_gradient - context_target_gradient) / np.maximum(context_target_gradient, 1)
        eoa_penalty = np.maximum(0, CAL["eoa_target"] - eoa) / CAL["eoa_target"]
        thin_penalty = np.maximum(0, 0.32 - thickness) / 0.32
        score = (
            0.42 * integration_12m
            + 0.18 * (residual_mw_12m / 100)
            + 0.15 * (eoa / 3.0)
            - 0.12 * gradient_penalty
            - 0.16 * failure_risk
            - 0.08 * material_risk
            - 0.10 * thin_penalty
        )
        score = np.clip(score, -1, 1)

        feasible = (
            (thickness >= 0.32)
            & (predicted_gradient <= context_target_gradient + 12)
            & (regurgitation_fraction <= 25)
            & (residual_mw_12m >= 65)
            & (integration_12m >= 0.45)
            & (failure_risk <= 0.35)
        )

        frames.append(pd.DataFrame({
            "scenario_group": g,
            "valve_context": valve_context,
            "thickness_mm": thickness,
            "pore_index_0_1": pore_index,
            "degradation_half_life_months": degradation_half_life,
            "material_risk_index_0_1": material_risk,
            "morphology_index_0_1": morphology,
            "target_gradient_mmHg": context_target_gradient,
            "predicted_gradient_mmHg": predicted_gradient,
            "effective_orifice_area_cm2": eoa,
            "residual_Mw_12m_percent": residual_mw_12m,
            "integration_12m_index_0_1": integration_12m,
            "region_specificity_hinge_belly_ratio": region_specificity_ratio,
            "regurgitation_fraction_percent": regurgitation_fraction,
            "failure_risk_0_1": failure_risk,
            "objective_score": score,
            "feasible": feasible.astype(int),
            "result_type": "scenario_model_generated_from_literature_calibrated_priors",
        }))
    return pd.concat(frames, ignore_index=True)

candidates = generate_candidate_scenarios()
safe_to_csv(candidates, "table_04_generated_candidate_scenarios")


def summarize_candidates(df: pd.DataFrame):
    metrics = [
        "objective_score", "predicted_gradient_mmHg", "effective_orifice_area_cm2",
        "residual_Mw_12m_percent", "integration_12m_index_0_1",
        "regurgitation_fraction_percent", "failure_risk_0_1", "feasible"
    ]
    summary = (
        df.groupby("scenario_group")
        .agg(**{
            f"{m}_mean": (m, "mean") for m in metrics
        })
        .reset_index()
    )
    # Add median and SD in a second pass for cleaner table
    for m in metrics:
        tmp = df.groupby("scenario_group")[m].agg(["median", "std", "count"]).reset_index()
        summary = summary.merge(tmp.rename(columns={"median": f"{m}_median", "std": f"{m}_sd", "count": f"{m}_n"}), on="scenario_group", how="left")
    safe_to_csv(summary, "table_05_candidate_group_summary")

    top = (
        df[df["feasible"] == 1]
        .sort_values("objective_score", ascending=False)
        .head(50)
        .reset_index(drop=True)
    )
    safe_to_csv(top, "table_06_top_50_feasible_designs")
    return summary, top

candidate_summary, top_designs = summarize_candidates(candidates)

# -----------------------------------------------------------------------------
# 9. STATISTICAL TESTS AND COMPARISONS
# -----------------------------------------------------------------------------
def run_statistical_tests(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    endpoints = [
        "objective_score", "predicted_gradient_mmHg", "integration_12m_index_0_1",
        "residual_Mw_12m_percent", "regurgitation_fraction_percent", "failure_risk_0_1", "feasible"
    ]
    groups = list(df["scenario_group"].unique())
    for endpoint in endpoints:
        arrays = [df.loc[df["scenario_group"] == g, endpoint].dropna().values for g in groups]
        if STATS is not None and all(len(a) > 2 for a in arrays):
            try:
                stat, p = STATS.kruskal(*arrays)
                rows.append({
                    "comparison": "global",
                    "test": "Kruskal-Wallis",
                    "endpoint": endpoint,
                    "group_a": "all",
                    "group_b": "all",
                    "statistic": stat,
                    "p_value": p,
                    "effect_size_or_delta": np.nan,
                    "mean_difference_bootstrap_low_high": "",
                })
            except Exception:
                pass
        for a, b in combinations(groups, 2):
            xa = df.loc[df["scenario_group"] == a, endpoint].dropna().values
            xb = df.loc[df["scenario_group"] == b, endpoint].dropna().values
            if len(xa) < 3 or len(xb) < 3:
                continue
            if STATS is not None:
                try:
                    stat, p = STATS.mannwhitneyu(xa, xb, alternative="two-sided", method="auto")
                except Exception:
                    stat, p = np.nan, np.nan
            else:
                stat, p = np.nan, np.nan
            delta = cliffs_delta(xa, xb)
            diff = xa - np.mean(xb)
            mean_diff, low, high = bootstrap_ci(diff, func=np.mean, seed=RANDOM_SEED)
            rows.append({
                "comparison": "pairwise",
                "test": "Mann-Whitney U" if STATS is not None else "effect-size-only",
                "endpoint": endpoint,
                "group_a": a,
                "group_b": b,
                "statistic": stat,
                "p_value": p,
                "effect_size_or_delta": delta,
                "mean_difference_bootstrap_low_high": f"{mean_diff:.4f} [{low:.4f}, {high:.4f}]",
            })
    tests = pd.DataFrame(rows)
    safe_to_csv(tests, "table_07_statistical_tests_and_effect_sizes")
    return tests

stat_tests = run_statistical_tests(candidates)

# -----------------------------------------------------------------------------
# 10. PARETO FRONT AND FEATURE IMPORTANCE
# -----------------------------------------------------------------------------
def pareto_efficient(points: np.ndarray) -> np.ndarray:
    """Return Boolean mask for nondominated points. Assumes all columns are to maximize."""
    n = points.shape[0]
    is_efficient = np.ones(n, dtype=bool)
    for i in range(n):
        if not is_efficient[i]:
            continue
        dominates_i = np.all(points >= points[i], axis=1) & np.any(points > points[i], axis=1)
        dominates_i[i] = False
        if np.any(dominates_i):
            is_efficient[i] = False
    return is_efficient


def build_pareto_table(df: pd.DataFrame):
    d = df.copy()
    # Maximize integration, residual Mw, EOA, inverse gradient deviation, inverse failure risk
    d["gradient_closeness"] = 1 - np.abs(d["predicted_gradient_mmHg"] - d["target_gradient_mmHg"]) / d["target_gradient_mmHg"]
    d["gradient_closeness"] = d["gradient_closeness"].clip(0, 1)
    d["safety_score"] = 1 - d["failure_risk_0_1"]
    # Use a manageable top subset by objective to avoid O(n^2) too large
    sub = d.sort_values("objective_score", ascending=False).head(1800).copy()
    pts = sub[["integration_12m_index_0_1", "residual_Mw_12m_percent", "gradient_closeness", "safety_score"]].copy()
    pts["residual_Mw_12m_percent"] = pts["residual_Mw_12m_percent"] / 100
    mask = pareto_efficient(pts.values)
    sub["pareto_efficient"] = mask.astype(int)
    safe_to_csv(sub, "table_08_pareto_candidate_subset")
    return sub

pareto_df = build_pareto_table(candidates)


def feature_importance_table(df: pd.DataFrame) -> pd.DataFrame:
    features = [
        "thickness_mm", "pore_index_0_1", "degradation_half_life_months",
        "material_risk_index_0_1", "morphology_index_0_1",
        "predicted_gradient_mmHg", "effective_orifice_area_cm2",
        "residual_Mw_12m_percent", "regurgitation_fraction_percent", "failure_risk_0_1",
    ]
    target = "objective_score"
    d = df.dropna(subset=features + [target]).copy()
    rows = []

    # Spearman fallback / baseline
    if STATS is not None:
        for f in features:
            rho, p = STATS.spearmanr(d[f], d[target])
            rows.append({"method": "spearman", "feature": f, "importance": abs(rho), "signed_value": rho, "p_value": p})
    else:
        for f in features:
            rho = np.corrcoef(d[f], d[target])[0, 1]
            rows.append({"method": "pearson_fallback", "feature": f, "importance": abs(rho), "signed_value": rho, "p_value": np.nan})

    # Optional random forest permutation importance. Disabled by default to keep the script fast.
    if RUN_RF_PERMUTATION_IMPORTANCE:
        try:
            from sklearn.ensemble import RandomForestRegressor
            from sklearn.inspection import permutation_importance
            from sklearn.model_selection import train_test_split
            from sklearn.metrics import r2_score
            X = d[features].values
            y = d[target].values
            # Subsample for speed while preserving reproducibility.
            if len(d) > 2500:
                idx = np.random.default_rng(RANDOM_SEED).choice(np.arange(len(d)), size=2500, replace=False)
                X = X[idx]
                y = y[idx]
            X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=RANDOM_SEED)
            rf = RandomForestRegressor(n_estimators=80, random_state=RANDOM_SEED, min_samples_leaf=4, n_jobs=1)
            rf.fit(X_train, y_train)
            pred = rf.predict(X_test)
            r2 = r2_score(y_test, pred)
            perm = permutation_importance(rf, X_test, y_test, n_repeats=5, random_state=RANDOM_SEED, n_jobs=1)
            for f, imp, sd in zip(features, perm.importances_mean, perm.importances_std):
                rows.append({"method": f"rf_permutation_R2_{r2:.3f}", "feature": f, "importance": imp, "signed_value": imp, "p_value": np.nan})
        except Exception as e:
            rows.append({"method": "rf_permutation_not_run", "feature": str(e), "importance": np.nan, "signed_value": np.nan, "p_value": np.nan})
    else:
        rows.append({"method": "rf_permutation_skipped_fast_mode", "feature": "set RUN_RF_PERMUTATION_IMPORTANCE=True to run", "importance": np.nan, "signed_value": np.nan, "p_value": np.nan})

    imp = pd.DataFrame(rows).sort_values(["method", "importance"], ascending=[True, False])
    safe_to_csv(imp, "table_09_feature_importance")
    return imp

feature_importance = feature_importance_table(candidates)

# -----------------------------------------------------------------------------
# 11. SCENARIO FIGURES
# -----------------------------------------------------------------------------
def plot_candidate_scenarios(df: pd.DataFrame, summary: pd.DataFrame, pareto: pd.DataFrame, importance: pd.DataFrame):
    order = ["synthetic_only", "hydrodynamic_only", "biology_only", "full_data_informed"]
    colors = [COL["gray"], COL["burgundy"], COL["forest"], COL["navy"]]
    color_map = dict(zip(order, colors))

    # Objective distributions boxplot
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    data = [df.loc[df["scenario_group"] == g, "objective_score"].values for g in order]
    bp = ax.boxplot(data, labels=[g.replace("_", "\n") for g in order], patch_artist=True, showfliers=False)
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_alpha(0.75)
        patch.set_edgecolor("black")
    for med in bp["medians"]:
        med.set_color("black")
        med.set_linewidth(1.2)
    ax.set_ylabel("Objective score")
    ax.grid(True, axis="y", linestyle=":")
    subplot_label(ax, "(a)")
    save_figure(fig, "fig_12_scenario_objective_score_comparison")

    # Feasibility rates
    feasible = df.groupby("scenario_group")["feasible"].mean().reindex(order)
    fig, ax = plt.subplots(figsize=(6.8, 4.0))
    ax.bar(np.arange(len(order)), feasible.values, color=colors, edgecolor="black", linewidth=0.6, alpha=0.88)
    ax.set_xticks(np.arange(len(order)))
    ax.set_xticklabels([g.replace("_", "\n") for g in order])
    ax.set_ylim(0, 1)
    ax.set_ylabel("Feasibility rate")
    ax.grid(True, axis="y", linestyle=":")
    subplot_label(ax, "(b)")
    save_figure(fig, "fig_13_scenario_feasibility_rates")

    # Design space: thickness vs gradient, colored by group; use sampled subset
    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    rng = np.random.default_rng(RANDOM_SEED)
    for g in order:
        dg = df[df["scenario_group"] == g]
        if len(dg) > 700:
            dg = dg.sample(700, random_state=RANDOM_SEED)
        ax.scatter(dg["thickness_mm"], dg["predicted_gradient_mmHg"], s=9, alpha=0.35, color=color_map[g], label=g.replace("_", " "))
    ax.axvline(0.32, color="black", linestyle=":", linewidth=1)
    ax.set_xlabel("Thickness (mm)")
    ax.set_ylabel("Predicted gradient (mmHg)")
    ax.grid(True, linestyle=":")
    ax.legend(loc="upper left", fontsize=8)
    subplot_label(ax, "(c)")
    save_figure(fig, "fig_14_design_space_thickness_vs_gradient")

    # Pareto front
    if not pareto.empty:
        fig, ax = plt.subplots(figsize=(6.5, 4.8))
        p0 = pareto[pareto["pareto_efficient"] == 0]
        p1 = pareto[pareto["pareto_efficient"] == 1]
        ax.scatter(p0["predicted_gradient_mmHg"], p0["integration_12m_index_0_1"], s=10, alpha=0.25, color=COL["gray"], label="Dominated")
        ax.scatter(p1["predicted_gradient_mmHg"], p1["integration_12m_index_0_1"], s=18, alpha=0.85, color=COL["burgundy"], label="Pareto-efficient")
        ax.set_xlabel("Predicted gradient (mmHg)")
        ax.set_ylabel("Integration index at 12 months")
        ax.set_ylim(0, 1.05)
        ax.grid(True, linestyle=":")
        ax.legend(loc="best")
        subplot_label(ax, "(d)")
        save_figure(fig, "fig_15_pareto_front_gradient_integration")

    # Endpoint ablation panel
    endpoints = ["integration_12m_index_0_1", "residual_Mw_12m_percent", "failure_risk_0_1", "regurgitation_fraction_percent"]
    labels = ["Integration", "Residual Mw (%)", "Failure risk", "Regurgitation (%)"]
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    x = np.arange(len(order))
    width = 0.18
    for i, endpoint in enumerate(endpoints):
        vals = df.groupby("scenario_group")[endpoint].mean().reindex(order).values
        # scale mixed-unit endpoints for a single comparison display
        if endpoint in ["residual_Mw_12m_percent", "regurgitation_fraction_percent"]:
            vals = vals / 100
        ax.bar(x + (i - 1.5) * width, vals, width=width, label=labels[i], color=[COL["navy"], COL["burgundy"], COL["forest"], COL["plum"]][i], alpha=0.84, edgecolor="black", linewidth=0.4)
    ax.set_xticks(x)
    ax.set_xticklabels([g.replace("_", "\n") for g in order])
    ax.set_ylabel("Scaled mean endpoint")
    ax.grid(True, axis="y", linestyle=":")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.18), ncol=2, fontsize=8)
    subplot_label(ax, "(e)")
    save_figure(fig, "fig_16_ablation_scaled_endpoint_comparison")

    # Correlation heatmap using imshow, no seaborn
    feat_cols = [
        "thickness_mm", "pore_index_0_1", "degradation_half_life_months", "material_risk_index_0_1",
        "morphology_index_0_1", "predicted_gradient_mmHg", "effective_orifice_area_cm2",
        "residual_Mw_12m_percent", "integration_12m_index_0_1", "failure_risk_0_1", "objective_score"
    ]
    corr = df[feat_cols].corr(method="spearman")
    fig, ax = plt.subplots(figsize=(7.2, 6.2))
    im = ax.imshow(corr.values, vmin=-1, vmax=1, cmap="Greys" if STYLE_MODE == "grayscale" else "bone_r")
    ax.set_xticks(np.arange(len(feat_cols)))
    ax.set_yticks(np.arange(len(feat_cols)))
    ax.set_xticklabels([c.replace("_", "\n") for c in feat_cols], rotation=90, fontsize=7)
    ax.set_yticklabels([c.replace("_", " ") for c in feat_cols], fontsize=7)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Spearman r")
    subplot_label(ax, "(f)")
    save_figure(fig, "fig_17_spearman_correlation_heatmap")

    # Feature importance
    if not importance.empty:
        imp = importance[importance["method"].astype(str).str.startswith("spearman")].copy()
        if imp.empty:
            imp = importance.dropna(subset=["importance"]).copy()
        imp = imp.sort_values("importance", ascending=True).tail(12)
        fig, ax = plt.subplots(figsize=(6.8, 4.8))
        ax.barh(imp["feature"].astype(str), imp["importance"], color=COL["navy"], alpha=0.88, edgecolor="black", linewidth=0.4)
        ax.set_xlabel("Importance magnitude")
        ax.set_ylabel("Feature")
        ax.grid(True, axis="x", linestyle=":")
        subplot_label(ax, "(g)")
        save_figure(fig, "fig_18_feature_importance_objective_score")

plot_candidate_scenarios(candidates, candidate_summary, pareto_df, feature_importance)

# -----------------------------------------------------------------------------
# 12. CALIBRATION-VS-GENERATED CHECKS
# -----------------------------------------------------------------------------
def plot_calibration_generated_checks(df: pd.DataFrame):
    full = df[df["scenario_group"] == "full_data_informed"].copy()
    if full.empty:
        return
    # Gradient target check by context
    contexts = full["valve_context"].unique()
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    data = [full.loc[full["valve_context"] == c, "predicted_gradient_mmHg"].values for c in contexts]
    bp = ax.boxplot(data, labels=[str(c).replace("_", "\n") for c in contexts], patch_artist=True, showfliers=False)
    for patch, c in zip(bp["boxes"], [COL["navy"], COL["burgundy"], COL["forest"]]):
        patch.set_facecolor(c)
        patch.set_alpha(0.75)
        patch.set_edgecolor("black")
    for c in contexts:
        target = full.loc[full["valve_context"] == c, "target_gradient_mmHg"].median()
        ax.axhline(target, color="black", linestyle=":", linewidth=0.9)
    ax.set_ylabel("Predicted gradient (mmHg)")
    ax.grid(True, axis="y", linestyle=":")
    subplot_label(ax, "(a)")
    save_figure(fig, "fig_19_calibration_check_gradient_by_context")

    # Generated vs target endpoint dashboard
    values = {
        "Integration\nmedian": full["integration_12m_index_0_1"].median(),
        "Residual Mw\nmedian/100": full["residual_Mw_12m_percent"].median() / 100,
        "Failure risk\nmedian": full["failure_risk_0_1"].median(),
        "Feasible\nrate": full["feasible"].mean(),
    }
    fig, ax = plt.subplots(figsize=(6.5, 4.0))
    ax.bar(range(len(values)), list(values.values()), color=[COL["navy"], COL["forest"], COL["burgundy"], COL["plum"]], alpha=0.88, edgecolor="black", linewidth=0.5)
    ax.set_xticks(range(len(values)))
    ax.set_xticklabels(list(values.keys()))
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Scaled value")
    ax.grid(True, axis="y", linestyle=":")
    subplot_label(ax, "(b)")
    save_figure(fig, "fig_20_full_data_informed_endpoint_dashboard")

plot_calibration_generated_checks(candidates)

# -----------------------------------------------------------------------------
# 13. SOURCE-LEVEL AND EXTRACTED METRIC VISUALIZATION
# -----------------------------------------------------------------------------
def plot_extracted_metric_distributions(extracted_df: pd.DataFrame):
    if extracted_df.empty or not {"metric_group", "mean"}.issubset(extracted_df.columns):
        return
    d = extracted_df.dropna(subset=["metric_group", "mean"]).copy()
    groups = d["metric_group"].value_counts().head(8).index.tolist()
    d = d[d["metric_group"].isin(groups)]
    if d.empty:
        return
    # Strip plot-like display without seaborn
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    rng = np.random.default_rng(RANDOM_SEED)
    for i, g in enumerate(groups):
        vals = d.loc[d["metric_group"] == g, "mean"].values.astype(float)
        x = rng.normal(i, 0.035, size=len(vals))
        ax.scatter(x, vals, s=18, alpha=0.65, color=COL[["navy", "burgundy", "forest", "plum", "brown", "slate", "gray", "charcoal"][i % 8]])
        if len(vals) > 0:
            ax.plot([i - 0.18, i + 0.18], [np.nanmedian(vals), np.nanmedian(vals)], color="black", linewidth=1.1)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([g.replace("_", "\n") for g in groups], rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Extracted numeric value (native units)")
    ax.grid(True, axis="y", linestyle=":")
    subplot_label(ax, "(a)")
    save_figure(fig, "fig_21_extracted_metric_value_distributions")

plot_extracted_metric_distributions(extracted)

# -----------------------------------------------------------------------------
# 14. MANUSCRIPT-READY OUTPUT EXCEL
# -----------------------------------------------------------------------------
def manuscript_claim_flags() -> pd.DataFrame:
    rows = [
        {
            "old_or_risky_claim": "The optimizer reduces scaffold persistence by ~79% vs Yacoub.",
            "data_informed_revision": "Yacoub should be used as a reference/calibration case, not as a superiority comparator.",
            "basis": "Yacoub molecular-weight retention and hydrodynamic values indicate persistent scaffold response over months.",
            "claim_strength": "moderate; literature-calibrated, not clinically validated",
        },
        {
            "old_or_risky_claim": "Colonization is uniformly high across the scaffold.",
            "data_informed_revision": "Biological integration is region-specific; hinge/belly/tip patterns should be reported separately.",
            "basis": "Region-specific cell density and hinge/belly ratios from extracted data.",
            "claim_strength": "stronger than uniform claim; still source-limited",
        },
        {
            "old_or_risky_claim": "Pressure gradient is an arbitrary surrogate.",
            "data_informed_revision": "Pressure-gradient endpoint is calibrated to extracted mmHg hydrodynamic targets.",
            "basis": "Yacoub HCCV and Vis aortic sheep gradient references.",
            "claim_strength": "improved; pulmonary/aortic contexts must not be pooled uncritically",
        },
        {
            "old_or_risky_claim": "RL gives a patient-specific scaffold prescription.",
            "data_informed_revision": "The framework performs real-data-informed exploratory design optimization under uncertainty.",
            "basis": "No patient-level clinical validation; data are literature-extracted and aggregate.",
            "claim_strength": "safe and defensible",
        },
    ]
    df = pd.DataFrame(rows)
    safe_to_csv(df, "table_10_claim_revision_flags")
    return df

claim_flags = manuscript_claim_flags()


def figure_index() -> pd.DataFrame:
    figs = sorted(OUT["figures"].glob("*.png"))
    rows = []
    for i, f in enumerate(figs, start=1):
        rows.append({
            "figure_no": i,
            "file_name": f.name,
            "suggested_use": "main_text" if i <= 12 else "supplementary",
            "notes": "PNG only; 600 DPI; academic dark style" if STYLE_MODE == "dark_academic" else "PNG only; strict grayscale style",
        })
    df = pd.DataFrame(rows)
    safe_to_csv(df, "table_11_figure_index")
    return df

fig_index = figure_index()

# Save final Excel workbook with core outputs
final_xlsx = OUT["excel"] / "TEHV_all_possible_generated_results.xlsx"
with pd.ExcelWriter(final_xlsx, engine="openpyxl") as writer:
    if not extracted.empty:
        extracted.to_excel(writer, sheet_name="Clean_Extracted_Numeric", index=False)
    if not source_log.empty:
        source_log.to_excel(writer, sheet_name="Source_Log", index=False)
    if not scaffold_response.empty:
        scaffold_response.to_excel(writer, sheet_name="Scaffold_Response", index=False)
    if not parameter_map.empty:
        parameter_map.to_excel(writer, sheet_name="Parameter_Map", index=False)
    candidate_summary.to_excel(writer, sheet_name="Scenario_Group_Summary", index=False)
    top_designs.to_excel(writer, sheet_name="Top_50_Feasible_Designs", index=False)
    stat_tests.to_excel(writer, sheet_name="Statistical_Tests", index=False)
    pareto_df.to_excel(writer, sheet_name="Pareto_Subset", index=False)
    feature_importance.to_excel(writer, sheet_name="Feature_Importance", index=False)
    claim_flags.to_excel(writer, sheet_name="Claim_Revision_Flags", index=False)
    fig_index.to_excel(writer, sheet_name="Figure_Index", index=False)
    pd.DataFrame([CAL]).to_excel(writer, sheet_name="Calibration_Constants", index=False)

# Copy input files into output log folder for traceability, if present
for p in [DATASET_PATH, RESULTS_PATH]:
    if p is not None and p.exists():
        try:
            shutil.copy2(p, OUT["logs"] / p.name)
        except Exception:
            pass

# Run earlier plotting functions
try:
    audit, metric_summary, source_summary = build_data_audit_tables()
    plot_data_audit(audit, metric_summary)
except Exception as e:
    print(f"Data audit plotting skipped: {e}")

try:
    plot_hydrodynamic_targets(hydro)
    plot_regional_colonization(regional)
    plot_integration_timecourse(integration, mc_summary)
    plot_polymer_degradation(polymer, trend_fits)
    plot_mc_distributions(mc_summary)
except Exception as e:
    print(f"Real-data figures skipped or partially completed: {e}")

# Regenerate figure index after all plots
fig_index = figure_index()
with pd.ExcelWriter(final_xlsx, engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
    fig_index.to_excel(writer, sheet_name="Figure_Index", index=False)

# Write final run log
write_log({
    "timestamp": datetime.now().isoformat(),
    "output_dir": str(OUTPUT_DIR),
    "dataset_path": str(DATASET_PATH) if DATASET_PATH else None,
    "results_path": str(RESULTS_PATH) if RESULTS_PATH else None,
    "style_mode": STYLE_MODE,
    "save_png_only": SAVE_PNG_ONLY,
    "n_candidate_rows": int(len(candidates)),
    "n_figures": len(list(OUT["figures"].glob("*.png"))),
    "final_excel": str(final_xlsx),
    "scientific_warning": "Scenario results are literature-calibrated and exploratory; not clinical validation.",
})

print("\nDone.")
print(f"Output directory: {OUTPUT_DIR}")
print(f"Final Excel: {final_xlsx}")
print(f"Figures: {OUT['figures']}")
print(f"Tables: {OUT['tables']}")
print(f"PNG figures generated: {len(list(OUT['figures'].glob('*.png')))}")
print("No PDF files were saved.")
