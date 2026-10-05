# Real-Data-Informed TEHV Scaffold Scenario Generation

This repository contains the curated dataset, reproducible Python code, and generated analysis outputs for a real-data-informed computational framework for exploratory bioresorbable tissue-engineered heart valve (TEHV) scaffold design.

**Manuscript text is intentionally not included in this repository.** The repository is limited to data, code, tables, and figures needed for transparency and reproducibility.

## Repository structure

```text
tehv-real-data-informed-scenario-generation/
├── code/
│   ├── tehv_generate_all_results.py
│   ├── tehv_generate_all_results_converted.ipynb
│   └── r1_robustness/
│       ├── tehv_core.py
│       └── run_r1_analyses.py
├── data/
│   ├── TEHV_real_data_extraction_dataset_v2.xlsx
│   └── TEHV_real_data_informed_generated_results_v1.xlsx
├── results/
│   ├── excel_outputs/
│   ├── figures_png/
│   ├── tables_csv/
│   └── r1_robustness_tables/
├── logs/
│   └── run_log.json
├── archive/
│   └── all_results_csv.zip
├── README.md
├── requirements.txt
├── .gitignore
└── LICENSE
```

## Contents

- `data/TEHV_real_data_extraction_dataset_v2.xlsx`: curated real-data-informed TEHV extraction workbook.
- `data/TEHV_real_data_informed_generated_results_v1.xlsx`: generated intermediate result workbook used by the result generator.
- `code/tehv_generate_all_results.py`: main reproducible script for generating tables and figures.
- `code/tehv_generate_all_results_converted.ipynb`: notebook version of the same workflow.
- `results/tables_csv/`: exported CSV tables.
- `results/figures_png/`: manuscript-ready PNG figures.
- `results/excel_outputs/`: combined Excel output workbook.
- `archive/all_results_csv.zip`: optional compressed result archive.

## Installation

Create a Python environment and install dependencies:

```bash
pip install -r requirements.txt
```

## Reproducing the results

From the repository root, run:

```bash
python code/tehv_generate_all_results.py
```

By default, regenerated outputs are written to:

```text
generated_outputs/
```

To choose another output folder, set `TEHV_OUTPUT_DIR` before running.

Windows PowerShell:

```powershell
$env:TEHV_OUTPUT_DIR="C:\path\to\output"
python code/tehv_generate_all_results.py
```

Linux/macOS:

```bash
export TEHV_OUTPUT_DIR=/path/to/output
python code/tehv_generate_all_results.py
```

## Revision R1: correction and added analyses

### Correction

The first public version of `code/tehv_generate_all_results.py` contained an error in `derive_calibration_constants()`. The pulmonary pressure-gradient target was assigned by matching text labels in the `Hydrodynamic_Targets` sheet, and the matching rule also captured the row `Vis/Yacoub_gradient_ratio` (value 1.40). The pulmonary target was therefore set to 1.40 instead of 23.37 mmHg, and pulmonary-context candidates were generated around the wrong target.

The rule now excludes the ratio row. The following outputs were regenerated:

- `results/tables_csv/table_04` to `table_09`
- `results/figures_png/fig_12` to `fig_20`
- `results/excel_outputs/TEHV_all_possible_generated_results.xlsx`
- `archive/all_results_csv.zip`

Effect of the correction: feasibility rates and the integration, residual molecular-weight, regurgitation, and failure-risk endpoints are unchanged. The objective score, the predicted gradient, the composition of the top-50 and Pareto sets, and the feature associations changed. For example, the mean objective score of the baseline and full-data-informed groups changed from 0.042 and 0.297 to 0.305 and 0.513.

Two further changes do not affect the numbers: the first scenario group is labelled `baseline` in the code, as in the result tables, and the box-plot calls were made compatible with current matplotlib releases.

### Added analyses

`code/r1_robustness/` contains a parameterised re-implementation of the generator (`tehv_core.py`), which reproduces the main script exactly at default settings, and the analysis script (`run_r1_analyses.py`). To run:

```bash
cd code/r1_robustness
python run_r1_analyses.py
```

Outputs are written to `code/r1_robustness/tables/`. The stored copies are in `results/r1_robustness_tables/`:

- seed replication (200 seeds);
- one-at-a-time and joint perturbation of the 24 endpoint coefficients;
- alternative and random objective weights;
- feasibility-threshold perturbation;
- ablation of the anchoring step and the group-specific noise levels;
- leave-one-time-point-out and cross-study held-out checks of the calibration targets.

`results/figures_png/fig_22_robustness_summary.png` summarises the coefficient and weight analyses.

## Scientific scope

This repository supports an exploratory real-data-informed computational/scenario-generation study. The generated scaffold candidates are virtual design scenarios conditioned on literature-derived priors and extracted quantitative endpoints. They are **not** patient-level observations, clinical validation results, or patient-specific treatment recommendations.

## Main generated outputs

The workflow exports:

- dataset audit tables;
- extracted numerical metric summaries;
- scenario-level generated candidate tables;
- statistical tests and effect-size summaries;
- top feasible candidate summaries;
- Pareto candidate subsets;
- feature-association summaries;
- PNG figures for dataset audit, biological remodeling, polymer persistence, hydrodynamic calibration, scenario comparisons, Pareto screening, and feature associations.

## Citation and data use

If using this repository, cite the associated manuscript once available and cite the original data sources referenced in the manuscript. The curated workbook contains extracted and transformed values from public literature and public/controlled-access resources; users are responsible for respecting the terms of the original sources.
