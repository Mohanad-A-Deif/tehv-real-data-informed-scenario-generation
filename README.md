# Literature-Calibrated TEHV Scaffold Simulation and Sensitivity Analysis

This repository contains the curated dataset, reproducible Python code, and generated analysis outputs for a literature-calibrated computational framework for simulation and sensitivity analysis of bioresorbable tissue-engineered heart valve (TEHV) scaffolds.

**Manuscript text is intentionally not included in this repository.** The repository is limited to data, code, tables, and figures needed for transparency and reproducibility.

## Repository structure

```text
tehv-real-data-informed-scenario-generation/
├── code/
│   ├── tehv_generate_all_results.py
│   ├── tehv_generate_all_results_converted.ipynb
│   └── robustness/
│       ├── tehv_core.py
│       ├── run_r2_analyses.py
│       ├── r2_extra.py
│       └── make_fig22.py
├── data/
│   ├── TEHV_real_data_extraction_dataset_v2.xlsx
│   └── TEHV_real_data_informed_generated_results_v1.xlsx
├── results/
│   ├── excel_outputs/
│   ├── figures_png/
│   ├── tables_csv/
│   └── robustness_tables/
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

- `data/TEHV_real_data_extraction_dataset_v2.xlsx`: curated TEHV evidence-extraction workbook.
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

## Scenario groups

The four scenario groups are named by the structure of their design priors. The code keeps short internal identifiers:

| Name used in the manuscript | Identifier in code and tables |
|---|---|
| Baseline | `baseline` |
| Intermediate-prior | `hydrodynamic_only` |
| Evidence-centered wide-prior | `biology_only` |
| Evidence-centered narrow-prior | `full_data_informed` |

The loading context enters the gradient model only through the context-specific calibration target; no additional aortic offset is applied.

## Robustness, sensitivity, and held-out analyses

`code/robustness/` contains a parameterised re-implementation of the generator (`tehv_core.py`), which reproduces the main script at default settings, and the analysis scripts. To run:

```bash
cd code/robustness
python run_r2_analyses.py
python make_fig22.py
```

Outputs are written to `code/robustness/tables/`. The stored copies are in `results/robustness_tables/`:

- seed replication (200 seeds);
- one-at-a-time and joint perturbation of the 23 endpoint coefficients, with a classification of the influence of each coefficient;
- one-at-a-time sign reversal of the coefficients;
- one-at-a-time and joint perturbation of the indicator cut-points and the thickness scale;
- alternative and random objective weights;
- feasibility-threshold perturbation;
- stepwise prior ablation and ablation of the anchoring step and group-specific noise levels;
- leave-one-time-point-out and cross-study held-out checks of the calibration targets.

`results/figures_png/fig_22_robustness_summary.png` summarises the coefficient and weight analyses.

## Scientific scope

This repository supports an exploratory literature-calibrated simulation study. The generated scaffold candidates are virtual design scenarios conditioned on literature-derived priors and extracted quantitative endpoints. They are **not** patient-level observations, clinical validation results, or patient-specific treatment recommendations.

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
