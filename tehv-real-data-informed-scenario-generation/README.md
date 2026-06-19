# Real-Data-Informed TEHV Scaffold Scenario Generation

This repository contains the curated dataset, reproducible Python code, and generated analysis outputs for a real-data-informed computational framework for exploratory bioresorbable tissue-engineered heart valve (TEHV) scaffold design.

**Manuscript text is intentionally not included in this repository.** The repository is limited to data, code, tables, and figures needed for transparency and reproducibility.

## Repository structure

```text
tehv-real-data-informed-scenario-generation/
├── code/
│   ├── tehv_generate_all_results.py
│   └── tehv_generate_all_results_converted.ipynb
├── data/
│   ├── TEHV_real_data_extraction_dataset_v2.xlsx
│   └── TEHV_real_data_informed_generated_results_v1.xlsx
├── results/
│   ├── excel_outputs/
│   ├── figures_png/
│   └── tables_csv/
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
