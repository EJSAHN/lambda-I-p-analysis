# λ–I–p Framework

Analysis code for the λ–I–p colony-morphology framework for mechanistic antifungal screening.

This repository generates a processed summary workbook (`SupplementaryData_S1.xlsx`) from locally supplied input spreadsheets and derives downstream tables and figures from that workbook.

## Overview

The pipeline uses three external input spreadsheets:

- a chemical fixed-dose screen with colony morphology measurements (including `Area` and `LWR`)
- a UV-C cacao dataset with exposure-time gradients
- a UV-C coffee dataset with exposure-time gradients

From these inputs, the code computes:

- **Potency**
  \[
  \lambda_{\mathrm{area}}=-\ln(S_{\mathrm{area}}), \qquad
  S_{\mathrm{area}}=\frac{\mathrm{Area}_{t}}{\mathrm{Area}_{c}}
  \]

- **Polarity**
  \[
  I_{\mathrm{LWR}}=\ln\left(\frac{\mathrm{LWR}_{t}}{\mathrm{LWR}_{c}}\right)
  \]

- **Ellipse-based axis reconstruction**
  \[
  a_{\mathrm{ratio}}=\sqrt{S_{\mathrm{area}}\cdot e^{I}}, \qquad
  b_{\mathrm{ratio}}=\sqrt{S_{\mathrm{area}}/e^{I}}
  \]

- **Event probability (chemical screen)**  
  event defined as replicate `LWR \geq \tau`, with exact Clopper–Pearson 95% confidence intervals when SciPy is available

- **Event probability (UV-C)**  
  event defined as `I_{rep} \geq \ln(\tau)`, with exact Clopper–Pearson 95% confidence intervals when SciPy is available

## What this repository includes

- analysis code
- dependency list
- instructions for generating the summary workbook

## What this repository does not include

This repository does **not** redistribute the raw input spreadsheets used in the analyses.

The output workbook `SupplementaryData_S1.xlsx` is a **processed summary workbook** containing the numerical tables used for figures and supplementary analyses. It is **not** a redistribution of the raw colony-morphology datasets.

## Quick start

### 1) Create a Python environment

```bash
python -m venv .venv
# Windows:
.venv\Scripts\activate
# macOS/Linux:
source .venv/bin/activate

pip install -r requirements.txt
```

### 2) Prepare the input spreadsheets locally

Store the required input spreadsheets in a local directory of your choice.

This repository does not require any specific input filenames if you provide explicit command-line arguments.

### 3) Generate the summary workbook

```bash
python scripts/chem_analysis_main.py \
  --chem-xlsx  "<path_to_chemical_input.xlsx>" \
  --cacao-xlsx "<path_to_cacao_uv_input.xlsx>" \
  --coffee-xlsx "<path_to_coffee_uv_input.xlsx>" \
  --out-xlsx   "outputs/SupplementaryData_S1.xlsx"
```

You may also use:

```bash
python scripts/chem_analysis_main.py \
  --chem-xlsx  "<path_to_chemical_input.xlsx>" \
  --cacao-xlsx "<path_to_cacao_uv_input.xlsx>" \
  --coffee-xlsx "<path_to_coffee_uv_input.xlsx>" \
  --outdir outputs
```

## Input data sources

Input datasets are not hosted in this repository. They should be obtained from the original published sources listed below.

- **Chemical screening source data**: available through the supplementary Excel file associated with the published *Smart Agricultural Technology* paper. DOI: `10.1016/j.atech.2026.101895`
- **UV-C coffee validation data**: available from Zenodo as `Coffee.xlsx`. DOI: `10.5281/zenodo.18133112`
- **UV-C cacao reference data**: available as supplementary material in Baek et al. (2025), *Scientific Reports*. DOI: `10.1038/s41598-025-20277-2`
- **UV-C dataset description and code provenance**: documented in the *Food Control* paper. DOI: `10.1016/j.foodcont.2026.111956`

## Reproducibility notes

- The analysis code computes the summary workbook from external input spreadsheets supplied by the user.
- The summary workbook contains processed tables used for manuscript figures and supplementary analyses.
- Exact Clopper–Pearson intervals are used when SciPy is available.

## Repository layout

- `scripts/chem_analysis_main.py` — main analysis pipeline
- `requirements.txt` — Python dependencies
- `LICENSE` — repository license information

## Citation and data use

Please cite the associated publications when using the code or source datasets.

Users are responsible for obtaining and using the external input datasets in accordance with the terms of the original data providers.
