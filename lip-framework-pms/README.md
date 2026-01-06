# λ–I–p framework (reproducible pipeline)

This repository provides a reproducible Python pipeline for the **λ–I–p** morphology framework:

- **λ (potency)** = −ln(S_area), where S_area = (treated mean area) / (control mean area)
- **I (polarity)** = ln(LWR_t / LWR_c)
- **p̂ (event probability)** with **exact 95% Clopper–Pearson** confidence bounds
- **Axis reconstruction** (ellipse model):  
  a_ratio = sqrt(S_area · exp(I)); b_ratio = sqrt(S_area / exp(I))

## What this repo generates
- `outputs/SupplementaryData_S1.xlsx` (multi-sheet workbook used for the manuscript)
- Supplementary Table S1 (CSV + DOCX)

## Install
```bash
python -m venv .venv
# Windows:
.venv\Scripts\activate
# macOS/Linux:
# source .venv/bin/activate
pip install -r requirements.txt
```

## Run
Place input files locally (not redistributed by default):
- Chemical screen Excel (e.g., `data/chemical.xlsx`)
- UV‑C cacao Excel (e.g., `data/Cacao.xlsx`)
- UV‑C coffee Excel (e.g., `data/Coffee.xlsx`)

Then:
```bash
python scripts/chem_analysis_main.py --chem-xlsx data/chemical.xlsx --cacao-xlsx data/Cacao.xlsx --coffee-xlsx data/Coffee.xlsx --outdir outputs
python scripts/make_figures.py --s1 outputs/SupplementaryData_S1.xlsx --outdir outputs
python scripts/make_table_S1.py --s1 outputs/SupplementaryData_S1.xlsx --outdir outputs
```

## Data availability / provenance (not redistributed here)
- Cacao UV‑C dataset: Baek et al. (2025) *Scientific Reports* (DOI: 10.1038/s41598-025-20277-2). **Action:** Download Supplementary Data and rename to **`Cacao.xlsx`**.
- Coffee UV‑C dataset: Ahn et al. (2026) Food Control* (DOI: 10.1016/j.foodcont.2026.111956) (Data available via Zenodo: DOI [10.5281/zenodo.18133371](https://doi.org/10.5281/zenodo.18133371)). 
  **Action:** Download and rename to **`Coffee.xlsx`**.
- Chemical morphology dataset: compiled screening dataset used in the manuscript; ensure you have redistribution rights before uploading any raw files. 
  **Action:** Download the supplementary Excel file and rename to **`Chem.xlsx`**.

## License
MIT (see `LICENSE`).
