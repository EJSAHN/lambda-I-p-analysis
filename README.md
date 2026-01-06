# λ–I–p Framework (Pest Management Science revision package)

Reproducible analysis code for the **λ–I–p** morphology-based geometric + probabilistic framework used for mechanistic antifungal screening and mode-of-action (MoA) classification.

This repository generates a single, self-contained workbook (**`SupplementaryData_S1.xlsx`**) from raw colony-morphology spreadsheets and then produces all downstream tables/figures **exclusively from that workbook** (numerical reproducibility).

## What this repo does

Given three input spreadsheets (not included here):

1) **Chemical single-dose screen** (colony morphology; includes *Area* and *LWR*), and  
2) **UV-C cacao dataset** (time gradient; includes *Area*, *LWR*, exposure time), and  
3) **UV-C coffee dataset** (time gradient; includes *Area*, *LWR*, exposure time),

the pipeline computes:

- **Potency**:  
  \[
  \lambda_{\mathrm{area}}=-\ln(S_{\mathrm{area}}),\quad S_{\mathrm{area}}=\frac{\mathrm{Area}_{t}}{\mathrm{Area}_{c}}
  \]

- **Polarity**:  
  \[
  I_{\mathrm{LWR}}=\ln\left(\frac{\mathrm{LWR}_{t}}{\mathrm{LWR}_{c}}\right)
  \]

- **Ellipse axis reconstruction** (major/minor axis ratios):  
  \[
  a_{\mathrm{ratio}}=\sqrt{S_{\mathrm{area}}\cdot e^{I}},\quad b_{\mathrm{ratio}}=\sqrt{S_{\mathrm{area}}/e^{I}}
  \]

- **Event probability** (chemical): event if replicate **LWR ≥ τ**, with exact **Clopper–Pearson 95% CI**  
- **Event probability** (UV-C): event if **I_rep ≥ ln(τ)** (log-ratio space), with exact CI

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

### 2) Put input spreadsheets under `./data/`

Default filenames expected:

- `data/Supplementary Data1_ISME Com.xlsx`  (chemical screen)
- `data/Cacao.xlsx`                         (UV-C cacao)
- `data/Coffee.xlsx`                        (UV-C coffee)

Or you can pass explicit paths via CLI flags.

### 3) Generate SupplementaryData_S1.xlsx (single command)

```bash
python scripts/chem_analysis_main.py \
  --chem-xlsx  "data/Supplementary Data1_ISME Com.xlsx" \
  --cacao-xlsx data/Cacao.xlsx \
  --coffee-xlsx data/Coffee.xlsx \
  --out-xlsx   outputs/SupplementaryData_S1.xlsx
```

You can also use `--outdir outputs` instead of `--out-xlsx ...`.

## Input data sources (not hosted here)

This repository does **not** redistribute third-party datasets. Please obtain them from the original sources:

- **Cacao UV-C reference dataset**: Baek et al. (2025), *Scientific Reports* (supplementary material). DOI: `10.1038/s41598-025-20277-2`
- **Coffee UV-C validation dataset**: Zenodo record `10.5281/zenodo.18133112` (file: `Coffee.xlsx`)
- Coffee UV-C dataset is described in: *Food Control* (2026) DOI: `10.1016/j.foodcont.2026.111956`
- Chemical cacao dataset related publication: *Journal of Fungi* (2026) DOI: `10.3390/jof12010033`

## Reproducibility notes

- The analysis is intentionally **manuscript-locked**: defaults (e.g., τ=1.10 primary; τ∈{1.10,1.15,1.20} sensitivity) match the paper text.
- The generated workbook includes a `metadata` sheet with file hashes (SHA-256), package versions, and key parameters.
- If SciPy is installed, **exact Clopper–Pearson** intervals are used; otherwise the script falls back to Wilson intervals and records this in metadata.

## Repository layout

- `scripts/chem_analysis_main.py` — main pipeline; writes `SupplementaryData_S1.xlsx`
- `scripts/make_figures.py` — UV-C figures from `SupplementaryData_S1.xlsx`
- `scripts/make_table_S3.py` — generates the “near-equal potency, divergent polarity” examples table
- `scripts/make_graphical_abstract_from_S1.py` — optional GA generator

## License

See `LICENSE`.
