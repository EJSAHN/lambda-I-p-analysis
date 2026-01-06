# -*- coding: utf-8 -*-
"""
chem_analysis_main.py

Reproduces the manuscript-locked λ–I–p analysis from raw morphology spreadsheets and
writes a single, self-contained Excel workbook (SupplementaryData_S1.xlsx) that downstream
plotting scripts read *exclusively*.

Key outputs (sheet names):
  - summary_metrics                 (MEAN; manuscript main)
  - summary_metrics_MEDIAN          (MEDIAN; sensitivity / QA)
  - Axis_Reconstruction             (MEAN; a_ratio, b_ratio, SPI)
  - Axis_Reconstruction_MEDIAN      (MEDIAN; sensitivity / QA)
  - pI_Chem_tauPrimary              (chemical events; primary tau)
  - pI_Chem_tauSensitivity          (chemical events; tau grid)
  - UV_Lambda_vs_Dose               (UV-C cacao; λ(d))
  - UV_I_vs_Dose                    (UV-C cacao; I(d))
  - UV_ShapeEvent_Coffee10_tau1_10  (UV-C coffee; p-hat @ ~10 min)
  - QA_checks
  - metadata

Manuscript-locked definitions:
  - S_area = agg(Area_treated) / agg(Area_control)   (agg = mean by default; median also saved)
  - λ_area = -ln(S_area)
  - I_LWR  = ln( agg(LWR_treated) / agg(LWR_control) )

Ellipse-based axis reconstruction (paper/ipynb exact):
  - lwr_ratio = exp(I_LWR)
  - a_ratio   = sqrt(S_area * lwr_ratio)
  - b_ratio   = sqrt(S_area / lwr_ratio)
  - SPI       = (lwr_ratio - 1) * (1 - S_area)    # optional scalar index

Event probability:
  - Chemical: event if replicate LWR >= tau
  - UV-C:     event if replicate I_rep >= ln(tau), where I_rep = ln(LWR_rep / LWR0_mean)
  - p-hat with exact 95% Clopper–Pearson CI (SciPy if available; Wilson fallback otherwise)

Known UV isolate label fix:
  - "P24-55" is a known typo; canonical is "P23-55" (normalized automatically).

Usage (recommended):
  python scripts/chem_analysis_main.py \
    --chem-xlsx  data/Supplementary\ Data1_ISME\ Com.xlsx \
    --cacao-xlsx data/Cacao.xlsx \
    --coffee-xlsx data/Coffee.xlsx \
    --out-xlsx   outputs/SupplementaryData_S1.xlsx

If you place the input spreadsheets under ./data/ using the default filenames above,
you can omit the CLI arguments and just run:
  python scripts/chem_analysis_main.py
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import math
import os
import platform
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

# Exact Clopper–Pearson if SciPy exists; otherwise Wilson fallback
try:
    from scipy.stats import beta as _beta  # type: ignore
    _HAS_SCIPY = True
except Exception:
    _HAS_SCIPY = False


# -------------------------
# Defaults (manuscript-locked)
# -------------------------
DEFAULT_CHEM_TAU_PRIMARY = 1.10
DEFAULT_CHEM_TAU_SENS = (1.10, 1.15, 1.20)

DEFAULT_UV_TAU_FOR_I = 1.10
DEFAULT_TEN_MIN_WINDOW = (9.0, 11.0)  # inclusive window for "10 min cohort"


# -------------------------
# Column candidates (robust loading)
# -------------------------
CONTROL_CANDS_RAW = {
    "CONTROL", "CTRL", "ETOH", "ETHANOL", "MOCK", "DMSO", "WATER", "BLANK", "UNTREATED", "NEGATIVE"
}
ISOLATE_CANDS = ["Isolate", "isolate", "Strain", "Sample"]
TREATMENT_CANDS = ["Treatment", "treatment", "Chemical", "Chem", "Compound", "CompoundName", "Detailed treatment", "Specific"]
LWR_CANDS = ["LWR", "Length-to-width ratio(LWR)", "Length-to-width ratio", "L/W"]
LENGTH_CANDS = ["Length", "Length(L)[mm]", "Length(L)"]
WIDTH_CANDS = ["Width", "Width(W)[mm]", "Width(W)"]
AREA_CANDS = ["Area", "Area size", "Area size(AS)[mm2]", "Area(mm2)", "Area size (AS) [mm2]"]
TIME_CANDS = ["time(min)", "Time(min)", "Time (min)", "time_min", "Time_min", "Dose(min)", "Dose (min)", "Dose", "Exposure", "Exposure time"]

# Optional compound mapping (idempotent)
MR_NAME_MAP: Dict[str, str] = {
    "MR1001AM": "PhSOAM", "MR1001OH": "PhSOFA", "MR1015AM": "BCBGAM", "MR1015OH": "BCBGFA",
    "PHSOAM": "PhSOAM", "PHSOFA": "PhSOFA", "BCBGAM": "BCBGAM", "BCBGFA": "BCBGFA",
    "PhSOAM": "PhSOAM", "PhSOFA": "PhSOFA", "BCBGAM": "BCBGAM", "BCBGFA": "BCBGFA",
}
AMIDE_SET = {"PhSOAM", "BCBGAM"}
ACID_SET = {"PhSOFA", "BCBGFA"}

# UV isolate typo fix
ISOLATE_ALIAS = {"P24-55": "P23-55"}
DASHES = ["–", "—", "‑", "−"]  # en-dash, em-dash, non-breaking hyphen, minus


# -------------------------
# Utilities
# -------------------------
def _ensure_parent_dir(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def _norm_cols(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [str(c).strip().replace("\r", " ").replace("\n", " ") for c in df.columns]
    return df


def _pick_col(df: pd.DataFrame, candidates: Iterable[str]) -> Optional[str]:
    for c in candidates:
        if c in df.columns:
            return c
    low = {str(c).lower(): c for c in df.columns}
    for c in candidates:
        key = str(c).lower()
        if key in low:
            return low[key]
    return None


def _to_num(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def _agg(s: pd.Series, how: str) -> float:
    x = _to_num(s)
    how = how.lower().strip()
    if how == "mean":
        return float(x.mean(skipna=True))
    if how == "median":
        return float(x.median(skipna=True))
    raise ValueError(f"Unknown agg method: {how}")


def _map_treatment_name(x):
    if pd.isna(x):
        return x
    s = str(x).strip()
    key = s.replace(" ", "").upper()
    return MR_NAME_MAP.get(key, MR_NAME_MAP.get(s, s))


def _norm_controls_series(s: pd.Series) -> pd.Series:
    def f(x):
        if pd.isna(x):
            return x
        t = str(x).strip().upper()
        return "Control" if (t in CONTROL_CANDS_RAW or t == "CONTROL") else x
    return s.apply(f)


def normalize_isolate(x) -> str:
    if pd.isna(x):
        return x
    s = str(x).strip()
    for d in DASHES:
        s = s.replace(d, "-")
    s = re.sub(r"\s+", "", s)  # remove spaces like "P23 - 55"
    return ISOLATE_ALIAS.get(s, s)


def _parse_time_minutes(text) -> Optional[float]:
    """Best-effort parser for UV 'exposure time' from free text."""
    if pd.isna(text):
        return None
    s = str(text)
    if re.search(r"\b(control|ctrl)\b", s, flags=re.IGNORECASE):
        return 0.0
    m = re.search(r"(\d+(?:\.\d+)?)\s*(?:min|mins|minute|minutes)\b", s, flags=re.IGNORECASE)
    if m:
        return float(m.group(1))
    m = re.search(r"(\d+(?:\.\d+)?)\s*m\b", s, flags=re.IGNORECASE)
    if m:
        return float(m.group(1))
    return None


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            b = f.read(chunk_size)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def clopper_pearson_ci(k: int, n: int, alpha: float = 0.05) -> Tuple[float, float]:
    """
    Exact Clopper–Pearson CI via Beta quantiles if SciPy is available.
    Wilson score interval fallback if SciPy is not installed.
    """
    if n <= 0:
        return (np.nan, np.nan)

    if _HAS_SCIPY:
        low = 0.0 if k == 0 else float(_beta.ppf(alpha / 2, k, n - k + 1))
        high = 1.0 if k == n else float(_beta.ppf(1 - alpha / 2, k + 1, n - k))
        return (low, high)

    # Wilson fallback
    p = k / n
    z = 1.959963984540054
    den = 1 + z * z / n
    center = p + z * z / (2 * n)
    adj = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (max(0.0, (center - adj) / den), min(1.0, (center + adj) / den))


# -------------------------
# Loaders
# -------------------------
def load_chemical_dataframe(xlsx_path: Path, sheet: Optional[str] = None) -> pd.DataFrame:
    """
    Prefer sheet 'Morphology data' if present. Otherwise, search for a sheet with Isolate/Treatment columns.
    """
    xl = pd.ExcelFile(str(xlsx_path), engine="openpyxl")

    if sheet is not None:
        if sheet not in xl.sheet_names:
            raise RuntimeError(f"Requested chemical sheet '{sheet}' not found. Available: {xl.sheet_names}")
        df = pd.read_excel(xl, sheet_name=sheet)
        df = _norm_cols(df)
    elif "Morphology data" in xl.sheet_names:
        df = pd.read_excel(xl, sheet_name="Morphology data")
        df = _norm_cols(df)
    else:
        chosen = None
        for sh in xl.sheet_names:
            tmp = pd.read_excel(xl, sheet_name=sh)
            tmp = _norm_cols(tmp)
            iso = _pick_col(tmp, ISOLATE_CANDS)
            trt = _pick_col(tmp, TREATMENT_CANDS)
            if iso and trt:
                chosen = (sh, tmp)
                break
        if chosen is None:
            raise RuntimeError(f"Could not find a suitable chemical sheet. Sheets: {xl.sheet_names}")
        _, df = chosen

    iso_col = _pick_col(df, ISOLATE_CANDS) or "Isolate"
    trt_col = _pick_col(df, TREATMENT_CANDS) or "Treatment"
    if iso_col != "Isolate":
        df = df.rename(columns={iso_col: "Isolate"})
    if trt_col != "Treatment":
        df = df.rename(columns={trt_col: "Treatment"})

    df["Isolate"] = df["Isolate"].apply(normalize_isolate).astype(str).str.strip()
    df["Treatment"] = df["Treatment"].astype(str).str.strip()
    df["Treatment"] = _norm_controls_series(df["Treatment"])
    df["Treatment"] = df["Treatment"].apply(_map_treatment_name)

    # LWR (prefer existing; derive if needed)
    lwr_col = _pick_col(df, LWR_CANDS)
    if lwr_col and lwr_col != "LWR":
        df = df.rename(columns={lwr_col: "LWR"})

    len_col = _pick_col(df, LENGTH_CANDS)
    wid_col = _pick_col(df, WIDTH_CANDS)
    if len_col and len_col != "Length":
        df = df.rename(columns={len_col: "Length"})
    if wid_col and wid_col != "Width":
        df = df.rename(columns={wid_col: "Width"})
    if "LWR" not in df.columns and {"Length", "Width"}.issubset(df.columns):
        with np.errstate(divide="ignore", invalid="ignore"):
            df["LWR"] = _to_num(df["Length"]) / _to_num(df["Width"])

    # Area
    area_col = _pick_col(df, AREA_CANDS)
    if area_col and area_col != "Area":
        df = df.rename(columns={area_col: "Area"})

    if "Area" not in df.columns or "LWR" not in df.columns:
        raise RuntimeError("Chemical sheet must contain (or allow deriving) 'Area' and 'LWR'.")

    return df


def load_uv_dataframe(xlsx_path: Path) -> pd.DataFrame:
    df = pd.read_excel(str(xlsx_path), engine="openpyxl")
    df = _norm_cols(df)

    iso = _pick_col(df, ISOLATE_CANDS)
    if not iso:
        raise RuntimeError(f"UV file must have an isolate column. File: {xlsx_path}")
    if iso != "Isolate":
        df = df.rename(columns={iso: "Isolate"})
    df["Isolate"] = df["Isolate"].apply(normalize_isolate).astype(str).str.strip()

    # LWR
    lwr_col = _pick_col(df, LWR_CANDS)
    if lwr_col and lwr_col != "LWR":
        df = df.rename(columns={lwr_col: "LWR"})
    len_col = _pick_col(df, LENGTH_CANDS)
    wid_col = _pick_col(df, WIDTH_CANDS)
    if len_col and len_col != "Length":
        df = df.rename(columns={len_col: "Length"})
    if wid_col and wid_col != "Width":
        df = df.rename(columns={wid_col: "Width"})
    if "LWR" not in df.columns and {"Length", "Width"}.issubset(df.columns):
        with np.errstate(divide="ignore", invalid="ignore"):
            df["LWR"] = _to_num(df["Length"]) / _to_num(df["Width"])

    # Area
    area_col = _pick_col(df, AREA_CANDS)
    if area_col and area_col != "Area":
        df = df.rename(columns={area_col: "Area"})

    # Time (min)
    time_col = _pick_col(df, TIME_CANDS)
    if time_col:
        if time_col != "time_min":
            df = df.rename(columns={time_col: "time_min"})
        df["time_min"] = pd.to_numeric(df["time_min"], errors="coerce")
    else:
        src = _pick_col(df, ["Detailed treatment", "Treatment", "Specific"])
        if src is None:
            df["time_min"] = np.nan
        else:
            df["time_min"] = df[src].apply(_parse_time_minutes)

    if "Area" not in df.columns or "LWR" not in df.columns:
        raise RuntimeError(f"UV file must contain (or allow deriving) 'Area' and 'LWR'. File: {xlsx_path}")

    return df


# -------------------------
# Chemical computations
# -------------------------
def compute_chem_summary_and_axis(df: pd.DataFrame, how: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Returns:
      - summary_metrics table (lambda_area, I_LWR, S_area, n)
      - Axis reconstruction (a_ratio, b_ratio, SPI, etc.)
    """
    how = how.lower().strip()
    if how not in ("mean", "median"):
        raise ValueError("how must be 'mean' or 'median'")

    sum_rows: List[dict] = []
    axis_rows: List[dict] = []

    for iso, sub in df.groupby("Isolate", dropna=False):
        ctrl = sub[sub["Treatment"].astype(str).str.strip().str.upper().eq("CONTROL")]
        if len(ctrl) == 0:
            # Keep going; this should not happen for curated datasets
            continue

        A_c = _agg(ctrl["Area"], how)
        LWR_c = _agg(ctrl["LWR"], how)
        n_ctrl = int(_to_num(ctrl["LWR"]).notna().sum())

        for trt, g in sub.groupby("Treatment", dropna=False):
            if str(trt).strip().upper() == "CONTROL":
                continue

            A_t = _agg(g["Area"], how)
            LWR_t = _agg(g["LWR"], how)
            n_trt = int(_to_num(g["LWR"]).notna().sum())

            S_area = (A_t / A_c) if (pd.notna(A_t) and pd.notna(A_c) and A_c > 0) else np.nan
            lam = float(-math.log(S_area)) if (pd.notna(S_area) and S_area > 0) else np.nan
            Ival = float(math.log(LWR_t / LWR_c)) if (pd.notna(LWR_t) and pd.notna(LWR_c) and LWR_t > 0 and LWR_c > 0) else np.nan

            sum_rows.append({
                "isolate": iso,
                "treatment": trt,
                "lambda_area": lam,
                "I_LWR": Ival,
                "S_area": S_area,
                "n_trt": n_trt,
                "n_ctrl": n_ctrl,
                "agg": how,
            })

            # Axis reconstruction (paper/ipynb exact)
            if pd.notna(S_area) and pd.notna(Ival) and S_area > 0:
                lwr_ratio = float(math.exp(Ival))
                a_ratio = float(math.sqrt(S_area * lwr_ratio))
                b_ratio = float(math.sqrt(S_area / lwr_ratio)) if lwr_ratio > 0 else np.nan
                SPI = float((lwr_ratio - 1.0) * (1.0 - S_area))
            else:
                lwr_ratio = np.nan
                a_ratio = np.nan
                b_ratio = np.nan
                SPI = np.nan

            axis_rows.append({
                "isolate": iso,
                "treatment": trt,
                "a_ratio": a_ratio,
                "b_ratio": b_ratio,
                "SPI": SPI,
                "lwr_ratio": lwr_ratio,
                "S_area": S_area,
                "agg": how,
            })

    sm = pd.DataFrame(sum_rows)
    ax = pd.DataFrame(axis_rows)

    if len(sm) > 0:
        sm = sm.sort_values(["isolate", "treatment"]).reset_index(drop=True)
    if len(ax) > 0:
        ax = ax.sort_values(["isolate", "treatment"]).reset_index(drop=True)

    return sm, ax


def chem_event_tables(df: pd.DataFrame, tau_primary: float, tau_sens: Iterable[float]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Chemical events defined on absolute LWR: event if LWR >= tau.
    """
    rows_primary: List[dict] = []
    rows_sens: List[dict] = []

    taus = [float(t) for t in tau_sens]

    for (iso, trt), g in df.groupby(["Isolate", "Treatment"], dropna=False):
        if str(trt).strip().upper() == "CONTROL":
            continue

        x = _to_num(g["LWR"])
        n = int(x.notna().sum())
        if n <= 0:
            continue

        for tau in taus:
            k = int((x >= tau).sum())
            p_hat = k / n
            lo, hi = clopper_pearson_ci(k, n, alpha=0.05)
            rows_sens.append({
                "isolate": iso, "treatment": trt, "tau": tau,
                "n": n, "k": k, "p_hat": p_hat, "ci_low": lo, "ci_high": hi,
                "definition": "LWR >= tau",
            })

        tau0 = float(tau_primary)
        k0 = int((x >= tau0).sum())
        p0 = k0 / n
        lo0, hi0 = clopper_pearson_ci(k0, n, alpha=0.05)
        rows_primary.append({
            "isolate": iso, "treatment": trt, "tau": tau0,
            "n": n, "k": k0, "p_hat": p0, "ci_low": lo0, "ci_high": hi0,
            "definition": "LWR >= tau",
        })

    primary = pd.DataFrame(rows_primary)
    sens = pd.DataFrame(rows_sens)
    if len(primary) > 0:
        primary = primary.sort_values(["isolate", "treatment"]).reset_index(drop=True)
    if len(sens) > 0:
        sens = sens.sort_values(["tau", "isolate", "treatment"]).reset_index(drop=True)

    return primary, sens


def make_moa_legend(df: pd.DataFrame) -> pd.DataFrame:
    uniq = sorted(df["Treatment"].dropna().unique())
    rows = []
    for t in uniq:
        t_up = str(t).strip().upper()
        if t_up == "CONTROL":
            cls = "control"
        else:
            canon = _map_treatment_name(t)
            if canon in AMIDE_SET:
                cls = "amide"
            elif canon in ACID_SET:
                cls = "acid"
            else:
                cls = "other"
        rows.append({"treatment": t, "canonical": _map_treatment_name(t), "class": cls})
    out = pd.DataFrame(rows)
    if len(out) > 0:
        out = out.sort_values(["class", "canonical", "treatment"]).reset_index(drop=True)
    return out


# -------------------------
# UV computations
# -------------------------
def uv_lambda_I_vs_dose(df_uv: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    For each isolate:
      baseline at time_min == 0
      lambda(d) = -ln( mean(Area_t) / mean(Area_0) )
      I(d)      = ln( mean(LWR_t) / mean(LWR_0) )
    """
    out_lam: List[dict] = []
    out_I: List[dict] = []

    for iso, sub in df_uv.groupby("Isolate", dropna=False):
        base = sub[_to_num(sub["time_min"]).fillna(-1).eq(0.0)]
        if len(base) == 0:
            continue

        A0 = float(_to_num(base["Area"]).mean())
        LWR0 = float(_to_num(base["LWR"]).mean())
        if not (A0 > 0 and LWR0 > 0):
            continue

        for t, g in sub.groupby("time_min", dropna=True):
            t = float(t)
            if t == 0.0:
                continue

            At = float(_to_num(g["Area"]).mean())
            LWRt = float(_to_num(g["LWR"]).mean())

            lam = float(-math.log(At / A0)) if (At > 0) else np.nan
            Ival = float(math.log(LWRt / LWR0)) if (LWRt > 0) else np.nan

            out_lam.append({"isolate": iso, "dose_min": t, "lambda_area": lam})
            out_I.append({"isolate": iso, "dose_min": t, "I_LWR": Ival})

    lam_df = pd.DataFrame(out_lam)
    I_df = pd.DataFrame(out_I)

    if len(lam_df) > 0:
        lam_df = lam_df.sort_values(["isolate", "dose_min"]).reset_index(drop=True)
    else:
        lam_df = pd.DataFrame(columns=["isolate", "dose_min", "lambda_area"])

    if len(I_df) > 0:
        I_df = I_df.sort_values(["isolate", "dose_min"]).reset_index(drop=True)
    else:
        I_df = pd.DataFrame(columns=["isolate", "dose_min", "I_LWR"])

    return lam_df, I_df


def uv_coffee_shape_event_at_10min(df_uv: pd.DataFrame, tau: float, win: Tuple[float, float]) -> pd.DataFrame:
    """
    UV coffee event:
      I_rep = ln(LWR_rep / LWR0_mean)
      event if I_rep >= ln(tau)
      exact Clopper–Pearson CI
    """
    lo, hi = float(win[0]), float(win[1])
    I_star = float(math.log(float(tau)))
    rows: List[dict] = []

    for iso, sub in df_uv.groupby("Isolate", dropna=False):
        base = sub[_to_num(sub["time_min"]).fillna(-1).eq(0.0)]
        if len(base) == 0:
            continue

        LWR0 = float(_to_num(base["LWR"]).mean())
        if not (LWR0 > 0):
            continue

        tvals = _to_num(sub["time_min"])
        at10 = sub[(tvals >= lo) & (tvals <= hi)]
        if len(at10) == 0:
            continue

        Irep = np.log(_to_num(at10["LWR"]) / LWR0)
        n = int(Irep.notna().sum())
        k = int((Irep >= I_star).sum())

        p_hat = (k / n) if n > 0 else np.nan
        lo_ci, hi_ci = clopper_pearson_ci(k, n, alpha=0.05) if n > 0 else (np.nan, np.nan)

        rows.append({
            "isolate": iso,
            "n": n, "k": k, "p_hat": p_hat, "ci_low": lo_ci, "ci_high": hi_ci,
            "tau": float(tau), "I_star": float(I_star),
            "time_window": f"{lo:g}–{hi:g} min",
            "definition": "I_rep >= ln(tau)",
        })

    out = pd.DataFrame(rows)
    if len(out) > 0:
        out = out.sort_values("isolate").reset_index(drop=True)
    else:
        out = pd.DataFrame(columns=["isolate", "n", "k", "p_hat", "ci_low", "ci_high", "tau", "I_star", "time_window", "definition"])
    return out


# -------------------------
# QA + metadata
# -------------------------
def qa_checks(sm_mean: pd.DataFrame) -> pd.DataFrame:
    if sm_mean is None or len(sm_mean) == 0:
        return pd.DataFrame([{"metric": "summary_metrics_rows", "value": 0}])

    df = sm_mean.copy()
    df["lambda_area"] = pd.to_numeric(df["lambda_area"], errors="coerce")
    df["I_LWR"] = pd.to_numeric(df["I_LWR"], errors="coerce")

    df2 = df.dropna(subset=["lambda_area", "I_LWR"])
    r = float(df2["lambda_area"].corr(df2["I_LWR"])) if len(df2) >= 2 else np.nan

    max_row = df.dropna(subset=["lambda_area"]).sort_values("lambda_area", ascending=False).head(1)
    if len(max_row) == 1:
        max_lam = float(max_row["lambda_area"].iloc[0])
        max_id = f"{max_row['isolate'].iloc[0]}–{max_row['treatment'].iloc[0]}"
    else:
        max_lam = np.nan
        max_id = ""

    return pd.DataFrame([
        {"metric": "rows_summary_metrics_mean", "value": int(len(df))},
        {"metric": "unique_isolates", "value": int(df["isolate"].nunique())},
        {"metric": "unique_treatments", "value": int(df["treatment"].nunique())},
        {"metric": "pearson_r(lambda,I)", "value": r},
        {"metric": "max_lambda", "value": max_lam},
        {"metric": "max_lambda_isolate_treatment", "value": max_id},
        {"metric": "scipy_exact_cp", "value": str(_HAS_SCIPY)},
    ])


def build_metadata(
    chem_path: Path,
    cacao_path: Optional[Path],
    coffee_path: Optional[Path],
    out_xlsx: Path,
    tau_primary: float,
    tau_sens: Iterable[float],
    uv_tau_for_I: Optional[float],
    uv_win: Optional[Tuple[float, float]],
) -> pd.DataFrame:
    rows = []

    def add_file(prefix: str, p: Optional[Path]) -> None:
        if p is None:
            rows.append((f"{prefix}_file", ""))
            rows.append((f"{prefix}_sha256", ""))
            return
        rows.append((f"{prefix}_file", str(p)))
        try:
            rows.append((f"{prefix}_sha256", sha256_file(p)))
        except Exception:
            rows.append((f"{prefix}_sha256", ""))

    rows.append(("created_at", dt.datetime.now().isoformat(timespec="seconds")))
    rows.append(("python", sys.version.replace("\n", " ")))
    rows.append(("platform", platform.platform()))
    rows.append(("pandas_version", pd.__version__))
    rows.append(("numpy_version", np.__version__))
    rows.append(("scipy_exact_cp", str(_HAS_SCIPY)))

    add_file("chem", chem_path)
    add_file("cacao_uv", cacao_path)
    add_file("coffee_uv", coffee_path)

    rows.append(("out_xlsx", str(out_xlsx)))
    rows.append(("chem_tau_primary", str(float(tau_primary))))
    rows.append(("chem_tau_sensitivity", ",".join(map(lambda x: f"{float(x):g}", tau_sens))))

    if uv_tau_for_I is not None:
        rows.append(("uv_tau_for_I", str(float(uv_tau_for_I))))
    if uv_win is not None:
        rows.append(("ten_min_window", f"{uv_win[0]:g}–{uv_win[1]:g}"))

    rows.append(("uv_isolate_alias_fix", "P24-55→P23-55"))

    return pd.DataFrame(rows, columns=["key", "value"])


# -------------------------
# CLI
# -------------------------
def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    here = Path(__file__).resolve()
    repo_root = here.parents[1]
    data_dir = repo_root / "data"
    outputs_dir = repo_root / "outputs"

    p = argparse.ArgumentParser(
        prog="chem_analysis_main.py",
        description="Generate SupplementaryData_S1.xlsx (λ–I–p framework) from raw spreadsheets.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    p.add_argument("--chem-xlsx", type=Path, default=(data_dir / "Supplementary Data1_ISME Com.xlsx"),
                   help="Chemical morphology spreadsheet (raw).")
    p.add_argument("--chem-sheet", type=str, default=None,
                   help="Optional: specify chemical sheet name explicitly. If omitted, auto-detect.")
    p.add_argument("--cacao-xlsx", type=Path, default=(data_dir / "Cacao.xlsx"),
                   help="UV-C cacao spreadsheet (raw).")
    p.add_argument("--coffee-xlsx", type=Path, default=(data_dir / "Coffee.xlsx"),
                   help="UV-C coffee spreadsheet (raw).")

    p.add_argument("--outdir", type=Path, default=outputs_dir,
                   help="Output directory (used if --out-xlsx is not provided).")

    p.add_argument("--out-xlsx", type=Path, default=None,
                   help="Output workbook path (overrides --outdir).")

    p.add_argument("--chem-tau-primary", type=float, default=DEFAULT_CHEM_TAU_PRIMARY,
                   help="Primary chemical event threshold tau (LWR >= tau).")
    p.add_argument("--chem-tau-sens", type=str, default=",".join(map(str, DEFAULT_CHEM_TAU_SENS)),
                   help="Comma-separated tau grid for sensitivity (chemical).")

    p.add_argument("--uv-tau", type=float, default=DEFAULT_UV_TAU_FOR_I,
                   help="UV-C event threshold tau for I_rep >= ln(tau). (Coffee ~10 min cohort)")
    p.add_argument("--uv-win", type=str, default=f"{DEFAULT_TEN_MIN_WINDOW[0]},{DEFAULT_TEN_MIN_WINDOW[1]}",
                   help="Coffee 10-min cohort window as 'lo,hi' (minutes).")

    p.add_argument("--skip-uv", action="store_true",
                   help="Skip UV-C reanalysis (writes chemical sheets only).")
    p.add_argument("--skip-median", action="store_true",
                   help="Skip median-based summary sheets (saves time; mean still written).")

    return p.parse_args(argv)


# -------------------------
# Main
# -------------------------
def main(argv: Optional[List[str]] = None) -> None:
    args = parse_args(argv)

    chem_path: Path = args.chem_xlsx
    cacao_path: Path = args.cacao_xlsx
    coffee_path: Path = args.coffee_xlsx
    out_xlsx: Path = (args.out_xlsx if args.out_xlsx is not None else (args.outdir / "SupplementaryData_S1.xlsx"))

    if not chem_path.exists():
        raise FileNotFoundError(
            f"Chemical spreadsheet not found: {chem_path}\n"
            f"Tip: place it under ./data/ or pass --chem-xlsx PATH"
        )

    tau_sens = [float(x) for x in str(args.chem_tau_sens).split(",") if str(x).strip() != ""]
    uv_win_parts = [float(x) for x in str(args.uv_win).split(",")]
    if len(uv_win_parts) != 2:
        raise ValueError("--uv-win must be 'lo,hi' (e.g., '9,11').")
    uv_win = (float(uv_win_parts[0]), float(uv_win_parts[1]))

    _ensure_parent_dir(out_xlsx)

    print(f"[INFO] Loading chemical data: {chem_path}")
    chem_df = load_chemical_dataframe(chem_path, sheet=args.chem_sheet)
    print(f"[INFO] Chemical rows: {len(chem_df):,} | isolates: {chem_df['Isolate'].nunique()} | treatments: {chem_df['Treatment'].nunique()}")

    print("[INFO] Summarizing chemical dataset (MEAN)...")
    sm_mean, ax_mean = compute_chem_summary_and_axis(chem_df, how="mean")

    if args.skip_median:
        sm_med = pd.DataFrame()
        ax_med = pd.DataFrame()
    else:
        print("[INFO] Summarizing chemical dataset (MEDIAN)...")
        sm_med, ax_med = compute_chem_summary_and_axis(chem_df, how="median")

    legend = make_moa_legend(chem_df)

    print("[INFO] Computing chemical event tables...")
    p_primary, p_sens = chem_event_tables(chem_df, tau_primary=float(args.chem_tau_primary), tau_sens=tau_sens)

    qa = qa_checks(sm_mean)

    # UV (optional)
    cacao_df = None
    coffee_df = None
    lam_uv = pd.DataFrame()
    I_uv = pd.DataFrame()
    coffee_event = pd.DataFrame()

    if not args.skip_uv:
        if not cacao_path.exists():
            raise FileNotFoundError(
                f"UV-C cacao spreadsheet not found: {cacao_path}\n"
                f"Tip: place it under ./data/ or pass --cacao-xlsx PATH, or use --skip-uv."
            )
        if not coffee_path.exists():
            raise FileNotFoundError(
                f"UV-C coffee spreadsheet not found: {coffee_path}\n"
                f"Tip: place it under ./data/ or pass --coffee-xlsx PATH, or use --skip-uv."
            )

        print(f"[INFO] Loading UV-C cacao data: {cacao_path}")
        cacao_df = load_uv_dataframe(cacao_path)
        print(f"[INFO] Loading UV-C coffee data: {coffee_path}")
        coffee_df = load_uv_dataframe(coffee_path)

        print("[INFO] Summarizing UV datasets...")
        lam_uv, I_uv = uv_lambda_I_vs_dose(cacao_df)
        coffee_event = uv_coffee_shape_event_at_10min(coffee_df, tau=float(args.uv_tau), win=uv_win)

    meta = build_metadata(
        chem_path=chem_path,
        cacao_path=None if args.skip_uv else cacao_path,
        coffee_path=None if args.skip_uv else coffee_path,
        out_xlsx=out_xlsx,
        tau_primary=float(args.chem_tau_primary),
        tau_sens=tau_sens,
        uv_tau_for_I=None if args.skip_uv else float(args.uv_tau),
        uv_win=None if args.skip_uv else uv_win,
    )

    # Write workbook
    with pd.ExcelWriter(str(out_xlsx), engine="openpyxl") as xw:
        # manuscript/main tables
        sm_mean.to_excel(xw, index=False, sheet_name="summary_metrics")
        ax_mean.to_excel(xw, index=False, sheet_name="Axis_Reconstruction")
        p_primary.to_excel(xw, index=False, sheet_name="pI_Chem_tauPrimary")
        p_sens.to_excel(xw, index=False, sheet_name="pI_Chem_tauSensitivity")

        # sensitivity/QA
        if not args.skip_median and len(sm_med) > 0:
            sm_med.to_excel(xw, index=False, sheet_name="summary_metrics_MEDIAN")
        if not args.skip_median and len(ax_med) > 0:
            ax_med.to_excel(xw, index=False, sheet_name="Axis_Reconstruction_MEDIAN")

        legend.to_excel(xw, index=False, sheet_name="moa_vector_index_legend")

        if not args.skip_uv:
            lam_uv.to_excel(xw, index=False, sheet_name="UV_Lambda_vs_Dose")
            I_uv.to_excel(xw, index=False, sheet_name="UV_I_vs_Dose")
            coffee_event.to_excel(xw, index=False, sheet_name="UV_ShapeEvent_Coffee10_tau1_10")

        qa.to_excel(xw, index=False, sheet_name="QA_checks")
        meta.to_excel(xw, index=False, sheet_name="metadata")

    print(f"[DONE] Wrote: {out_xlsx}")
    if len(qa) > 0:
        print("[QA] Key values:")
        for _, r in qa.iterrows():
            print(f" - {r['metric']}: {r['value']}")


if __name__ == "__main__":
    main()
