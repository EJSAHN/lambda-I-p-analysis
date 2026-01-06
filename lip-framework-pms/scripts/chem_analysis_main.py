# -*- coding: utf-8 -*-
"""
chem_analysis_main.py (FINAL)

- Reproduces manuscript λ–I–p framework tables into outputs/SupplementaryData_S1.xlsx
- Chemical metrics saved as:
    - summary_metrics (MEAN; manuscript main)
    - summary_metrics_MEDIAN (sensitivity / QA)
- Axis reconstruction saved as (paper/ipynb algebraic reconstruction):
    - Axis_Reconstruction (MEAN; a_ratio, b_ratio, SPI)
    - Axis_Reconstruction_MEDIAN
- Chemical event probabilities:
    - pI_Chem_tauPrimary (tau=1.10)
    - pI_Chem_tauSensitivity (tau grid)
- UV-C reanalysis:
    - UV_Lambda_vs_Dose
    - UV_I_vs_Dose
    - UV_ShapeEvent_Coffee10_tau1_10  (I >= ln(tau), exact Clopper–Pearson CI)
- IMPORTANT: UV isolate label fix:
    P24-55 is a known typo; canonical is P23-55. We normalize automatically.

Run:
    # 1. Place data files (Chem.xlsx, Cacao.xlsx, Coffee.xlsx) in the 'data/' folder
    # 2. Run the script:
    python scripts/chem_analysis_main.py
"""

from __future__ import annotations

import os
import re
import math
import datetime as dt
from typing import Optional, Iterable, Dict, Tuple, List

import numpy as np
import pandas as pd

# Exact Clopper–Pearson if scipy exists; otherwise Wilson fallback
try:
    from scipy.stats import beta as _beta  # type: ignore
    _HAS_SCIPY = True
except Exception:
    _HAS_SCIPY = False

# =========================
# Configuration
# =========================
# 1. Setup Paths (Relative to this script)
# Assumes structure: repo/scripts/chem_analysis_main.py -> repo/data/
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT  = os.path.dirname(SCRIPT_DIR) 

DATA_DIR = os.path.join(REPO_ROOT, "data")
OUTDIR   = os.path.join(REPO_ROOT, "outputs")

# 2. Input Files (Clean names for GitHub)
# NOTE: Please rename your local file "Supplementary Data1_ISME Com.xlsx" to "Chem.xlsx"
CHEM_XLSX   = os.path.join(DATA_DIR, "Chem.xlsx")
CACAO_XLSX  = os.path.join(DATA_DIR, "Cacao.xlsx")
COFFEE_XLSX = os.path.join(DATA_DIR, "Coffee.xlsx")

# 3. Output File
OUT_XLSX = os.path.join(OUTDIR, "SupplementaryData_S1.xlsx")

# Chemical event thresholds (absolute LWR)
CHEM_TAU_PRIMARY = 1.10
CHEM_TAU_SENS    = [1.10, 1.15, 1.20]

# UV-C: coffee ~10 min event defined on I >= ln(tau) with tau=1.10
UV_TAU_FOR_I   = 1.10
TEN_MIN_WINDOW = (9.0, 11.0)  # inclusive

# Control candidates
CONTROL_CANDS_RAW = [
    "CONTROL", "CTRL", "ETOH", "ETHANOL", "MOCK", "DMSO", "WATER", "BLANK", "UNTREATED", "NEGATIVE"
]

# Column candidates (robust loading)
ISOLATE_CANDS   = ["Isolate", "isolate", "Strain", "Sample"]
TREATMENT_CANDS = ["Treatment", "treatment", "Chemical", "Chem", "Compound", "CompoundName", "Detailed treatment", "Specific"]
LWR_CANDS       = ["LWR", "Length-to-width ratio(LWR)", "Length-to-width ratio", "L/W"]
LENGTH_CANDS    = ["Length", "Length(L)[mm]", "Length(L)"]
WIDTH_CANDS     = ["Width", "Width(W)[mm]", "Width(W)"]
AREA_CANDS      = ["Area", "Area size", "Area size(AS)[mm2]", "Area(mm2)", "Area size (AS) [mm2]"]

TIME_CANDS      = ["time(min)", "Time(min)", "Time (min)", "time_min", "Time_min", "Dose(min)", "Dose (min)", "Dose", "Exposure", "Exposure time"]

# MR code mapping (optional; safe/idempotent)
MR_NAME_MAP: Dict[str, str] = {
    "MR1001AM":"PhSOAM", "MR1001OH":"PhSOFA", "MR1015AM":"BCBGAM", "MR1015OH":"BCBGFA",
    "PHSOAM":"PhSOAM", "PHSOFA":"PhSOFA", "BCBGAM":"BCBGAM", "BCBGFA":"BCBGFA",
    "PhSOAM":"PhSOAM", "PhSOFA":"PhSOFA", "BCBGAM":"BCBGAM", "BCBGFA":"BCBGFA"
}
AMIDE_SET = {"PhSOAM", "BCBGAM"}
ACID_SET  = {"PhSOFA", "BCBGFA"}

# ---- UV isolate typo fix ----
ISOLATE_ALIAS = {
    "P24-55": "P23-55",  # known typo in UV dataset
}
DASHES = ["–", "—", "‑", "−"]  # en-dash, em-dash, non-breaking hyphen, minus

# =========================
# Utilities
# =========================
def _ensure_outdir(p: str) -> None:
    os.makedirs(p, exist_ok=True)

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
        if str(c).lower() in low:
            return low[str(c).lower()]
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

def clopper_pearson_ci(k: int, n: int, alpha: float = 0.05) -> Tuple[float, float]:
    if n <= 0:
        return (np.nan, np.nan)
    if _HAS_SCIPY:
        low  = 0.0 if k == 0 else float(_beta.ppf(alpha/2, k, n-k+1))
        high = 1.0 if k == n else float(_beta.ppf(1-alpha/2, k+1, n-k))
        return (low, high)
    # Wilson fallback (rare)
    p = k / n
    z = 1.959963984540054
    den = 1 + z*z/n
    center = p + z*z/(2*n)
    adj = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n))
    return (max(0.0, (center-adj)/den), min(1.0, (center+adj)/den))

# =========================
# Loaders
# =========================
def load_chemical_dataframe(xlsx_path: str) -> pd.DataFrame:
    """
    Prefer sheet 'Morphology data' if present. Otherwise, search for a sheet with Isolate/Treatment columns.
    """
    xl = pd.ExcelFile(xlsx_path, engine="openpyxl")

    if "Morphology data" in xl.sheet_names:
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
            raise RuntimeError(f"Could not find suitable chemical sheet. Sheets: {xl.sheet_names}")
        sh, df = chosen

    # normalize isolate + treatment columns
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

    # Ensure LWR exists
    lwr_col = _pick_col(df, LWR_CANDS)
    if lwr_col and lwr_col != "LWR":
        df = df.rename(columns={lwr_col: "LWR"})
    len_col = _pick_col(df, LENGTH_CANDS)
    wid_col = _pick_col(df, WIDTH_CANDS)
    if len_col and len_col != "Length":
        df = df.rename(columns={len_col: "Length"})
    if wid_col and wid_col != "Width":
        df = df.rename(columns={wid_col: "Width"})
    if "LWR" not in df.columns and {"Length","Width"}.issubset(df.columns):
        with np.errstate(divide="ignore", invalid="ignore"):
            df["LWR"] = _to_num(df["Length"]) / _to_num(df["Width"])

    # Ensure Area exists
    area_col = _pick_col(df, AREA_CANDS)
    if area_col and area_col != "Area":
        df = df.rename(columns={area_col: "Area"})

    if "Area" not in df.columns or "LWR" not in df.columns:
        raise RuntimeError("Chemical sheet must contain (or allow deriving) 'Area' and 'LWR'.")

    print(f"[INFO] Chemical sheet selected. Columns: {list(df.columns)}")
    return df

def load_uv_dataframe(xlsx_path: str) -> pd.DataFrame:
    df = pd.read_excel(xlsx_path, engine="openpyxl")
    df = _norm_cols(df)

    iso = _pick_col(df, ISOLATE_CANDS)
    if not iso:
        raise RuntimeError(f"UV file must have isolate column. File: {xlsx_path}")
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
    if "LWR" not in df.columns and {"Length","Width"}.issubset(df.columns):
        with np.errstate(divide="ignore", invalid="ignore"):
            df["LWR"] = _to_num(df["Length"]) / _to_num(df["Width"])

    # Area
    area_col = _pick_col(df, AREA_CANDS)
    if area_col and area_col != "Area":
        df = df.rename(columns={area_col: "Area"})

    # Time
    time_col = _pick_col(df, TIME_CANDS)
    if time_col:
        if time_col != "time_min":
            df = df.rename(columns={time_col: "time_min"})
        df["time_min"] = pd.to_numeric(df["time_min"], errors="coerce")
    else:
        # parse from likely string columns
        src = _pick_col(df, ["Detailed treatment", "Treatment", "Specific"])
        if src is None:
            df["time_min"] = np.nan
        else:
            df["time_min"] = df[src].apply(_parse_time_minutes)

    # for any remaining NaNs, attempt fallback scan for control keyword
    if df["time_min"].isna().any():
        mask = df["time_min"].isna()
        df.loc[mask, "time_min"] = df.loc[mask].apply(
            lambda r: 0.0 if re.search(r"(control|ctrl)", " ".join(map(str, r.values)), re.I) else np.nan,
            axis=1
        )

    if "LWR" not in df.columns or "Area" not in df.columns:
        raise RuntimeError(f"UV file must contain (or allow deriving) 'LWR' and 'Area'. File: {xlsx_path}")

    return df

# =========================
# Chemical computations
# =========================
def compute_chem_summary_and_axis(df: pd.DataFrame, how: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Manuscript-locked definitions:
      S_area = mean(treated Area) / mean(control Area)   (or median if how='median')
      lambda = -ln(S_area)
      I      = ln(LWR_t / LWR_c)

    Axis reconstruction (paper/ipynb):
      lwr_ratio = exp(I)
      a_ratio = sqrt(S_area * lwr_ratio)
      b_ratio = sqrt(S_area / lwr_ratio)
      SPI = (lwr_ratio - 1) * (1 - S_area)
    """
    how = how.lower().strip()
    if how not in ("mean", "median"):
        raise ValueError("how must be 'mean' or 'median'")

    sum_rows: List[dict] = []
    axis_rows: List[dict] = []

    for iso, sub in df.groupby("Isolate", dropna=False):
        ctrl = sub[sub["Treatment"].astype(str).str.strip().str.upper().eq("CONTROL")]
        if len(ctrl) == 0:
            print(f"[WARN] Chemical isolate '{iso}': no control rows; skipped.")
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

            # Survival ratio
            S_area = (A_t / A_c) if (not np.isnan(A_t) and not np.isnan(A_c) and A_c > 0) else np.nan

            # lambda
            lam = float(-math.log(S_area)) if (not np.isnan(S_area) and S_area > 0) else np.nan

            # I
            Ival = float(math.log(LWR_t / LWR_c)) if (not np.isnan(LWR_t) and not np.isnan(LWR_c) and LWR_t > 0 and LWR_c > 0) else np.nan

            sum_rows.append({
                "isolate": iso,
                "treatment": trt,
                "lambda_area": lam,
                "I_LWR": Ival,
                "S_area": S_area,
                "n_trt": n_trt,
                "n_ctrl": n_ctrl,
                "agg": how
            })

            # Axis reconstruction (paper/ipynb exact)
            if (not np.isnan(S_area)) and (not np.isnan(Ival)) and (S_area > 0):
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
                "agg": how
            })

    sm = pd.DataFrame(sum_rows)
    ax = pd.DataFrame(axis_rows)

    if len(sm) > 0:
        sm = sm.sort_values(["isolate", "treatment"]).reset_index(drop=True)
    if len(ax) > 0:
        ax = ax.sort_values(["isolate", "treatment"]).reset_index(drop=True)

    return sm, ax

def chem_event_tables(df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Chemical events defined on absolute LWR:
      event if LWR >= tau
    """
    rows_primary: List[dict] = []
    rows_sens: List[dict] = []

    for (iso, trt), g in df.groupby(["Isolate","Treatment"], dropna=False):
        if str(trt).strip().upper() == "CONTROL":
            continue

        x = _to_num(g["LWR"])
        n = int(x.notna().sum())
        if n <= 0:
            continue

        for tau in CHEM_TAU_SENS:
            k = int((x >= float(tau)).sum())
            p_hat = k / n
            lo, hi = clopper_pearson_ci(k, n, alpha=0.05)
            rows_sens.append({
                "isolate": iso, "treatment": trt,
                "tau": float(tau),
                "n": n, "k": k,
                "p_hat": p_hat,
                "ci_low": lo, "ci_high": hi,
                "definition": "LWR >= tau"
            })

        # primary tau
        tau0 = float(CHEM_TAU_PRIMARY)
        k0 = int((x >= tau0).sum())
        p0 = k0 / n
        lo0, hi0 = clopper_pearson_ci(k0, n, alpha=0.05)
        rows_primary.append({
            "isolate": iso, "treatment": trt,
            "tau": tau0,
            "n": n, "k": k0,
            "p_hat": p0,
            "ci_low": lo0, "ci_high": hi0,
            "definition": "LWR >= tau"
        })

    primary = pd.DataFrame(rows_primary)
    sens = pd.DataFrame(rows_sens)

    if len(primary) > 0:
        primary = primary.sort_values(["isolate","treatment"]).reset_index(drop=True)
    if len(sens) > 0:
        sens = sens.sort_values(["tau","isolate","treatment"]).reset_index(drop=True)

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
        out = out.sort_values(["class","canonical","treatment"]).reset_index(drop=True)
    return out

# =========================
# UV computations
# =========================
def uv_lambda_I_vs_dose(df_uv: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    For each isolate:
      baseline at time_min == 0
      lambda(d) = -ln( mean(Area_t) / mean(Area_0) )
      I(d)      =  ln( mean(LWR_t)  / mean(LWR_0)  )
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
        lam_df = lam_df.sort_values(["isolate","dose_min"]).reset_index(drop=True)
    else:
        lam_df = pd.DataFrame(columns=["isolate","dose_min","lambda_area"])
    if len(I_df) > 0:
        I_df = I_df.sort_values(["isolate","dose_min"]).reset_index(drop=True)
    else:
        I_df = pd.DataFrame(columns=["isolate","dose_min","I_LWR"])
    return lam_df, I_df

def uv_coffee_shape_event_at_10min(df_uv: pd.DataFrame, tau: float, win=(9.0,11.0)) -> pd.DataFrame:
    """
    UV coffee event:
      I_rep = ln(LWR_rep / LWR0_mean)
      event if I_rep >= ln(tau)
      exact Clopper–Pearson CI
    """
    lo, hi = float(win[0]), float(win[1])
    I_star = float(math.log(tau))

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
            "n": n, "k": k, "p_hat": p_hat,
            "ci_low": lo_ci, "ci_high": hi_ci,
            "tau": float(tau),
            "I_star": float(I_star),
            "time_window": f"{lo:g}–{hi:g} min"
        })

    out = pd.DataFrame(rows)
    if len(out) > 0:
        out = out.sort_values("isolate").reset_index(drop=True)
    else:
        out = pd.DataFrame(columns=["isolate","n","k","p_hat","ci_low","ci_high","tau","I_star","time_window"])
    return out

# =========================
# QA + Metadata
# =========================
def qa_checks(sm_mean: pd.DataFrame) -> pd.DataFrame:
    """
    Quick sanity checks to help you confirm manuscript consistency (e.g., max lambda ~0.67).
    """
    if sm_mean is None or len(sm_mean) == 0:
        return pd.DataFrame([{"metric": "summary_metrics_rows", "value": 0}])

    df = sm_mean.copy()
    df["lambda_area"] = pd.to_numeric(df["lambda_area"], errors="coerce")
    df["I_LWR"] = pd.to_numeric(df["I_LWR"], errors="coerce")

    df2 = df.dropna(subset=["lambda_area","I_LWR"])
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

# =========================
# Main
# =========================
def main() -> None:
    _ensure_outdir(OUTDIR)

    print(f"[INFO] Loading chemical data: {CHEM_XLSX}")
    chem_df = load_chemical_dataframe(CHEM_XLSX)

    print(f"[INFO] Loading UV-C Cacao data: {CACAO_XLSX}")
    cacao_df = load_uv_dataframe(CACAO_XLSX)

    print(f"[INFO] Loading UV-C Coffee data: {COFFEE_XLSX}")
    coffee_df = load_uv_dataframe(COFFEE_XLSX)

    print("[INFO] Summarizing chemical dataset (MEAN and MEDIAN)...")
    sm_mean, ax_mean = compute_chem_summary_and_axis(chem_df, how="mean")
    sm_med,  ax_med  = compute_chem_summary_and_axis(chem_df, how="median")

    legend = make_moa_legend(chem_df)

    print("[INFO] Computing chemical event tables...")
    p_primary, p_sens = chem_event_tables(chem_df)

    print("[INFO] Summarizing UV datasets...")
    lam_uv, I_uv = uv_lambda_I_vs_dose(cacao_df)
    coffee_event = uv_coffee_shape_event_at_10min(coffee_df, tau=UV_TAU_FOR_I, win=TEN_MIN_WINDOW)

    qa = qa_checks(sm_mean)

    meta = pd.DataFrame({
        "key": [
            "created_at",
            "chem_file", "cacao_file", "coffee_file",
            "chem_tau_primary", "chem_tau_sensitivity",
            "uv_tau_for_I", "ten_min_window",
            "uv_isolate_alias_fix",
            "scipy_exact_CP"
        ],
        "value": [
            dt.datetime.now().isoformat(timespec="seconds"),
            CHEM_XLSX, CACAO_XLSX, COFFEE_XLSX,
            str(CHEM_TAU_PRIMARY),
            ",".join(map(str, CHEM_TAU_SENS)),
            str(UV_TAU_FOR_I),
            f"{TEN_MIN_WINDOW[0]}–{TEN_MIN_WINDOW[1]}",
            "P24-55→P23-55",
            str(_HAS_SCIPY),
        ]
    })

    with pd.ExcelWriter(OUT_XLSX, engine="openpyxl") as xw:
        # manuscript/main tables
        sm_mean.to_excel(xw, index=False, sheet_name="summary_metrics")
        ax_mean.to_excel(xw, index=False, sheet_name="Axis_Reconstruction")
        p_primary.to_excel(xw, index=False, sheet_name="pI_Chem_tauPrimary")
        p_sens.to_excel(xw, index=False, sheet_name="pI_Chem_tauSensitivity")

        # sensitivity/QA
        sm_med.to_excel(xw, index=False, sheet_name="summary_metrics_MEDIAN")
        ax_med.to_excel(xw, index=False, sheet_name="Axis_Reconstruction_MEDIAN")

        # legend
        legend.to_excel(xw, index=False, sheet_name="moa_vector_index_legend")

        # UV sheets used for figures
        lam_uv.to_excel(xw, index=False, sheet_name="UV_Lambda_vs_Dose")
        I_uv.to_excel(xw, index=False, sheet_name="UV_I_vs_Dose")
        coffee_event.to_excel(xw, index=False, sheet_name="UV_ShapeEvent_Coffee10_tau1_10")

        # QA/meta
        qa.to_excel(xw, index=False, sheet_name="QA_checks")
        meta.to_excel(xw, index=False, sheet_name="metadata")

    print(f"[DONE] Wrote: {OUT_XLSX}")
    if len(qa) > 0:
        print("[QA] Key values:")
        for _, r in qa.iterrows():
            print(f"  - {r['metric']}: {r['value']}")

if __name__ == "__main__":
    # Safety Check: Input files exist?
    if not os.path.exists(CHEM_XLSX):
        print(f"[ERROR] Input file not found: {CHEM_XLSX}")
        print("Please place 'Chem.xlsx', 'Cacao.xlsx', and 'Coffee.xlsx' in the 'data/' folder.")
        exit(1)
        
    main()