# -*- coding: utf-8 -*-
"""
chem_analysis_main.py

Generate a single summary workbook (SupplementaryData_S1.xlsx) from externally supplied
colony-morphology spreadsheets and derive all downstream summary tables used by the
λ–I–p analysis framework.

Key outputs (sheet names):
  - summary_metrics                 (MEAN; manuscript main)
  - summary_metrics_MEDIAN          (MEDIAN; sensitivity / QA)
  - Axis_Reconstruction             (MEAN; a_ratio, b_ratio, SPI)
  - Axis_Reconstruction_MEDIAN      (MEDIAN; sensitivity / QA)
  - pI_Chem_tauPrimary              (chemical events; primary tau)
  - pI_Chem_tauSensitivity          (chemical events; tau grid)
  - pI_Chem_tauWide                 (chemical events; wide-form key indicators across tau)
  - pI_Control_tauSensitivity       (control baseline events; tau grid)
  - pI_Control_tauWide              (control baseline events; wide-form key indicators across tau)
  - Supplementary_Table_S1          (near-equal potency, divergent polarity examples)
  - Axis_Reconstruction_Validation  (reconstructed vs direct axis ratios)
  - Axis_Recon_Val_Summary          (validation summary metrics)
  - UV_Lambda_vs_Dose               (UV-C cacao; λ(d))
  - UV_I_vs_Dose                    (UV-C cacao; I(d))
  - UV_ShapeEvent_Coffee10_tau1_10  (UV-C coffee; p-hat @ ~10 min)
  - QA_checks
  - metadata

Definitions:
  - S_area = agg(Area_treated) / agg(Area_control)   (agg = mean by default; median also saved)
  - λ_area = -ln(S_area)
  - I_LWR  = ln( agg(LWR_treated) / agg(LWR_control) )

Ellipse-based axis reconstruction:
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

Usage:
  python scripts/chem_analysis_main.py \
    --chem-xlsx  "<path_to_chemical_input.xlsx>" \
    --cacao-xlsx "<path_to_cacao_uv_input.xlsx>" \
    --coffee-xlsx "<path_to_coffee_uv_input.xlsx>" \
    --out-xlsx   "outputs/SupplementaryData_S1.xlsx"

All three input spreadsheet paths are required explicitly.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import math
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
# Defaults (analysis defaults)
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


def control_event_tables(df: pd.DataFrame, tau_sens: Iterable[float]) -> pd.DataFrame:
    rows: List[dict] = []
    taus = [float(t) for t in tau_sens]
    ctrl = df[df["Treatment"].astype(str).str.upper() == "CONTROL"].copy()

    for iso, g in ctrl.groupby("Isolate", dropna=False):
        x = _to_num(g["LWR"])
        n = int(x.notna().sum())
        if n <= 0:
            continue

        for tau in taus:
            k = int((x >= tau).sum())
            p_hat = k / n
            lo, hi = clopper_pearson_ci(k, n, alpha=0.05)
            rows.append({
                "isolate": iso, "treatment": "Control", "tau": tau,
                "n": n, "k": k, "p_hat": p_hat, "ci_low": lo, "ci_high": hi,
                "definition": "LWR >= tau",
            })

    out = pd.DataFrame(rows)
    if len(out) > 0:
        out = out.sort_values(["tau", "isolate"]).reset_index(drop=True)
    return out


def _tau_label(tau: float) -> str:
    s = f"{float(tau):.2f}"
    return "tau" + s.replace(".", "_")


def make_tau_wide_table(df_long: pd.DataFrame) -> pd.DataFrame:
    if df_long is None or len(df_long) == 0:
        return pd.DataFrame()

    out_rows: List[dict] = []
    for keys, g in df_long.groupby(["isolate", "treatment"], dropna=False):
        iso, trt = keys
        row = {"isolate": iso, "treatment": trt}
        g = g.sort_values("tau")
        for _, r in g.iterrows():
            lbl = _tau_label(float(r["tau"]))
            row[f"n_{lbl}"] = int(r["n"])
            row[f"k_{lbl}"] = int(r["k"])
            row[f"p_hat_{lbl}"] = float(r["p_hat"])
            row[f"ci_low_{lbl}"] = float(r["ci_low"])
            row[f"ci_high_{lbl}"] = float(r["ci_high"])

        taus_sorted = sorted([float(x) for x in g["tau"].dropna().unique().tolist()])
        if len(taus_sorted) >= 2:
            t_first = _tau_label(taus_sorted[0])
            t_last = _tau_label(taus_sorted[-1])
            p_first = row.get(f"p_hat_{t_first}", np.nan)
            p_last = row.get(f"p_hat_{t_last}", np.nan)
            k_first = row.get(f"k_{t_first}", np.nan)
            k_last = row.get(f"k_{t_last}", np.nan)
            row[f"delta_p_hat_{t_first}_to_{t_last}"] = (
                float(p_last) - float(p_first)
                if pd.notna(p_first) and pd.notna(p_last) else np.nan
            )
            row[f"delta_k_{t_first}_to_{t_last}"] = (
                int(k_last) - int(k_first)
                if pd.notna(k_first) and pd.notna(k_last) else np.nan
            )
        out_rows.append(row)

    out = pd.DataFrame(out_rows)
    if len(out) > 0:
        out = out.sort_values(["isolate", "treatment"]).reset_index(drop=True)
    return out


def build_near_equal_potency_pairs(
    sm: pd.DataFrame,
    p_primary: pd.DataFrame,
    max_abs_d_lambda: float = 0.05,
    min_abs_d_I: float = 0.02,
    top_n: int = 30,
) -> pd.DataFrame:
    if sm is None or len(sm) == 0:
        return pd.DataFrame()

    rows: List[dict] = []
    p_lookup = None
    if p_primary is not None and len(p_primary) > 0:
        p_lookup = p_primary.set_index(["isolate", "treatment"])[["p_hat", "ci_low", "ci_high"]]

    for iso, sub in sm.groupby("isolate", dropna=False):
        by_trt = sub.set_index("treatment")
        trts = list(by_trt.index.tolist())
        for i in range(len(trts)):
            for j in range(i + 1, len(trts)):
                tA = trts[i]
                tB = trts[j]
                rA = by_trt.loc[tA]
                rB = by_trt.loc[tB]

                lamA = float(rA["lambda_area"])
                lamB = float(rB["lambda_area"])
                IA = float(rA["I_LWR"])
                IB = float(rB["I_LWR"])
                dlam = abs(lamA - lamB)
                dI = abs(IA - IB)

                if dlam <= max_abs_d_lambda and dI >= min_abs_d_I:
                    rec = {
                        "isolate": iso,
                        "treatment_A": tA,
                        "treatment_B": tB,
                        "lambda_A": lamA,
                        "lambda_B": lamB,
                        "abs_delta_lambda": dlam,
                        "I_A": IA,
                        "I_B": IB,
                        "abs_delta_I": dI,
                    }
                    if p_lookup is not None:
                        for suffix, trt in [("A", tA), ("B", tB)]:
                            try:
                                rec[f"p_hat_{suffix}_tau1_10"] = float(p_lookup.loc[(iso, trt), "p_hat"])
                                rec[f"ci_low_{suffix}_tau1_10"] = float(p_lookup.loc[(iso, trt), "ci_low"])
                                rec[f"ci_high_{suffix}_tau1_10"] = float(p_lookup.loc[(iso, trt), "ci_high"])
                            except Exception:
                                rec[f"p_hat_{suffix}_tau1_10"] = np.nan
                                rec[f"ci_low_{suffix}_tau1_10"] = np.nan
                                rec[f"ci_high_{suffix}_tau1_10"] = np.nan
                    rows.append(rec)

    out = pd.DataFrame(rows)
    if len(out) > 0:
        out = out.sort_values(["abs_delta_I", "abs_delta_lambda"], ascending=[False, True]).head(int(top_n)).reset_index(drop=True)
    return out


def axis_reconstruction_validation(
    df: pd.DataFrame,
    axis_recon: pd.DataFrame,
    how: str = "mean",
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    needed = {"Length", "Width", "Isolate", "Treatment"}
    if df is None or len(df) == 0 or not needed.issubset(set(df.columns)):
        return pd.DataFrame(), pd.DataFrame()

    direct_rows: List[dict] = []
    for iso, g_iso in df.groupby("Isolate", dropna=False):
        g_ctrl = g_iso[g_iso["Treatment"].astype(str).str.upper() == "CONTROL"]
        if len(g_ctrl) == 0:
            continue

        Lc = _agg(g_ctrl["Length"], how)
        Wc = _agg(g_ctrl["Width"], how)
        if pd.isna(Lc) or pd.isna(Wc) or float(Lc) == 0 or float(Wc) == 0:
            continue

        for trt, g_trt in g_iso.groupby("Treatment", dropna=False):
            if str(trt).strip().upper() == "CONTROL":
                continue

            Lt = _agg(g_trt["Length"], how)
            Wt = _agg(g_trt["Width"], how)
            if pd.isna(Lt) or pd.isna(Wt):
                continue

            direct_rows.append({
                "isolate": iso,
                "treatment": trt,
                "length_ratio_direct": float(Lt) / float(Lc),
                "width_ratio_direct": float(Wt) / float(Wc),
            })

    direct_df = pd.DataFrame(direct_rows)
    if len(direct_df) == 0 or axis_recon is None or len(axis_recon) == 0:
        return pd.DataFrame(), pd.DataFrame()

    merged = axis_recon.merge(direct_df, on=["isolate", "treatment"], how="left")
    if "a_ratio" not in merged.columns or "b_ratio" not in merged.columns:
        return pd.DataFrame(), pd.DataFrame()

    merged["a_abs_err"] = (pd.to_numeric(merged["a_ratio"], errors="coerce") - pd.to_numeric(merged["length_ratio_direct"], errors="coerce")).abs()
    merged["b_abs_err"] = (pd.to_numeric(merged["b_ratio"], errors="coerce") - pd.to_numeric(merged["width_ratio_direct"], errors="coerce")).abs()

    def _pearson_or_nan(x: pd.Series, y: pd.Series) -> float:
        x = pd.to_numeric(x, errors="coerce")
        y = pd.to_numeric(y, errors="coerce")
        ok = x.notna() & y.notna()
        if int(ok.sum()) < 2:
            return np.nan
        return float(np.corrcoef(x[ok].to_numpy(dtype=float), y[ok].to_numpy(dtype=float))[0, 1])

    summary = pd.DataFrame([
        {"metric": "n_pairs", "value": int(len(merged))},
        {"metric": "agg", "value": how},
        {"metric": "pearson_r_a_vs_length_ratio", "value": _pearson_or_nan(merged["a_ratio"], merged["length_ratio_direct"])},
        {"metric": "pearson_r_b_vs_width_ratio", "value": _pearson_or_nan(merged["b_ratio"], merged["width_ratio_direct"])},
        {"metric": "mean_abs_err_a", "value": float(pd.to_numeric(merged["a_abs_err"], errors="coerce").mean())},
        {"metric": "mean_abs_err_b", "value": float(pd.to_numeric(merged["b_abs_err"], errors="coerce").mean())},
        {"metric": "max_abs_err_a", "value": float(pd.to_numeric(merged["a_abs_err"], errors="coerce").max())},
        {"metric": "max_abs_err_b", "value": float(pd.to_numeric(merged["b_abs_err"], errors="coerce").max())},
    ])

    merged = merged.sort_values(["isolate", "treatment"]).reset_index(drop=True)
    return merged, summary


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
    p = argparse.ArgumentParser(
        prog="chem_analysis_main.py",
        description="Generate SupplementaryData_S1.xlsx from external input spreadsheets.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    p.add_argument("--chem-xlsx", type=Path, required=True,
                   help="Chemical morphology spreadsheet (required).")
    p.add_argument("--chem-sheet", type=str, default=None,
                   help="Optional: specify chemical sheet name explicitly. If omitted, auto-detect.")
    p.add_argument("--cacao-xlsx", type=Path, required=True,
                   help="UV-C cacao spreadsheet (required).")
    p.add_argument("--coffee-xlsx", type=Path, required=True,
                   help="UV-C coffee spreadsheet (required).")

    p.add_argument("--outdir", type=Path, default=(Path.cwd() / "outputs"),
                   help="Output directory (used if --out-xlsx is not provided). Defaults to ./outputs under the current working directory.")

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

    chem_path = Path(args.chem_xlsx)
    cacao_path = Path(args.cacao_xlsx)
    coffee_path = Path(args.coffee_xlsx)

    if not chem_path.exists():
        raise FileNotFoundError(f"Chemical spreadsheet not found: {chem_path}")
    if not args.skip_uv and not cacao_path.exists():
        raise FileNotFoundError(f"Cacao UV spreadsheet not found: {cacao_path}")
    if not args.skip_uv and not coffee_path.exists():
        raise FileNotFoundError(f"Coffee UV spreadsheet not found: {coffee_path}")

    out_xlsx: Path = (args.out_xlsx if args.out_xlsx is not None else (args.outdir / "SupplementaryData_S1.xlsx"))

    print(f"[INFO] Chemical input: {chem_path}")
    print(f"[INFO] Cacao UV input: {cacao_path}")
    print(f"[INFO] Coffee UV input: {coffee_path}")
    print(f"[INFO] Output workbook: {out_xlsx}")

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
    p_tau_wide = make_tau_wide_table(p_sens)
    p_ctrl_sens = control_event_tables(chem_df, tau_sens=tau_sens)
    p_ctrl_wide = make_tau_wide_table(p_ctrl_sens)

    print("[INFO] Building near-equal potency / divergent polarity support table...")
    supp_table_s1 = build_near_equal_potency_pairs(sm_mean, p_primary)

    print("[INFO] Validating axis reconstruction against direct Length/Width ratios...")
    axis_val, axis_val_summary = axis_reconstruction_validation(chem_df, ax_mean, how="mean")

    qa = qa_checks(sm_mean)

    cacao_df = None
    coffee_df = None
    lam_uv = pd.DataFrame()
    I_uv = pd.DataFrame()
    coffee_event = pd.DataFrame()

    if not args.skip_uv:
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

    with pd.ExcelWriter(str(out_xlsx), engine="openpyxl") as xw:
        sm_mean.to_excel(xw, index=False, sheet_name="summary_metrics")
        ax_mean.to_excel(xw, index=False, sheet_name="Axis_Reconstruction")
        p_primary.to_excel(xw, index=False, sheet_name="pI_Chem_tauPrimary")
        p_sens.to_excel(xw, index=False, sheet_name="pI_Chem_tauSensitivity")
        if len(p_tau_wide) > 0:
            p_tau_wide.to_excel(xw, index=False, sheet_name="pI_Chem_tauWide")
        if len(p_ctrl_sens) > 0:
            p_ctrl_sens.to_excel(xw, index=False, sheet_name="pI_Control_tauSensitivity")
        if len(p_ctrl_wide) > 0:
            p_ctrl_wide.to_excel(xw, index=False, sheet_name="pI_Control_tauWide")
        if len(supp_table_s1) > 0:
            supp_table_s1.to_excel(xw, index=False, sheet_name="Supplementary_Table_S1")
        if len(axis_val) > 0:
            axis_val.to_excel(xw, index=False, sheet_name="Axis_Reconstruction_Validation")
        if len(axis_val_summary) > 0:
            axis_val_summary.to_excel(xw, index=False, sheet_name="Axis_Recon_Val_Summary")

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
