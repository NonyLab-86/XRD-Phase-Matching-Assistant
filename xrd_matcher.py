"""Core functions for XRD peak matching.

This module intentionally has no Streamlit dependency so it can be tested independently.
"""
from __future__ import annotations

import io
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.signal import find_peaks
from scipy.sparse.linalg import spsolve


def safe_stem(filename: str) -> str:
    """Return a filesystem-safe filename stem."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", Path(filename).stem).strip("._-") or "XRD"


def infer_output_name(experimental_names: list[str]) -> str:
    """Infer a stable output basename from uploaded experimental filenames."""
    stems = [safe_stem(name) for name in experimental_names]
    if not stems:
        return "XRD_phase_match"
    if len(stems) == 1:
        core = stems[0]
    else:
        prefix = os.path.commonprefix(stems).rstrip("._-")
        if len(prefix) >= 5:
            core = prefix
        else:
            core = "_".join(stems[:3])
            if len(stems) > 3:
                core += f"_plus{len(stems) - 3}"
    return f"{core}_XRD_phase_match"


def decode_text(data: bytes, suffix: str = "") -> str:
    """Decode text and perform a lightweight conversion of common RTF exports."""
    text = data.decode("utf-8", errors="ignore")
    if suffix.lower() == ".rtf" or text.lstrip().startswith(r"{\rtf"):
        text = re.sub(r"\\par[d]?\b", "\n", text)
        text = re.sub(r"\\tab\b", "\t", text)
        text = re.sub(r"\\'[0-9a-fA-F]{2}", " ", text)
        text = re.sub(r"\\[a-zA-Z]+-?\d* ?", " ", text)
        text = re.sub(r"[{}]", " ", text)
        text = re.sub(r"\\[{}\\]", " ", text)
    return text


def read_experimental_xy(data: bytes, filename: str) -> pd.DataFrame:
    """Read 2theta and intensity from a delimited text file or spreadsheet bytes."""
    suffix = Path(filename).suffix.lower()
    if suffix in (".xlsx", ".xls"):
        frame = pd.read_excel(io.BytesIO(data), header=None)
        frame = frame.apply(pd.to_numeric, errors="coerce")
        frame = frame.dropna(how="all").dropna(axis=1, how="all")
        if frame.shape[1] < 2:
            raise ValueError("Spreadsheet must contain two numeric columns: 2θ and intensity.")
        frame = frame.iloc[:, :2].copy()
        frame.columns = ["two_theta", "intensity"]
    else:
        rows = []
        for line in decode_text(data, suffix).splitlines():
            line = line.strip()
            if not line or line.startswith(("#", ";", "//")):
                continue
            parts = re.split(r"[,;\s]+", line)
            values = []
            for part in parts:
                try:
                    values.append(float(part))
                except ValueError:
                    continue
            if len(values) >= 2:
                rows.append(values[:2])
        if not rows:
            raise ValueError("No two-column numeric data could be read.")
        frame = pd.DataFrame(rows, columns=["two_theta", "intensity"])
    frame = frame.replace([np.inf, -np.inf], np.nan).dropna()
    frame = frame.sort_values("two_theta").drop_duplicates("two_theta").reset_index(drop=True)
    if len(frame) < 5:
        raise ValueError("Fewer than five valid data points were found.")
    return frame


def extract_reference_metadata(text: str, filename: str = "reference") -> dict:
    """Extract common metadata labels from a text/RTF reference export."""
    def find(patterns, default="Not stated"):
        for pattern in patterns:
            match = re.search(pattern, text, flags=re.I)
            if match:
                return re.sub(r"\s+", " ", match.group(1)).strip(" :;,.")
        return default

    formula = find([
        r"Empirical formula\s*:\s*([A-Za-z0-9₀-₉().+\- ]+)",
        r"Chemical formula\s*:\s*([A-Za-z0-9₀-₉().+\- ]+)",
        r"Formula\s*:\s*([A-Za-z0-9₀-₉().+\- ]+)",
    ])
    compound = find([
        r"Compound name\s*:\s*([^\n\r]+)",
        r"PDF index name\s*:\s*([^\n\r]+)",
        r"ICSD name\s*:\s*([^\n\r]+)",
    ], Path(filename).stem)
    crystal = find([r"Crystal system\s*:\s*([A-Za-z -]+)"])
    space_group = find([r"Space group\s*:\s*([A-Za-z0-9_/\- ]+)"])
    pdf = find([
        r"Reference code\s*:\s*([0-9]{2}[-–][0-9]{3}[-–][0-9]{4})",
        r"PDF(?:\s+index)?(?:\s+file)?(?:\s+number)?\s*(?:no\.?|number|code)?\s*:\s*([0-9]{2}[-–][0-9]{3}[-–][0-9]{4})",
        r"PDF\s*#?\s*([0-9]{2}[-–][0-9]{3}[-–][0-9]{4})",
    ])
    phase = find([
        r"Phase\s*/\s*classification\s*:\s*([^\n\r]+)",
        r"Phase\s*:\s*([^\n\r]+)",
        r"Structure type\s*:\s*([^\n\r]+)",
    ])
    if formula == "Not stated":
        formula = compound
    combined = f"{crystal} {space_group} {compound} {formula} {phase}".lower()
    if phase == "Not stated":
        if "laves" in combined:
            phase = "Laves phase"
        elif "bcc" in combined or "im-3m" in combined:
            phase = "BCC solid solution"
        elif "fcc" in combined or "fm-3m" in combined:
            phase = "FCC / face-centred cubic"
        elif "hcp" in combined or "p63/mmc" in combined or "p6_3/mmc" in combined:
            phase = "HCP / hexagonal"
        elif crystal.lower() == "hexagonal":
            phase = "Hexagonal phase (specific type not established)"
        elif crystal.lower() == "tetragonal":
            phase = "Tetragonal phase"
        elif crystal.lower() == "cubic":
            phase = "Cubic phase (specific type not established)"
        else:
            phase = "Unclassified reference phase"
    return {
        "chemical_formula_or_name": formula,
        "compound_name": compound,
        "crystal_system": crystal,
        "space_group": space_group,
        "phase_classification": phase,
        "pdf_reference_number": pdf,
    }


def parse_reference_peak_table(text: str) -> pd.DataFrame:
    """Parse common text exports with columns ending in 2theta and relative intensity.

    Expected common row form:
    peak_no h k l d_A 2theta_deg relative_intensity
    """
    rows = []
    for line in text.splitlines():
        if not re.match(r"^\s*\d+\s", line):
            continue
        nums = re.findall(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?", line)
        if len(nums) < 3:
            continue
        try:
            values = [float(value) for value in nums]
        except ValueError:
            continue
        theta, intensity = values[-2], values[-1]
        if 10 <= theta <= 160 and 0 <= intensity <= 10000:
            rows.append((theta, intensity))
    return (pd.DataFrame(rows, columns=["two_theta", "reference_intensity"])
            .drop_duplicates("two_theta").sort_values("two_theta").reset_index(drop=True))


def asymmetric_least_squares(y: np.ndarray, lam: float = 1e7,
                             p: float = 0.001, iterations: int = 20) -> np.ndarray:
    """Estimate an asymmetric least-squares baseline."""
    y = np.asarray(y, dtype=float)
    n = len(y)
    if n < 5:
        return np.zeros_like(y)
    difference = sparse.diags([1, -2, 1], [0, 1, 2], shape=(n - 2, n), format="csc")
    weights = np.ones(n)
    for _ in range(iterations):
        weight_matrix = sparse.spdiags(weights, 0, n, n)
        baseline = spsolve(weight_matrix + lam * (difference.T @ difference), weights * y)
        weights = p * (y > baseline) + (1 - p) * (y < baseline)
    return baseline


def detect_experimental_peaks(frame: pd.DataFrame, lower_limit: float = 20.0,
                              baseline_correction: bool = True,
                              prominence: float = 0.0,
                              min_separation_deg: float = 0.15) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Detect local maxima and return (peak table, processed intensity table)."""
    data = frame[frame["two_theta"] >= lower_limit].copy().reset_index(drop=True)
    if len(data) < 5:
        raise ValueError("Fewer than five data points remain after applying the 2θ lower limit.")
    x = data["two_theta"].to_numpy(dtype=float)
    raw = data["intensity"].to_numpy(dtype=float)
    baseline = asymmetric_least_squares(raw) if baseline_correction else np.zeros_like(raw)
    corrected = raw - baseline if baseline_correction else raw.copy()
    signal = np.clip(corrected, 0, None)
    step = float(np.median(np.diff(x)))
    distance = max(1, int(round(min_separation_deg / step))) if step > 0 else 1
    minimum_prominence = prominence if prominence > 0 else max(float(np.ptp(signal)) * 0.01, 1e-9)
    indices, properties = find_peaks(signal, prominence=minimum_prominence, distance=distance)
    peaks = pd.DataFrame({
        "peak_2theta": x[indices],
        "peak_intensity": corrected[indices],
        "peak_prominence": properties.get("prominences", np.zeros(len(indices))),
    }).sort_values("peak_2theta").reset_index(drop=True)
    processed = pd.DataFrame({
        "two_theta": x,
        "raw_intensity": raw,
        "estimated_baseline": baseline,
        "corrected_intensity": corrected,
    })
    return peaks, processed


def match_peaks(experimental_peaks: pd.DataFrame, reference_peaks: pd.DataFrame,
                tolerance: float = 0.60) -> tuple[int, pd.DataFrame]:
    """Count experimental peaks whose nearest reference peak is within tolerance.

    Each experimental peak is counted at most once for each reference pattern.
    """
    columns = ["experimental_2theta", "experimental_intensity", "nearest_reference_2theta",
               "reference_intensity", "delta_2theta", "matched"]
    if experimental_peaks.empty or reference_peaks.empty:
        return 0, pd.DataFrame(columns=columns)
    reference_x = reference_peaks["two_theta"].to_numpy(dtype=float)
    rows = []
    for _, peak in experimental_peaks.iterrows():
        experimental_x = float(peak["peak_2theta"])
        index = int(np.argmin(np.abs(reference_x - experimental_x)))
        reference_x_nearest = float(reference_x[index])
        delta = reference_x_nearest - experimental_x
        rows.append({
            "experimental_2theta": experimental_x,
            "experimental_intensity": float(peak.get("peak_intensity", np.nan)),
            "nearest_reference_2theta": reference_x_nearest,
            "reference_intensity": float(reference_peaks.iloc[index]["reference_intensity"]),
            "delta_2theta": delta,
            "matched": abs(delta) <= tolerance + 1e-12,
        })
    details = pd.DataFrame(rows)
    return int(details["matched"].sum()), details
