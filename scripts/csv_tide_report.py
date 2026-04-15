#!/usr/bin/env python3
"""Generate tide time-series plots and harmonic-constant tables from CSV input.

Expected CSV columns (exact names):
- Timestamp
- TZ (optional)
- PRS1 (m)
- PRS2 (m)
- RAD1 (m)
- Solar (V) (optional)
- Battery (V) (optional)
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt

FOUR_CONSTITUENTS: Dict[str, float] = {
    "M2": 28.9841042,
    "S2": 30.0000000,
    "K1": 15.0410686,
    "O1": 13.9430356,
}

NINE_CONSTITUENTS: Dict[str, float] = {
    "M2": 28.9841042,
    "S2": 30.0000000,
    "K1": 15.0410686,
    "O1": 13.9430356,
    "N2": 28.4397295,
    "K2": 30.0821373,
    "P1": 14.9589314,
    "Q1": 13.3986609,
    "M4": 57.9682084,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot sea-level time domain and export harmonic constants (4 and 9 constituents)."
    )
    parser.add_argument("--input", required=True, help="Path to input CSV file.")
    parser.add_argument(
        "--column",
        default="PRS1 (m)",
        choices=["PRS1 (m)", "PRS2 (m)", "RAD1 (m)"],
        help="Water-level column to analyze.",
    )
    parser.add_argument("--output-dir", default="output", help="Output directory.")
    parser.add_argument(
        "--offset",
        type=float,
        default=0.0,
        help="Offset added after outlier correction (meters).",
    )
    parser.add_argument(
        "--hampel-window",
        type=int,
        default=7,
        help="Hampel filter window size for outlier detection.",
    )
    parser.add_argument(
        "--hampel-nsigma",
        type=float,
        default=3.0,
        help="Hampel threshold in sigma units.",
    )
    parser.add_argument(
        "--cutoff-cph",
        type=float,
        default=1 / 6,
        help="Low-pass cutoff in cycles per hour (default = 1/6 => 6-hour period).",
    )
    parser.add_argument(
        "--filter-order", type=int, default=4, help="Butterworth low-pass filter order."
    )
    return parser.parse_args()


def parse_timestamp_series(ts_series: pd.Series) -> pd.Series:
    """Parse Timestamp into timedelta.

    Supports:
    - HH:MM:SS(.f)
    - MM:SS(.f) -> interpreted as 00:MM:SS(.f)
    """
    cleaned = ts_series.astype(str).str.strip()
    mmss_mask = cleaned.str.match(r"^\d{1,2}:\d{2}(\.\d+)?$")
    cleaned.loc[mmss_mask] = "00:" + cleaned.loc[mmss_mask]
    return pd.to_timedelta(cleaned, errors="coerce")


def hampel_filter(x: np.ndarray, window_size: int, n_sigma: float) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(x, dtype=float)
    n = len(x)
    x_clean = x.copy()
    outlier_mask = np.zeros(n, dtype=bool)

    half = max(window_size // 2, 1)
    k = 1.4826  # MAD -> std factor

    for i in range(n):
        start = max(i - half, 0)
        stop = min(i + half + 1, n)
        window = x[start:stop]
        median = np.nanmedian(window)
        mad = np.nanmedian(np.abs(window - median))
        if not np.isfinite(mad) or mad == 0:
            continue
        threshold = n_sigma * k * mad
        if np.abs(x[i] - median) > threshold:
            outlier_mask[i] = True
            x_clean[i] = median

    return x_clean, outlier_mask


def lowpass_filter(x: np.ndarray, dt_hours: float, cutoff_cph: float, order: int) -> np.ndarray:
    fs = 1.0 / dt_hours
    nyquist = 0.5 * fs
    wn = cutoff_cph / nyquist
    if wn >= 1:
        raise ValueError(
            f"cutoff-cph={cutoff_cph} is too high for sample interval {dt_hours:.6f} h"
        )
    b, a = butter(order, wn, btype="low")
    return filtfilt(b, a, x)


def harmonic_fit(
    t_hours: np.ndarray,
    y: np.ndarray,
    constituents_degph: Dict[str, float],
) -> tuple[float, pd.DataFrame]:
    valid = np.isfinite(t_hours) & np.isfinite(y)
    t = np.asarray(t_hours[valid], dtype=float)
    yy = np.asarray(y[valid], dtype=float)

    names = list(constituents_degph.keys())
    freq_degph = np.array([constituents_degph[name] for name in names], dtype=float)
    omega = np.deg2rad(freq_degph)  # rad/hour

    cols = [np.ones_like(t)]
    for w in omega:
        cols.extend([np.cos(w * t), np.sin(w * t)])
    design = np.column_stack(cols)

    beta, *_ = np.linalg.lstsq(design, yy, rcond=None)
    mean_level = float(beta[0])

    rows = []
    idx = 1
    for i, name in enumerate(names):
        a_cos = beta[idx]
        b_sin = beta[idx + 1]
        idx += 2

        amplitude = float(np.hypot(a_cos, b_sin))
        phase_deg = float(np.degrees(np.arctan2(b_sin, a_cos)))
        rows.append(
            {
                "constituent": name,
                "freq_degph": float(freq_degph[i]),
                "amplitude_m": amplitude,
                "phase_deg": phase_deg,
            }
        )

    result = pd.DataFrame(rows).sort_values("amplitude_m", ascending=False).reset_index(drop=True)
    return mean_level, result


def main() -> None:
    args = parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input)
    required_cols = {"Timestamp", args.column}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    td = parse_timestamp_series(df["Timestamp"])
    y_raw = pd.to_numeric(df[args.column], errors="coerce")

    valid = td.notna() & y_raw.notna()
    t_hours = td[valid].dt.total_seconds().to_numpy() / 3600.0
    y_raw_np = y_raw[valid].to_numpy(dtype=float)

    if len(t_hours) < 10:
        raise ValueError("Not enough valid rows after parsing. Need at least 10 points.")

    dt_hours = float(np.median(np.diff(t_hours)))
    if not np.isfinite(dt_hours) or dt_hours <= 0:
        raise ValueError("Invalid time spacing. Ensure Timestamp is sorted and parseable.")

    y_no_outlier, outlier_mask = hampel_filter(
        y_raw_np, window_size=args.hampel_window, n_sigma=args.hampel_nsigma
    )
    y_offset = y_no_outlier + args.offset
    y_filtered = lowpass_filter(
        y_offset, dt_hours=dt_hours, cutoff_cph=args.cutoff_cph, order=args.filter_order
    )

    # Harmonic constants
    mean4, tbl4 = harmonic_fit(t_hours, y_filtered, FOUR_CONSTITUENTS)
    mean9, tbl9 = harmonic_fit(t_hours, y_filtered, NINE_CONSTITUENTS)

    summary = pd.DataFrame(
        [
            {"set": "4_constituents", "mean_level_m": mean4},
            {"set": "9_constituents", "mean_level_m": mean9},
        ]
    )

    tbl4_path = output_dir / "harmonic_constants_4.csv"
    tbl9_path = output_dir / "harmonic_constants_9.csv"
    summary_path = output_dir / "harmonic_summary.csv"

    tbl4.to_csv(tbl4_path, index=False)
    tbl9.to_csv(tbl9_path, index=False)
    summary.to_csv(summary_path, index=False)

    # Plot
    fig, ax = plt.subplots(figsize=(14, 7))
    ax.plot(t_hours, y_raw_np, label="Raw", linewidth=1.2, alpha=0.65)
    ax.scatter(
        t_hours[outlier_mask],
        y_raw_np[outlier_mask],
        c="red",
        s=20,
        label="Outlier",
        zorder=3,
    )
    ax.plot(t_hours, y_no_outlier, label="After outlier correction", linewidth=1.4)
    ax.plot(t_hours, y_filtered, label="Low-pass + offset", linewidth=2.0)

    ax.set_title(f"Sea-Level Time Series ({args.column})")
    ax.set_xlabel("Time (hours from start)")
    ax.set_ylabel("Water level (m)")
    ax.grid(True, alpha=0.3)
    ax.legend()

    png_path = output_dir / "tide_timeseries.png"
    fig.tight_layout()
    fig.savefig(png_path, dpi=150)
    plt.close(fig)

    print("Done. Generated files:")
    print(f"- {png_path}")
    print(f"- {tbl4_path}")
    print(f"- {tbl9_path}")
    print(f"- {summary_path}")


if __name__ == "__main__":
    main()
