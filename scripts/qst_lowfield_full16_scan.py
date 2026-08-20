#!/usr/bin/env python
"""
QST extension: full 16-level low-field robustness scan for surface 49Ti.

This script reconstructs the deterministic parameter/model sensitivity grid stated
in the final manuscript/SI and validates the implementation against the archived
1.35 T and 1.714 T envelope values before scanning the experimentally accessible
low-field maximum.

Outputs are written to --out (default: analysis_v3/01_lowfield).

Important:
- kappa_49 is an isotope-scaled same-EFG prediction, not a direct 49Ti measurement.
- 47Ti tip/g parameters are contextual sensitivity inputs.
- The output is a deterministic sensitivity envelope, not a statistical confidence interval.
"""

from __future__ import annotations
import argparse
import itertools
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Frequency-unit constants, MHz/T
MU_B_OVER_H = 13996.24555
MU_N_OVER_H = 7.622593285

I = 3.5
S = 0.5
G_N_ABS = 0.314686

LOCKED = {
    "B_1p35_T": 0.054196,
    "B_1p714_T": 1.092741,
}
VALIDATION_TOL_MHZ = 5e-5


def spin_matrices(j: float):
    """Jx, Jy, Jz in descending-m basis |j>, |j-1>, ..., |-j>."""
    m = np.arange(j, -j - 1, -1, dtype=float)
    n = len(m)
    jp = np.zeros((n, n), dtype=complex)
    for col, mi in enumerate(m):
        if col > 0:
            jp[col - 1, col] = np.sqrt(j * (j + 1) - mi * (mi + 1))
    jm = jp.conj().T
    jx = (jp + jm) / 2.0
    jy = (jp - jm) / (2.0j)
    jz = np.diag(m)
    return jx, jy, jz


sx, sy, sz = spin_matrices(S)
ix, iy, iz = spin_matrices(I)

E2 = np.eye(2)
E8 = np.eye(8)
Sx = np.kron(sx, E8)
Sy = np.kron(sy, E8)
Sz = np.kron(sz, E8)
Ix = np.kron(E2, ix)
Iy = np.kron(E2, iy)
Iz = np.kron(E2, iz)


def hamiltonian_mhz(B_T: float, p: dict) -> np.ndarray:
    """Full 16-level Hamiltonian H/h in MHz."""
    phi = np.deg2rad(p["phi_deg"])
    A_perp = p["Az_MHz"] * p["Aperp_over_Az"]
    g_perp = p["gz"] * p["gperp_over_gz"]

    H = MU_B_OVER_H * (
        p["gz"] * (B_T + p["Btip_T"] * np.cos(phi)) * Sz
        + g_perp * p["Btip_T"] * np.sin(phi) * Sx
    )
    H = H + MU_N_OVER_H * G_N_ABS * B_T * Iz
    H = H + p["Az_MHz"] * (Sz @ Iz)
    H = H + A_perp * (Sx @ Ix + Sy @ Iy)
    Hq = p["kappa_MHz"] / (2.0 * I * (2.0 * I - 1.0)) * (
        Iz @ Iz - 0.5 * (Ix @ Ix + Iy @ Iy)
    )
    return H + Hq


def nuclear_transition_sets(B_T: float, p: dict):
    """
    Seven adjacent nuclear-like transitions in each electron-like manifold.

    States are classified by sign(<Sz>) and ordered by <Iz>.
    """
    evals, evecs = np.linalg.eigh(hamiltonian_mhz(B_T, p))

    sz_exp = np.real(np.diag(evecs.conj().T @ Sz @ evecs))
    iz_exp = np.real(np.diag(evecs.conj().T @ Iz @ evecs))

    result = []
    for sign in (-1, +1):
        idx = np.where(sz_exp * sign > 0.0)[0]
        if len(idx) != 8:
            raise RuntimeError(
                f"Expected 8 states in manifold sign {sign}, found {len(idx)} at B={B_T}"
            )
        idx = idx[np.argsort(iz_exp[idx])]
        freqs = np.abs(np.diff(evals[idx]))
        result.append(freqs)
    return result[0], result[1]


def delta_inter_min(B_T: float, p: dict) -> float:
    fm, fp = nuclear_transition_sets(B_T, p)
    return float(np.min(np.abs(fm[:, None] - fp[None, :])))


def scenario_grid():
    # Explicit deterministic grid reconstructed from final manuscript/SI.
    values = {
        "kappa_MHz": [-49.0292, -46.3738, -43.7184],
        "gz": [0.54, 0.56, 0.58],
        "Btip_T": [0.0658, 0.0679, 0.0700],
        "phi_deg": [3.0, 5.0, 7.0],
        "Az_MHz": [130.0, 132.1],
        "Aperp_over_Az": [0.0, 0.05, 0.10, 0.20],
        "gperp_over_gz": [0.5, 1.0, 1.5],
    }
    keys = list(values)
    return [dict(zip(keys, vals)) for vals in itertools.product(*(values[k] for k in keys))]


def envelope_at_field(B_T: float, scenarios):
    vals = np.empty(len(scenarios), dtype=float)
    for i, p in enumerate(scenarios):
        vals[i] = delta_inter_min(B_T, p)
    iw = int(np.argmin(vals))
    return {
        "B_T": float(B_T),
        "worst_MHz": float(vals[iw]),
        "p05_MHz": float(np.quantile(vals, 0.05)),
        "median_MHz": float(np.median(vals)),
        "p95_MHz": float(np.quantile(vals, 0.95)),
        "best_MHz": float(np.max(vals)),
        "limiting_index": iw,
        "limiting_scenario": scenarios[iw],
    }


def scan(fields, scenarios, label="scan"):
    rows = []
    t0 = time.time()
    n = len(fields)
    for k, B in enumerate(fields, 1):
        r = envelope_at_field(float(B), scenarios)
        row = {key: val for key, val in r.items() if key != "limiting_scenario"}
        for key, val in r["limiting_scenario"].items():
            row[f"limit_{key}"] = val
        rows.append(row)
        if k == 1 or k == n or k % max(1, n // 10) == 0:
            print(
                f"[{label}] {k:4d}/{n}: B={B:.6f} T, "
                f"worst={r['worst_MHz']:.6f} MHz"
            )
    print(f"[{label}] elapsed: {time.time() - t0:.1f} s")
    return pd.DataFrame(rows)


def contiguous_window(df, threshold=1.0):
    """Window containing the global robust maximum with worst >= threshold."""
    j = int(df["worst_MHz"].to_numpy().argmax())
    good = df["worst_MHz"].to_numpy() >= threshold
    if not good[j]:
        return None
    lo = j
    hi = j
    while lo > 0 and good[lo - 1]:
        lo -= 1
    while hi + 1 < len(df) and good[hi + 1]:
        hi += 1
    return float(df.iloc[lo]["B_T"]), float(df.iloc[hi]["B_T"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--out", default="analysis_v3/01_lowfield",
        help="output directory"
    )
    ap.add_argument(
        "--quick", action="store_true",
        help="faster diagnostic grid before the publication-quality run"
    )
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    scenarios = scenario_grid()
    print(f"Scenario count: {len(scenarios)}")
    if len(scenarios) != 1944:
        raise RuntimeError("Scenario grid is not 1944 rows; stop.")

    # Validate against two locked publication-facing envelope values.
    validations = {}
    ok = True
    for name, B, target in [
        ("B_1p35_T", 1.35, LOCKED["B_1p35_T"]),
        ("B_1p714_T", 1.714, LOCKED["B_1p714_T"]),
    ]:
        r = envelope_at_field(B, scenarios)
        err = r["worst_MHz"] - target
        passed = abs(err) <= VALIDATION_TOL_MHZ
        ok = ok and passed
        validations[name] = {
            "B_T": B,
            "reconstructed_worst_MHz": r["worst_MHz"],
            "archived_worst_MHz": target,
            "difference_MHz": err,
            "tolerance_MHz": VALIDATION_TOL_MHZ,
            "PASS": passed,
            "limiting_scenario": r["limiting_scenario"],
        }
        print(
            f"VALIDATE {B:.3f} T: reconstructed={r['worst_MHz']:.9f} MHz, "
            f"archived={target:.9f}, diff={err:+.3e} -> "
            f"{'PASS' if passed else 'FAIL'}"
        )

    (out / "validation.json").write_text(
        json.dumps(validations, indent=2), encoding="utf-8"
    )
    if not ok:
        raise SystemExit(
            "VALIDATION FAIL. Do not use the low-field results until the mismatch is resolved."
        )

    # The analytically identified target is 1.0357 T. Search a wider local region.
    if args.quick:
        coarse = np.arange(0.95, 1.1100001, 0.002)
        refine_step = 0.0005
        refine_half_width = 0.010
    else:
        coarse = np.arange(0.95, 1.1100001, 0.001)
        refine_step = 0.0002
        refine_half_width = 0.012

    dfc = scan(coarse, scenarios, "coarse")
    dfc.to_csv(out / "lowfield_coarse.csv", index=False)

    jc = int(dfc["worst_MHz"].to_numpy().argmax())
    Bc = float(dfc.iloc[jc]["B_T"])
    refine = np.arange(
        Bc - refine_half_width,
        Bc + refine_half_width + refine_step / 2,
        refine_step
    )
    dfr = scan(refine, scenarios, "refine")
    dfr.to_csv(out / "lowfield_refined.csv", index=False)

    jr = int(dfr["worst_MHz"].to_numpy().argmax())
    best = dfr.iloc[jr].to_dict()

    analytic = envelope_at_field(1.0357, scenarios)
    window = contiguous_window(dfc, threshold=1.0)

    summary = {
        "scenario_count": len(scenarios),
        "analytic_candidate_B_T": 1.0357,
        "analytic_candidate_worst_full16_MHz": analytic["worst_MHz"],
        "robust_lowfield_B_T": float(best["B_T"]),
        "robust_lowfield_worst_MHz": float(best["worst_MHz"]),
        "robust_lowfield_median_MHz": float(best["median_MHz"]),
        "robust_lowfield_limiting_scenario": {
            k.replace("limit_", ""): float(v)
            for k, v in best.items()
            if k.startswith("limit_")
        },
        "coarse_window_worst_ge_1MHz_T": list(window) if window else None,
        "grid_note": (
            "Deterministic sensitivity envelope; not a probability distribution or confidence interval."
        ),
    }
    (out / "lowfield_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    # Plot coarse and refined envelope.
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(dfc["B_T"], dfc["worst_MHz"], label="worst sampled, coarse")
    ax.plot(dfr["B_T"], dfr["worst_MHz"], label="worst sampled, refined")
    ax.axvline(1.0357, linestyle="--", label="analytic 1.0357 T")
    ax.axvline(float(best["B_T"]), linestyle=":", label="full-ensemble robust point")
    ax.axhline(1.0, linestyle="--", label="1 MHz threshold")
    ax.set_xlabel("External field B (T)")
    ax.set_ylabel("Worst sampled inter-manifold separation (MHz)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "fig_lowfield_envelope.png", dpi=220)
    plt.close(fig)

    print("\n=== LOW-FIELD SUMMARY ===")
    print(json.dumps(summary, indent=2))
    print(f"\nOutputs: {out.resolve()}")


if __name__ == "__main__":
    main()
