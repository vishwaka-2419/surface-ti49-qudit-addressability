#!/usr/bin/env python
"""
QST extension: analytic B-|kappa| materials/control design map for 49Ti-like I=7/2 qudits.

Leading-order axial result:
    S = |kappa| / 14
    Delta_inter_min(B,kappa) = min_{ell=0,...,6} |ell*S - 2*gamma_N*B|
    Delta_max(kappa) = |kappa| / 28

The map is intentionally analytic. It is not the 16-level sensitivity envelope.
Use it to show:
- where inter-manifold crossings and maxima occur;
- how field selection moves operating windows;
- how |kappa| sets the hard addressability ceiling.

Outputs:
  design_map.csv
  fig_B_kappa_addressability.png
  fig_kappa_ceiling.png
  selected_cross_sections.csv
  design_map_summary.txt
"""

from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

MU_N_OVER_H_MHZ_T = 7.622593285
G_N_ABS_49TI = 0.314686
GAMMA_N = MU_N_OVER_H_MHZ_T * G_N_ABS_49TI  # MHz/T

KAPPA_CURRENT = 46.3738
KAPPA_3MHZ = 84.0
KAPPA_5MHZ = 140.0


def delta_leading(B_T, kappa_abs_MHz):
    B = np.asarray(B_T, dtype=float)
    K = np.asarray(kappa_abs_MHz, dtype=float)
    S = K / 14.0
    ell = np.arange(0, 7, dtype=float)
    # Broadcast to (..., ell)
    d = np.abs(S[..., None] * ell - (2.0 * GAMMA_N * B)[..., None])
    return np.min(d, axis=-1)


def ceiling(kappa_abs_MHz):
    return np.asarray(kappa_abs_MHz, dtype=float) / 28.0


def maxima_fields(kappa_abs_MHz, Bmin=0.0, Bmax=3.0):
    """
    Midpoints between adjacent crossing fields:
        B_cross(ell) = ell*S/(2 gamma)
        B_max(ell+1/2) = (ell+1/2)*S/(2 gamma)
    """
    S = float(kappa_abs_MHz) / 14.0
    rows = []
    for ell in range(0, 6):
        Bm = (ell + 0.5) * S / (2.0 * GAMMA_N)
        if Bmin <= Bm <= Bmax:
            rows.append((ell, Bm))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="analysis_v3/03_designmap")
    ap.add_argument("--Bmin", type=float, default=0.6)
    ap.add_argument("--Bmax", type=float, default=1.8)
    ap.add_argument("--dB", type=float, default=0.002)
    ap.add_argument("--Kmin", type=float, default=20.0)
    ap.add_argument("--Kmax", type=float, default=160.0)
    ap.add_argument("--dK", type=float, default=0.5)
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    B = np.arange(args.Bmin, args.Bmax + args.dB / 2, args.dB)
    K = np.arange(args.Kmin, args.Kmax + args.dK / 2, args.dK)
    BB, KK = np.meshgrid(B, K)
    DD = delta_leading(BB, KK)

    df = pd.DataFrame({
        "B_T": BB.ravel(),
        "kappa_abs_MHz": KK.ravel(),
        "delta_inter_min_leading_MHz": DD.ravel(),
        "ceiling_MHz": ceiling(KK).ravel(),
    })
    df.to_csv(out / "design_map.csv", index=False)

    # Main B-|kappa| map
    fig, ax = plt.subplots(figsize=(8.4, 5.8))
    im = ax.pcolormesh(B, K, DD, shading="auto")
    cb = fig.colorbar(im, ax=ax)
    cb.set_label(r"Leading-order $\Delta_{\mathrm{inter}}^{\min}$ (MHz)")

    levels = [0.5, 1.0, 2.0, 3.0, 5.0]
    cs = ax.contour(B, K, DD, levels=levels, linewidths=0.9)
    ax.clabel(cs, inline=True, fontsize=8, fmt=lambda x: f"{x:g} MHz")

    ax.axhline(KAPPA_CURRENT, linestyle="--", label=r"current $|\kappa_{49}|=46.37$ MHz")
    ax.axhline(KAPPA_3MHZ, linestyle=":", label=r"$|\kappa|=84$ MHz ($\Delta_{\max}=3$ MHz)")
    ax.axhline(KAPPA_5MHZ, linestyle="-.", label=r"$|\kappa|=140$ MHz ($\Delta_{\max}=5$ MHz)")
    ax.axvline(1.0292, linestyle="--", label="full-model robust field 1.0292 T")
    ax.axvline(1.35, linestyle=":", label="readout field 1.35 T")
    ax.axvline(1.714, linestyle="-.", label="historical high-field point 1.714 T")

    ax.set_xlabel("External field B (T)")
    ax.set_ylabel(r"$|\kappa|$ (MHz)")
    ax.set_title(r"Field–quadrupole design map for electron-conditioned $I=7/2$ nuclear addressability")
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(out / "fig_B_kappa_addressability.png", dpi=240)
    plt.close(fig)

    # Ceiling plot
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.plot(K, ceiling(K))
    ax.axvline(KAPPA_CURRENT, linestyle="--", label="current 49Ti scenario")
    ax.axhline(3.0, linestyle=":", label="3 MHz ceiling")
    ax.axhline(5.0, linestyle="-.", label="5 MHz ceiling")
    ax.set_xlabel(r"$|\kappa|$ (MHz)")
    ax.set_ylabel(r"Hard ceiling $\Delta_{\max}=|\kappa|/28$ (MHz)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "fig_kappa_ceiling.png", dpi=240)
    plt.close(fig)

    # Cross-sections at three useful kappa values.
    rows = []
    for kval, label in [
        (KAPPA_CURRENT, "current_49Ti"),
        (KAPPA_3MHZ, "3MHz_ceiling_target"),
        (KAPPA_5MHZ, "5MHz_ceiling_target"),
    ]:
        vals = delta_leading(B, np.full_like(B, kval))
        for b, d in zip(B, vals):
            rows.append({
                "label": label,
                "kappa_abs_MHz": kval,
                "B_T": b,
                "delta_inter_min_leading_MHz": d,
            })
    cross = pd.DataFrame(rows)
    cross.to_csv(out / "selected_cross_sections.csv", index=False)

    # Compact numerical summary.
    lines = []
    lines.append(f"gamma_N(49Ti) = {GAMMA_N:.9f} MHz/T")
    lines.append("")
    for kval, label in [
        (KAPPA_CURRENT, "current 49Ti"),
        (KAPPA_3MHZ, "3 MHz ceiling target"),
        (KAPPA_5MHZ, "5 MHz ceiling target"),
    ]:
        lines.append(f"{label}: |kappa| = {kval:.4f} MHz")
        lines.append(f"  hard ceiling = {ceiling(kval):.6f} MHz")
        mx = maxima_fields(kval, args.Bmin, args.Bmax)
        if mx:
            lines.append("  leading-order maxima in plotted field range:")
            for ell, bm in mx:
                lines.append(f"    ell={ell}+1/2 : B = {bm:.6f} T")
        else:
            lines.append("  no leading-order maxima in plotted field range")
        lines.append("")

    for b in [1.0292, 1.0357, 1.35, 1.714]:
        d = float(delta_leading(np.array([b]), np.array([KAPPA_CURRENT]))[0])
        lines.append(
            f"current |kappa| at B={b:.4f} T -> leading-order delta = {d:.6f} MHz"
        )

    (out / "design_map_summary.txt").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))
    print(f"\nOutputs written to: {out.resolve()}")


if __name__ == "__main__":
    main()
