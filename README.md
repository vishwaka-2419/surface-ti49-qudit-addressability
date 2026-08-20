# surface-ti49-qudit-addressability

Reproducibility material for **“Quadrupole-limited control windows for an atomically addressable 49Ti nuclear-spin qudit.”**

Version 2.0.0 extends the original release with the experimentally accessible low-field robustness calculation, dense dual-manifold control thresholds, the field-quadrupole design map, and a one-at-a-time sensitivity decomposition.

## Main scientific outputs in v2.0.0

- closed-form rigid-comb addressability, including the explicit absolute-frequency/order-reversal derivation;
- general maximum-field rule for a hyperfine-coupled quadrupolar nucleus;
- validated 16-level sensitivity scan around the low-field maximum;
- minimax field `B = 1.0292 T`, with worst sampled separation `1.187772 MHz`;
- contiguous `>= 1 MHz` low-field interval `0.992-1.069 T`;
- static dual-manifold 99% first-crossing rates of `4.815 kHz` at 1.35 T and `21.9289 kHz` at 1.0292 T;
- one-at-a-time sensitivity showing dominant leverage from `kappa` and `A_perp/A_z`;
- field-quadrupole design map and hard ceiling `Delta_max = |kappa|/28` for 49Ti.

The control propagation conditions on a static electron manifold. The manuscript explicitly treats electron switching as a separate experimental requirement and does not interpret the static result as a complete driven-dissipative gate prediction for the present STM junction.

## Layout

- `scripts/` - calculation scripts and the historical 16-level/control model required by the revised analysis.
- `data/v2/lowfield/` - validation and low-field field-scan outputs.
- `data/v2/control/` - dense control sweeps and publication-facing threshold summary.
- `data/v2/sensitivity/` - one-at-a-time sensitivity tables.
- `data/v2/design/` - compact field-quadrupole design outputs.
- `figures/v2/` - final publication figures.
- `manuscript/` - revised manuscript and Supplementary Information source/PDF.

## Environment

The v2 extension requires Python 3 with NumPy, SciPy, pandas and Matplotlib. The original repository environment remains suitable for the manuscript-level scripts.

## Reproduction sequence

From the repository root:

```bash
python scripts/qst_lowfield_full16_scan.py --out data/v2/lowfield_recomputed
python scripts/qst_oat_sensitivity.py
python scripts/qst_B_kappa_design_map.py --out data/v2/design_recomputed
```

The time-dependent fourteen-tone propagation is substantially more expensive. The exact dense outputs used in the manuscript are supplied in `data/v2/control/`, and the corresponding source scripts are retained in `scripts/`.

## Versioning and Zenodo

The existing Zenodo concept DOI is `10.5281/zenodo.21969249`. After updating the GitHub repository, create a GitHub release tagged `v2.0.0` so Zenodo archives the revised release and mints a new version-specific DOI. Insert that version DOI into the final publication proof.

## Citation

Please cite the version-specific Zenodo release together with the accompanying manuscript.

## License

Code: MIT. Derived data and original figures: CC BY 4.0 unless a third-party source is explicitly identified.
