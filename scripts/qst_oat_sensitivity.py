#!/usr/bin/env python
from pathlib import Path
import importlib.util
import pandas as pd

HERE = Path(__file__).resolve().parent
LOW = HERE / 'qst_lowfield_full16_scan.py'
spec = importlib.util.spec_from_file_location('lowfield', LOW)
low = importlib.util.module_from_spec(spec)
spec.loader.exec_module(low)

B = 1.0292
base = dict(kappa_MHz=-46.3738, gz=0.56, Btip_T=0.0679, phi_deg=5.0,
            Az_MHz=130.0, Aperp_over_Az=10.0/130.0, gperp_over_gz=1.0)
ranges = {
    'kappa_MHz': [-49.0292, -46.3738, -43.7184],
    'gz': [0.54, 0.56, 0.58],
    'Btip_T': [0.0658, 0.0679, 0.0700],
    'phi_deg': [3.0, 5.0, 7.0],
    'Az_MHz': [130.0, 132.1],
    'Aperp_over_Az': [0.0, 0.05, 0.10, 0.20],
    'gperp_over_gz': [0.5, 1.0, 1.5],
}
out = Path('data/v2/sensitivity')
out.mkdir(parents=True, exist_ok=True)
base_delta = low.delta_inter_min(B, base)
rows, summary = [], []
for par, vals in ranges.items():
    ds = []
    for value in vals:
        p = dict(base); p[par] = value
        d = low.delta_inter_min(B, p)
        ds.append(d)
        rows.append({'parameter': par, 'value': value, 'delta_MHz': d,
                     'change_from_source_MHz': d-base_delta})
    summary.append({'parameter': par, 'span_MHz': max(ds)-min(ds),
                    'min_MHz': min(ds), 'max_MHz': max(ds)})
pd.DataFrame(rows).to_csv(out/'oat_sensitivity_points.csv', index=False)
pd.DataFrame(summary).sort_values('span_MHz', ascending=False).to_csv(out/'oat_sensitivity_summary.csv', index=False)
print(pd.DataFrame(summary).sort_values('span_MHz', ascending=False).to_string(index=False))
