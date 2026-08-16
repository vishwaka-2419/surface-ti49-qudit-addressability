from pathlib import Path
import sys, numpy as np, pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.addressability import *
OUT=Path('outputs'); OUT.mkdir(exist_ok=True)
B=np.linspace(0.2,2.6,2401)
df=pd.DataFrame({'B_T':B,'Delta_inter_min_MHz':dmin_mhz(B)})
df.to_csv(OUT/'addressability_closed_form.csv',index=False)
print('S_MHz',comb_spacing_mhz())
print('period_T',period_t())
print('ceiling_MHz',ceiling_mhz())
print('crossings_T',crossing_fields_t())
print('local_maxima_T',local_max_fields_t())
print('dmin_1p35T_MHz',dmin_mhz(1.35))
