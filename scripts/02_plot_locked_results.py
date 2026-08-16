from pathlib import Path
import sys, numpy as np, pandas as pd, matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.addressability import *
D=Path('data/frozen'); O=Path('outputs'); O.mkdir(exist_ok=True)
fig,ax=plt.subplots(figsize=(6.4,4.0))
B=np.linspace(0.2,2.6,2000); ax.plot(B,dmin_mhz(B))
for b in crossing_fields_t(3): ax.axvline(b,ls=':',lw=.8)
for b in local_max_fields_t(4): ax.plot(b,dmin_mhz(b),'o')
ax.axvspan(.2,1.4,alpha=.12); ax.set(xlabel='B (T)',ylabel='leading-order minimum separation (MHz)')
fig.tight_layout(); fig.savefig(O/'addressability_comb.png',dpi=300); plt.close(fig)
# Hardware paired echo figure
p=pd.read_csv(D/'paired_seed_comparison.csv')
fig,ax=plt.subplots(figsize=(4.7,3.7))
for _,r in p.iterrows(): ax.plot([0,1],[r.binary,r.gray],marker='o',alpha=.7)
ax.set_xticks([0,1],['Binary','Gray']); ax.set_ylabel('matched inversion/echo return probability'); ax.set_ylim(0,1)
fig.tight_layout(); fig.savefig(O/'hardware_paired_echo.png',dpi=300); plt.close(fig)
print('wrote',O)
