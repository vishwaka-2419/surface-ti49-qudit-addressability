from pathlib import Path
import hashlib, pandas as pd
D=Path('data/frozen'); M=D/'PUBLIC_SHA256_MANIFEST.csv'
if not M.exists(): raise SystemExit('PUBLIC_SHA256_MANIFEST.csv missing')
df=pd.read_csv(M); bad=[]
for _,r in df.iterrows():
    p=D/r['file']; h=hashlib.sha256(p.read_bytes()).hexdigest() if p.exists() else 'MISSING'
    if h!=r['sha256']: bad.append((r['file'],h,r['sha256']))
print('PASS' if not bad else 'FAIL', 'files',len(df));
for x in bad: print(x)
raise SystemExit(1 if bad else 0)
