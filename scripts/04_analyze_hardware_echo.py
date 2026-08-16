from pathlib import Path
import pandas as pd
D=Path('data/frozen'); h=pd.read_csv(D/'hardware_summary.csv'); p=pd.read_csv(D/'paired_seed_comparison.csv')
print(h.to_string(index=False)); print('mean paired Gray-binary =',p.gray_minus_binary.mean())
print('gray wins all seeds =',bool((p.gray_minus_binary>0).all()))
