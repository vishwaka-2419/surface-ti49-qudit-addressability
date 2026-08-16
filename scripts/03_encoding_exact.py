from pathlib import Path
import numpy as np
from scipy.linalg import expm
try:
    from qiskit import QuantumCircuit, transpile
    from qiskit.circuit.library import UnitaryGate
except ImportError as exc:
    raise SystemExit(
        "Qiskit is required for this script. Create the pinned environment with "
        "`conda env create -f environment.yml` and activate `ti49-qudit`."
    ) from exc

def spin_matrices(I):
    m=np.arange(I,-I-1,-1,dtype=float); d=len(m)
    Ip=np.zeros((d,d),complex)
    for j,mj in enumerate(m):
        mp=mj+1
        if mp<=I and mp in m:
            i=int(np.where(np.isclose(m,mp))[0][0]); Ip[i,j]=np.sqrt(I*(I+1)-mj*(mj+1))
    Im=Ip.conj().T; Iy=(Ip-Im)/(2j)
    return Iy
I=3.5; theta=np.pi/2; Iy=spin_matrices(I); U=expm(-1j*theta*Iy)
def gray(k): return k^(k>>1)
P=np.zeros((8,8));
for k in range(8): P[gray(k),k]=1
for name,V in [('binary',U),('gray',P@U@P.T)]:
    qc=QuantumCircuit(3); qc.append(UnitaryGate(V),range(3))
    tq=transpile(qc,basis_gates=['rz','sx','x','cx'],optimization_level=3,seed_transpiler=7)
    print(name,'cx',tq.count_ops().get('cx',0),'depth',tq.depth())
print('Frozen release gives seed-ensemble local medians 16 CX (binary) and 9 CX (Gray); exact compiler counts may change with software versions.')
