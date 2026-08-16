"""Optional backend-targeted compilation only. This script does NOT submit a QPU job."""
import argparse, numpy as np
from scipy.linalg import expm
from qiskit import QuantumCircuit
from qiskit.circuit.library import UnitaryGate
from qiskit.transpiler import generate_preset_pass_manager
from qiskit_ibm_runtime import QiskitRuntimeService
parser=argparse.ArgumentParser(); parser.add_argument('--backend',default='ibm_marrakesh'); a=parser.parse_args()
# Keep repository release reproducible without credentials; uncomment only in an authenticated environment.
service=QiskitRuntimeService(); backend=service.backend(a.backend)
print('Backend target loaded:',backend.name,'- no job submitted')
