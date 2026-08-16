import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from src.addressability import *
def test_constants():
    assert abs(comb_spacing_mhz()-3.31241721854)<1e-9
    assert abs(ceiling_mhz()-1.65620860927)<1e-9
    assert abs(dmin_mhz(1.35)-0.148)<0.002
def test_crossings_and_maxima():
    for b in crossing_fields_t(3): assert dmin_mhz(b)<1e-10
    for b in local_max_fields_t(4): assert abs(dmin_mhz(b)-ceiling_mhz())<1e-10
