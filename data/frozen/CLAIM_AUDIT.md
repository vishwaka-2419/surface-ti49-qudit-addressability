# FINAL CLAIM AUDIT

## Allowed central claims

1. The surface 49Ti nucleus has I=7/2 and therefore provides an eight-dimensional nuclear-spin Hilbert space.

2. The experimentally demonstrated seconds-scale quantity is a nuclear population lifetime, not a nuclear T2.

3. In the uncertainty/model sensitivity calculation, the 1.35 T readout-derived field has a worst sampled inter-manifold nuclear-transition separation of 0.0542 MHz.

4. The uncertainty-aware field search predicts a robust operating point near 1.714 T with worst sampled inter-manifold separation 1.093 MHz.

5. The proposed robust field is a model prediction and extrapolates beyond the approximately 1.4 T upper field of the spectroscopy dataset used for contextual calibration.

6. Simultaneous 14-tone control implements the same logical Iy(pi/2) operation in both static electron manifolds substantially more accurately near the predicted robust operating field than at 1.35 T over the studied control-rate range.

7. Nuclear T2 values obtained here are design requirements under a specified pure-dephasing model, not measured 49Ti coherence times.

8. Electron-relaxation calculations bound sensitivity to bath distinguishability. The phenomenological parameter eta is not experimentally measured.

9. Exact Gray encoding reduces the local three-qubit synthesis cost of the same d=8 Iy(pi/2) operation from 16 to 9 CX gates.

10. On a three-qubit linear coupling map, Gray encoding reduces median CX count by 59.1% relative to binary encoding for the same logical operation.

11. On ibm_marrakesh, the matched binary/Gray echo benchmark gives mean return probabilities 0.8490 and 0.9104, respectively.

12. Gray improves the matched measured return probability by 0.0614 on average across five paired transpiler seeds.

13. The IBM observable is inversion/echo return probability. It is not process fidelity, average gate fidelity, or a direct native-Ti-versus-IBM fidelity comparison.


## Claims NOT allowed

- Do not call the work the first nuclear qudit.
- Do not claim the predicted ~1.714 T operating field has been experimentally validated.
- Do not call the 5.3 s population lifetime a coherence time.
- Do not call 150 kHz a measured 49Ti Rabi frequency.
- Do not state that electron relaxation causes a known number of nuclear errors per gate.
- Do not describe eta as a measured physical parameter.
- Do not claim Gray encoding is globally optimal.
- Do not use the historical ~56-two-qubit-gate IBM number as the cost of a single Iy(pi/2) gate.
- Do not call IBM return probability a gate/process fidelity.
- Do not claim a direct fidelity advantage of Ti hardware over IBM hardware.
