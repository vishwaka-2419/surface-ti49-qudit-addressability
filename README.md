# surface-ti49-qudit-addressability

Reproducibility release for **“Quadrupole-limited spectral addressability of an atomically addressable 49Ti nuclear-spin qudit.”**

This repository contains the manuscript-level analysis scripts, frozen derived result tables and plotting code used for the submission. It intentionally distinguishes source-backed measurements, isotope-scaled/model inputs, modeled control quantities and measured IBM echo-return data.

## What is reproducible here

- the closed-form 49Ti addressability comb, crossings, local maxima and quadrupole ceiling;
- publication plots from the locked derived tables;
- the exact local binary/Gray encoding comparison for the target spin-7/2 rotation;
- the paired statistical summary of the frozen IBM echo results;
- SHA-256 verification of the frozen release files.

The full historical simulation workspace contained additional exploratory and sensitivity scripts. This public release focuses on the frozen manuscript-level results and does not contain IBM credentials. The optional backend script performs compilation only unless explicitly enabled.

## Install

```bash
conda env create -f environment.yml
conda activate ti49-qudit
pytest -q
```

## Reproduce manuscript-level outputs

```bash
python scripts/01_addressability_comb.py
python scripts/02_plot_locked_results.py
python scripts/03_encoding_exact.py
python scripts/04_analyze_hardware_echo.py
python scripts/05_verify_frozen_release.py
```

## Zenodo DOI

The repository is prepared for GitHub-Zenodo integration using `CITATION.cff` and `.zenodo.json`. After pushing to GitHub, connect the repository in Zenodo and create a GitHub release tagged `v1.0.0`. Zenodo will ingest the release and mint a version DOI. See `ZENODO_RELEASE_STEPS.md`.

After Zenodo creates the DOI, add the DOI badge here and replace the DOI placeholder in the manuscript data-availability statement if desired.

## Citation

Please use the Zenodo version DOI for the exact release and cite the accompanying manuscript.

## License

Code: MIT (`LICENSE`). Derived data and original figures: CC BY 4.0 (`LICENSE_DATA`) unless a third-party source is explicitly identified.


## DOI links

- Zenodo GitHub integration: https://zenodo.org/account/settings/github/
- Official enable-repository guide: https://help.zenodo.org/docs/github/enable-repository/
- Official archive-a-GitHub-release guide: https://help.zenodo.org/docs/github/archive-software/github-upload/
