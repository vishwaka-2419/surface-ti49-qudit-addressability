# GitHub -> Zenodo DOI steps

1. Create a public GitHub repository named `surface-ti49-qudit-addressability` under your account and push the contents of this folder.
2. Log in to Zenodo and connect your GitHub account. In the Zenodo profile menu choose **GitHub**, sync repositories, and enable this repository.
3. In GitHub, create a release/tag `v1.0.0`. Zenodo will automatically ingest the enabled release and mint a version DOI.
4. Open the Zenodo record, copy the **version DOI** and concept DOI. Add the DOI badge to `README.md`.
5. Insert the version DOI in the QST Data and code availability statement and cite the software record in the bibliography if desired.
6. Do not change the frozen result files after release. For corrections, create `v1.0.1` or a later semantic version so the archival record remains auditable.

Official documentation:
- https://help.zenodo.org/docs/github/enable-repository/
- https://help.zenodo.org/docs/github/archive-software/github-upload/
- https://help.zenodo.org/docs/github/

Note: with GitHub integration, Zenodo mints the DOI when the release is ingested; a DOI cannot be pre-reserved through the GitHub integration itself.
