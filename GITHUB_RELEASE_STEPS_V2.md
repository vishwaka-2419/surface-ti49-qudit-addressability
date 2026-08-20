# GitHub and Zenodo release steps for v2.0.0

1. Open the existing repository `vishwaka-2419/surface-ti49-qudit-addressability`.
2. Upload the contents of this folder into the repository root, preserving the paths shown here. Existing root files `README.md`, `CITATION.cff`, `.zenodo.json` and `RELEASE_NOTES.md` should be replaced by the v2 versions.
3. Commit with a message such as: `Add v2.0.0 low-field control-window analysis`.
4. Verify that Zenodo GitHub integration remains enabled for the repository.
5. Create a GitHub release with tag `v2.0.0`.
6. Wait for Zenodo to archive the GitHub release and mint the new version DOI under the existing concept DOI `10.5281/zenodo.21969249`.
7. Update the manuscript/proof data-availability statement with the new version DOI.
8. Do not delete the v1.0.0 Zenodo record. Versioned releases should remain independently citable.
