# Changelog

## Unreleased

- Added the missing SciPy dependency and pinned the supervisor's psutil dependency for isolated installations.

- Clarified the tagged citation target and how to avoid overwriting historical reports during checks.

- Moved 46 root-level Python files into `src/` and added editable-install metadata.
- Moved five earlier notebooks to `notebooks/legacy/`, three historical guides to `docs/legacy/`, and the initial import manifest to `docs/provenance/`.
- Updated imports, repository-relative paths, notebook setup and documentation links for the new layout.
- Preserved algorithm implementations, input data, recorded outputs and the `v1.0.0` citation tag.

## 1.0.0 — 2026-09-20

First tagged research snapshot for citation. This version packages the saved experiments and numerical reproduction entry point; the isolated-environment checks documented on `main` were performed after this tag.

- Added paper-focused documentation and a standard-library summary command.
- Checked the common cohort (43 pairs), paper counts and median solver times.
- Removed runtime font caches and visual QA preview PNGs; retained final figures and scientific evidence.
- Repaired historical documentation links.
- Moved 12 historical maintenance helpers into `scripts/legacy/` and adjusted their root-path handling and the notebook builder import.
- Added an MIT license for the original software and citation metadata.

Core optimization algorithms and recorded experimental results are unchanged. Third-party data and research publications retain their upstream rights. At this release, clean-environment installation, automated CI and portability had not been validated.
