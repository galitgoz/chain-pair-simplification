# Accuracy review — 2026-09-22

This review covers the saved primary comparison and the reorganized repository. It is not a proof of algorithm correctness or a full solver rerun.

## Numerical and input checks

- Confirmed 45 observed eligible pairs: 12 proteins and 33 hurricanes. Two further protein pairs are preparation-blocked, with eight method rows in the combined table.
- Confirmed 43 pairs with validated outputs for all four methods: 12 proteins and 31 hurricanes. The other two eligible hurricane pairs have CPS-2F transition-budget exits.
- Verified SHA-256 and vertex counts for all 45 saved input pairs.
- Independently recomputed all 135 thresholds from the saved coordinates: each fidelity threshold is half the corresponding mean edge length, and each coupling threshold is the input discrete Fréchet distance. Comparisons used relative tolerance `1e-12` and absolute tolerance `1e-9`.
- Cross-checked all 172 matched primary observations against the frozen table used by the supplementary plotting scripts: input/output sizes, thresholds and solver times agree.

Median solver times, in seconds, round to:

| Domain | CPS-2F | CPS-3F |
| --- | ---: | ---: |
| Proteins | 1.858 | 0.032 |
| Hurricanes | 7.442 | 0.039 |

These medians use the 43-pair common cohort, excluding resource exits. They describe recorded solver times, not end-to-end execution on a new machine. Run the numerical summary command in the [reproduction guide](REPRODUCING.md) to regenerate the paper counts and medians.

## Repository checks and limits

All tracked Python files parsed successfully, notebook JSON was readable, and local Markdown/HTML links resolved. The six saved-summary and repository-layout tests passed. The source modules were discoverable through an editable installation in a separate environment that reused existing system packages; this was not a clean installation of all historical dependencies.

The reorganization preserved core algorithm content, input data, recorded results, historical notebook outputs and the initial import manifest. The `v1.0.0` tag remains unchanged. Citation metadata explicitly refer to that tag; `main` contains later changes.

Full solver execution, all historical notebook workflows, external source availability, and installation of the entire pinned dependency set remain unverified. Some historical runners write fixed output paths, and some require Windows, US-align or unbundled research PDFs. See the [execution limitations](REPRODUCING.md#solver-execution-current-limitations).
