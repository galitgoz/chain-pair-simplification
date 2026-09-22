# Repository hygiene

The 2026-09-20 cleanup reviewed 8,038 tracked files and exact-content duplicate groups. It removed 40 generated files (3,655,602 bytes): four Matplotlib font caches and 36 raster previews used for visual inspection. Their exact paths and pre-removal hashes are in [the cleanup inventory](cleanup-2026-09-20.csv). Git history retains them if needed.

Final figures, input coordinates, cached PDB source files, raw solver results, validation records, timing/resource logs, frozen configuration files and historical evidence were retained. In particular, `data/pdb_cache` is source evidence, not disposable runtime cache. Identical residue mappings and reference-chain copies in different experiment records are intentional provenance copies; removing them would require changing consumers and archived references.

Generated font caches and QA preview PNGs are now ignored by Git. Visual-check scripts may regenerate them locally. Machine-readable QA reports remain versioned.

Eighteen pre-existing local links were repaired: historical notebook PDF links now lead to the research-source guide, and seven links in the archived expanded-comparison report now resolve to their original parent directory. Scientific content and numerical result records were not changed.

The original `docs/provenance/snapshot_manifest.json` describes the initial import and is deliberately not rewritten. It is not a current-file manifest.

The repository still includes historical experiment stages for traceability. This cleanup does not claim they are a minimal portable execution package; see the [repository map](REPOSITORY_MAP.md) and [reproduction limitations](REPRODUCING.md).
