# Prepared inputs

`inputs.json` indexes 45 binary NumPy input files. Each file contains arrays `A` and `B`, with the exact coordinates used in the recorded experiments. Hashes and original archive paths are preserved. The 12 protein pairs use three-dimensional coordinates in angstroms; 33 hurricane pairs use projected two-dimensional coordinates in kilometres. The common validated comparison includes 31 of those hurricane pairs.

Protein data originate from [RCSB PDB](https://www.rcsb.org/). `protein_sources.json` preserves source URLs, retrieval information and identifiers from the original preparation. The experiments use prepared C-alpha backbone curves, alignment and endpoint-preserving sampling; they share reference chain `1o7j.a`. Raw structures, residue mappings and detailed audits are retained in the [research archive](https://github.com/galitgoz/chain-pair-simplification/tree/v1.0.0/data).

Hurricane tracks originate from [NOAA HURDAT2](https://www.nhc.noaa.gov/data/): the [Atlantic 1851–2015 release dated July 6, 2016](https://www.aoml.noaa.gov/hrd/hurdat/hurdat2-1851-2015-070616.txt) and [Northeast/North Central Pacific 1949–2015 release dated May 9, 2016](https://www.aoml.noaa.gov/hrd/hurdat/hurdat2-nepac-1949-2015-050916.txt). The saved inputs are processed planar coordinates, not raw latitude/longitude. Filtering, duplicate removal and basin projection are documented in the [archived data guide](https://github.com/galitgoz/chain-pair-simplification/blob/v1.0.0/data/hurricanes/README.md).

Source data retain their upstream rights. This package supports rerunning from prepared inputs; rebuilding them requires the archived preparation workflow and its external dependencies.
