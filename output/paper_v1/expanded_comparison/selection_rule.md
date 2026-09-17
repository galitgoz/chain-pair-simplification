# Frozen selection rule

Use supplied aligned w=16 protein pairs listed in decimated_pairs.csv and supplied
hurricane pairs listed in hurricanes/pairs.csv; retain the four original pairs.
No solver outcomes enter selection. Compute features from input coordinates only:
log(max(nA,nB)), log(max(nA,nB)/min(nA,nB)), mean chord/path-length ratio,
and centroid separation divided by mean path length. Standardize each feature
over the eligible candidate pool (zero-variance features contribute zero).
Select three additions per domain by farthest-first distance to the retained and
already selected pairs. First minimize reuse of individual structures/storms;
then maximize the minimum feature distance; break ties by pair identifier.
Protein identity for reuse is PDB structure, not chain: all supplied aligned
protein pairs share 1o7j.a, so that reuse is unavoidable. Only one additional
4eca chain will be selected when unused structures remain. All supplied protein
curves have 325 prepared points before the existing w=16 sampling, yielding 22
points; this is not full-resolution backbone data. New hurricanes use every
supplied row, with no alignment, cropping, or subsampling. Pair orientation is
the supplied orientation. Duplicate unordered pairs are removed deterministically.
Freeze all arrays, source hashes, numerical thresholds and configurations before
any new optimization. No selection by success, compression, or runtime.
