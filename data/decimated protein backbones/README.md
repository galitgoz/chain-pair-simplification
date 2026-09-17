# Decimated protein backbone dataset

Eight Cα backbone chains, each exactly 325 vertices, decimated at five stride
levels, in two coordinate frames. Built for chain pair simplification experiments
where δ3 is index-wise.

## Files

| file | rows | what it is |
|---|---|---|
| `decimated_curves.csv` | 10 144 | every vertex of every curve, one row per vertex |
| `decimated_chains.csv` | 80 | one row per (chain, frame, w): counts, edge statistics, gaps |
| `decimated_pairs.csv` | 70 | one row per (pair, frame, w): distances and thresholds |
| `load_decimated.py` | | loader returning numpy arrays; the decimator is embedded |

## The decimator

```python
def keep_indices(n, w):
    if w == 1:
        return np.arange(n)
    idx = list(range(0, n, w))
    if idx[-1] != n - 1:
        idx.append(n - 1)          # the last vertex is always kept
    return np.array(idx)
```

Three properties, all of them load-bearing:

**Uniform stride.** No Douglas-Peucker, no curvature sampling. A decimator that
picks out "important" vertices is itself a simplification, and it would take the
cheap compression before the algorithm under test ever runs.

**Lockstep.** The same index set is applied to every chain. δ3 is index-wise,
`max_t d(A'_t, B'_t) ≤ δ3`, so the index correspondence *is* the correspondence
between the chains. Decimating each chain on its own destroys it.

**Endpoints.** Row 0 and row 324 survive at every w, in every chain. Every
simplification keeps its endpoints, so the floor of δ3 must not move between
stride levels.

Verified on this build: all 80 (chain, frame, w) groups keep both endpoints, all
eight chains share one index set at every w, and the two frames retain identical
source rows.

## Coordinate frames

| frame | what it is |
|---|---|
| `deposited` | the coordinates exactly as in the PDB entry, no superposition |
| `aligned` | the same chains after a rigid Kabsch fit onto 1o7j.a |

Decimation is purely index-based, so the retained source rows are identical in
both frames and only the coordinates differ. Coordinates round-trip against the
source files to 0.0000 Å at 3 decimal places.

Use `aligned` for anything involving δ3. In the `deposited` frame the chains sit
wherever the crystallographers left them, so the pairwise distances there are
tens of ångström and carry no structural meaning. The frame is included so the
superposition can be redone (US-align, or any other transform) without going back
to the PDB: replace the coordinates, keep the index columns.

The Kabsch fit is index-wise least squares. It is meaningful for the homologues
and not for the three unrelated folds, where no rigid motion brings the chains
close. See "Known limits" below.

## Vertex counts and edge geometry

| w | vertices | mean edge (Å) | min edge | max edge | growth |
|---|---|---|---|---|---|
| 1 | 325 | 3.82 | 3.74 | 3.91 | 1.00 |
| 2 | 163 | 6.07 | 5.12 | 7.34 | 1.59 |
| 4 | 82 | 9.23 | 4.94 | 13.63 | 2.42 |
| 8 | 42 | 14.59 | 5.77 | 26.35 | 3.82 |
| 16 | 22 | 18.95 | 5.77 | 38.39 | 4.96 |

(1o7j.a, aligned frame; `decimated_chains.csv` has the figures for every chain.)

Two things to read off this table.

**The mean edge grows sublinearly.** From w=1 to w=16 the stride grows 16× but
the mean edge grows 5×, because the chain is folded and the chord is shorter than
the arc it replaces.

**Past roughly 10 Å the curve is no longer a protein backbone.** It is a coarse
trace that has lost the helix. `coarse_trace` flags it. w=8 and w=16 are over the
line: use them for runtime scaling only, and take headline results from w ≤ 4.

`appended_stub` marks the levels where the stride missed the last vertex and it
was appended, making the final index step shorter than w (w=8 and w=16 here, since
324 is divisible by 4 but not by 8). `last_stride` gives the actual final step.

## Thresholds in `decimated_pairs.csv`

```
delta3         = ceil(dF_discrete) on the curves ACTUALLY used, after decimation
delta1_alphaN  = N * mean edge after decimation, and delta2 = delta1
```

Both are recomputed at every stride level. α is held constant, not the absolute
tolerance in ångström, so these are bond-scaled experiments: a row at w=8 with
α=1 is not "the same tolerance" as a row at w=1 with α=1. The ceiling convention
follows Fan et al., arXiv:1409.2457, §7.

`feasibility_floor` is the exact lower bound on any δ3 that admits a solution: the
bottleneck of the best strictly increasing pairing from (0,0) to (n−1,m−1). δ3
does **not** have to reach the discrete Fréchet distance, because it constrains
only the retained pairs and a simplification may drop badly matched residues for
free. What it cannot drop is the endpoints, so `endpoint_floor` is a cheap lower
bound on the floor. Check `feasible` before blaming a solver.

## Known limits

**The three unrelated folds are endpoint-pinned.** For 1qd1.b, 1toh.a and 1d9q.d
the discrete Fréchet distance and the feasibility floor are equal and do not move
at all across w (35.98, 37.13, 47.96). Their maximum sits at an endpoint, which
decimation always keeps, so these instances are entirely determined by two
vertices. They are useful as a far-apart control and not much else.

**A smaller distance after decimation is not better agreement.** Decimation can
remove exactly the residues where the chains disagree, so `dF_discrete` falls
without the curves having moved closer. 1o7j.a/4eca.c drops from 6.98 at w=1 to
6.02 at w=2 for this reason, then rises again as the coarsening dominates.

**Gaps.** `n_gaps` and `gaps` list edges longer than 1.3× the median. At w=1 these
are crystallographic: 1toh.a has one bridging edge of 24.76 Å, 1d9q.d two of 15.16
and 14.23 Å, 4eca.b/c/d one each of 5.3 to 5.6 Å. The 4eca one is suspect: for
those entries the polymer sequence length equals the number of deposited Cα, so
there is no disorder in the structure, and the gap is a residue the extraction
dropped rather than one the crystal is missing. Resolve that before relying on
4eca gap statistics. At w > 1 the column no longer separates gaps from ordinary
long chords.

**This is the eight-chain set.** The seven additional chains (3ntx.a, 1wls.a,
2eq5.a, 2zsk.a, 1zq1.a, 3jq0.a, 2fep.a) are not here, and unlike these eight they
do not share a common length, so pairing them needs a stated policy for unequal
lengths before they can be decimated in lockstep.

## Usage

```python
from load_decimated import load, load_pair, thresholds, keep_indices

A, B = load_pair("4eca.c", w=4)             # (82, 3) each, same source rows
d1, d2, d3 = thresholds("4eca.c", w=4, alpha=1.0)
idxA, idxB = my_cps3f(A, B, d1, d2, d3)

rows = keep_indices(325, 4)                 # trace a vertex back to its residue
print(rows[idxA])                           # rows of the original 325-vertex chain
```

`python load_decimated.py` prints the instance table as a smoke test.

## Provenance

Extraction: ATOM records only, atom name `" CA "`, first model, altLoc blank or
`A`, coordinates copied verbatim from columns 31 to 54. Chain identifiers are the
author-assigned ones. Three identifiers in the source table are typographic
errors: `107j.a` is 1o7j.a, `4cea.b/.d` are 4eca.b/.d, and `1toh` is chain A. For
1hfj the author chain C is label chain B.

Checked against the pipeline in `full_run.csv`: 35 overlapping rows, identical
vertex counts and δ3, and discrete Fréchet agreeing to within the 0.005 Å
half-width of that file's two-decimal rounding.
