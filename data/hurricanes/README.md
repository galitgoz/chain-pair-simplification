# Curve data with measured packedness, for c-packed CPS experiments

Four families of polygonal curves, each curve annotated with its packedness under
both definitions in the literature, and each family supplying curve PAIRS with
their Frechet distance so the instances are usable for chain pair simplification
as they stand.

## 1. The two definitions are not the same number

**Ball definition** (Driemel, Har-Peled and Wenk): a curve P is c-packed if for
every point p and every r > 0, the length of P inside the Euclidean ball B(p, r)
is at most c r.

**Cube definition** (Gudmundsson, Sha and Wong, ISAAC 2020): the same with
S(p, r), the axis aligned cube centred at p whose *half* side length is r. Their
paper states the definition as "the length of the portion of pi contained in S is
at most c r, where r is half the side length of S", and chooses cubes because the
resulting function is piecewise linear, which balls are not.

Since B(p, r) is contained in S(p, r), which is contained in B(p, r sqrt(d)):

    c_ball <= c_cube <= sqrt(d) * c_ball

The paper states the right hand inequality in the same direction: a curve that is
c-packed under the ball definition is sqrt(d) c-packed under the cube definition.
In the plane the gap is a factor 1.414, in space 1.732. So a c value quoted from
that paper is a cube value and must not be substituted into a bound proved for
balls without the conversion, and the conversion is one sided: it bounds, it does
not translate.

The measurements here give both numbers for every curve, which removes the need
to convert at all. They also show the bracket is tight in practice, not a
pessimistic artefact:

| family | dim | max measured c_cube / c_ball | theoretical ceiling |
|---|---|---|---|
| hurricane tracks | 2 | 1.411 | 1.414 |
| character trajectories | 2 | 1.412 | 1.414 |
| protein backbones | 3 | 1.205 | 1.732 |

The planar families reach the ceiling almost exactly, which is expected: a curve
with a long straight stretch at 45 degrees to the axes realises the worst case,
since the cube then holds 2 sqrt2 r of it against 2r for the ball, a full extra
2(sqrt2 - 1)r of length. An axis parallel stretch gives no gap at all, so the
cube value depends on how the curve is oriented relative to the axes, and the
ball value does not.

## 2. What was measured, and what the numbers guarantee

`packedness.py` maximises f_p(r)/r, where f_p(r) is the exact length of the curve
inside the ball or cube of radius r centred at p, over centres at the vertices
of the curve and over candidate radii. Segment intersections are solved in
closed form, not sampled.

Restricting centres to vertices costs at most a factor 2 (Lemma 7 of the ISAAC
2020 paper, numbered as in arXiv:2009.07789v2). If an optimal cube S(q, r)
contains a vertex v, then S(v, 2r) contains S(q, r). If it contains no vertex,
the curve inside it consists of chords, the cube can be grown about q until its
boundary meets a vertex without lowering the ratio (Lemma 6), and the same
doubling applies. Both steps hold for balls as well, since the length of a chord
divided by the radius does not decrease as the ball grows.

Candidate radii are the distance to each vertex, the closest approach of each
segment, and a geometric grid of ratio 1.02 from 1/1000 of the extent of the
curve up to the extent. These are not all the breakpoints of f_p(r)/r: for cubes
the ratio can also turn where a corner of the cube crosses a segment (case (ii)
of Lemma 3 in the paper), and for balls the ratio is not monotone between these
radii at all. The grid covers both cases: since f_p is non decreasing, a grid
point just above the optimal radius loses at most 2 per cent. For curves with
more than 600 vertices the vertex and segment radii are skipped and only the
grid is used.

Every value is attained by an actual cube or ball, so it is a **certified lower
bound**, and the true packedness lies in [value, 2 * 1.02 * value]. On the
hurricane tracks the grid loss was checked directly. Recomputing c_cube with the
complete set of event radii, which gives exactly the maximum of the
2-approximation in Theorem 1 of the ISAAC 2020 paper, raised the values by 0.001
per cent at the median and 0.86 per cent at most (all 1295 tracks). Recomputing
c_ball with a grid of ratio 1.0005 changed the values by at most 0.006 per cent
(80 tracks).

The procedure differs from the paper's: the paper handles cubes only, sweeps the
event radii and updates the length incrementally in O(n^2 log n) time, while this
code recomputes the length at every candidate radius. The guarantee is the same
up to the grid factor, but the published values are for the authors' own
preprocessing of each data set, so they are comparable in kind, not value for
value (see the hurricane and character sections).

Validation on curves whose packedness is known in closed form:

| curve | measured cube | measured ball | exact |
|---|---|---|---|
| straight segment | 2.0000 | 2.0000 | 2 |
| zigzag of slope 1 | 2.8284 | 2.8284 | 2 sqrt2 = 2.8284 |
| circle, centres on the curve | 3.6685 | 3.1001 | pi = 3.1416 for balls |
| spiral, k turns | grows linearly in k | grows linearly in k | linear in k |

The circle is also the clean illustration of the factor 2: the true ball
packedness of a circle is 2 pi, attained by a ball centred at the centre, which
is not a point of the curve, and the on-curve optimum is exactly half of it. The
measured ball value is 1.3 per cent below pi because the test circle has 601
vertices, so only the grid is used, and the optimal radius r = 2 falls between
two grid points.

## 3. The data

### hurricanes/ — lowest packedness, 1295 curves, 2D

Atlantic and Northeast and North Central Pacific storm tracks from NOAA HURDAT2,
in the versions hurdat2-1851-2015-070616.txt (Atlantic, seasons 1851 to 2015) and
hurdat2-nepac-1949-2015-050916.txt (Pacific, seasons 1949 to 2015). A track is the
complete best track of one storm, with all its records in time order: the six
hourly positions, plus 577 records at other times in 336 of the tracks (531 of
them landfalls), and every stage of the system, including extratropical and
remnant low stages (12.6 per cent of the records). Tracks with fewer than 25
positions are dropped, as are tracks crossing the antimeridian, and repeated
consecutive positions are removed: 840 Atlantic and 455 Pacific tracks remain,
with 25 to 133 vertices (median 36).

Each basin is projected with its own azimuthal equidistant projection on a sphere
of radius 6371.0088 km, centred at 25N 60W for the Atlantic and 18N 120W for the
Pacific, so distances from the centre are exact and every track in a basin shares
one coordinate frame. Packedness is scale free but not projection free, so the
projection is part of the data definition. c_cube in addition depends on the
orientation of the axes, which here point east and north in the projection:
rotating 60 tracks by 15, 30 and 45 degrees changed c_cube by 10.6 per cent at
the median and by up to 39 per cent, while rotating them by 45 degrees changed
c_ball by at most 0.01 per cent.

    c_cube:  min 2.01   median 2.98   max 9.01
    c_ball:  min 2.01   median 2.41   max 8.00

These are as close to the theoretical floor as real curves get: any curve has
packedness at least 2, because a ball centred on it already holds about 2r of
length. This is the family to use when the point is to show that a c-packed
running time is an advantage.

The Hurdat row of the ISAAC 2020 paper describes a different population: Atlantic
only, and about 1788 tracks, close to every storm in the file rather than only
the long ones (the 2015 file has 1813 Atlantic storms, 1782 of them with at least
two records, and the paper's longest curve has 133 vertices, as does the longest
track in the file, AL031899). The paper does not state its coordinate system. Its
average of 3.24 is therefore not directly comparable; the mean c_cube of the 840
Atlantic tracks here is 3.57.

`pairs.csv` holds all 2222 pairs that pass a screen (same basin, both endpoints
within 500 km, length ratio between 0.7 and 1.43), with their continuous and
discrete Frechet distance. 277 of them have c_cube below 3 on both curves and
Frechet distance under 400 km; in those the shorter track is 1594 to 6894 km
long (median 3119). The closest of them is 164 km (Dalilia 1983 and Priscilla
1967, 5.7 per cent of the shorter track), the closest relative to length is 3.4
per cent (Charlie 1951 and an unnamed storm of 1903, 206 km), and the median is
10.8 per cent.

### characters/ — moderate packedness, 1429 curves, 2D

UCI Character Trajectories, 20 characters, pen tip velocities integrated into
pen paths. Only the training file of the dataset is used (trajectories_train.mat,
1429 curves), and coordinates are stored only for the 163 curves that appear in
pairs. Two samples of the same character are a natural chain pair.

    c_cube:  min 3.08   median 6.74   max 10.45
    c_ball:  min 2.63   median 5.96   max 10.09

For comparison, the ISAAC 2020 paper reports its "Pen" dataset at min 4.07,
average 8.79, max 20.82 under the cube definition. Same order, and the spread
overlaps; the difference can come from the data, since their Pen set has 2858
curves, twice as many as here, and is not necessarily reconstructed the same way.
Quote their numbers for their dataset and these for this one, not one for the
other.

`pairs.csv` holds 120 pairs, six per character: the six closest pairs under the
continuous Frechet distance among 14 samples of the character drawn at random
(seed 7), not among all samples.

### synthetic/ — packedness as a dial, 20 pairs, 2D

The serpentine family is the usable knob: k parallel passes across a unit square,
spaced 1/k apart. A ball of radius 1/2 in the middle meets about k passes, each
contributing about 2r, so c is about 2k. Measured:

    k        1     2     3     4     6     8    12    16    24    32
    c_cube  2.12  4.95  7.15  9.30 13.37 17.36 25.19 32.90 48.33 64.07
    c_ball  2.05  3.60  6.09  7.21 10.55 13.72 20.20 26.64 39.54 52.47

Each value is the larger of the two curves of the pair (`packedness_table.tex`
lists the base curve A instead). c_cube / 2k is 1.06 at k = 1 and 1.24 at k = 2,
then falls to 1.00 at k = 32, and is within 5 per cent of 2k from k = 12 on.
Spirals with t turns give c_cube from 3.9t at t = 1 down to 3.25t at t = 16.
Zigzags of slope s give exactly 2 sqrt(1 + s^2) for the base curve; the perturbed
copy is 10 to 30 per cent higher.

Every family member ships as a pair: a base curve and a smoothed random
perturbation of it. Along each ladder more than c changes: the number of vertices
grows with k (40 to 1311) and with t (49 to 769), and the perturbation is scaled
by 1/k and 1/t, so the Frechet distance between the two curves falls from 0.0167
to 0.0019 along the serpentine ladder. A running time against c plot on these
pairs does not hold n or the Frechet distance fixed.

### proteins/ and proteins_aligned/ — highest packedness, 8 curves, 3D

The alpha carbon backbones of the CPS-3F protein benchmark.

    c_cube:  49.28 to 50.77
    c_ball:  41.20 to 43.20

All eight sit at essentially the same value, which makes sense: each folds 1232
to 1256 A of backbone into a globule 59 to 70 A in diameter. In absolute terms c
is far higher than in the other families. Relative to size it is not: c/n is
0.155, the same ratio the ISAAC 2020 paper reports for its Hurdat set (0.154).

**Superposition matters, and the benchmark assumes it.** As deposited, the
chains sit in unrelated crystal frames, and the Frechet distance between
`1o7j.a` and the other chains is 49 to 84 A, which makes the delta values of the
published table impossible to satisfy. After a rigid Kabsch superposition
matching vertex i to vertex i:

| pair | RMSD | Frechet, aligned | delta3 in the table |
|---|---|---|---|
| 1o7j.a / 1hfj.c | 0.27 | 1.17 | 1 |
| 1o7j.a / 4eca.b | 9.89 | 6.89 | 5 |
| 1o7j.a / 4eca.d | 9.89 | 6.89 | 5 |
| 1o7j.a / 4eca.c | 9.88 | 7.02 | 6 |
| 1o7j.a / 1qd1.b | 20.45 | 36.03 | 21 |
| 1o7j.a / 1toh.a | 22.94 | 37.19 | 21 |
| 1o7j.a / 1d9q.d | 19.21 | 48.00 | 20 |

For the four similar pairs (1hfj.c and the three 4eca chains) the aligned
Frechet distance is within 2 A of the published delta3, which suggests the
benchmark is built on superimposed chains. For the three unrelated folds the
published delta3 (20 and 21) is below the distance between the first residues of
the superimposed chains (36.0, 37.1 and 48.0 A). With the endpoints fixed, as in
CPS, no pair of simplifications can meet it, so this superposition does not
reproduce those three rows of the table. `proteins_aligned/` holds the
superimposed curves.

## 4. Files

```
summary.csv                  one row per family, c ranges under both definitions
packedness_table.tex         LaTeX table of the low c curves and the serpentine ladder
hurricanes/<id>.txt          one track per file, x y in km
hurricanes/packedness.csv    per track: n, length, extent, c_cube, c_ball
hurricanes/pairs.csv         screened same basin pairs with Frechet distances
characters/<char>_<id>.txt   one pen path per file, only the 163 curves in pairs
characters/packedness.csv    per curve measurements
characters/pairs.csv         closest same character pairs
synthetic/<name>_A.txt, _B.txt   base curve and perturbed copy
synthetic/pairs.csv          per pair: c under both definitions, Frechet distance
synthetic/packedness.json    the same rows as synthetic/pairs.csv, as JSON
proteins/packedness.csv      the eight backbone curves
proteins/pairs.csv           as deposited Frechet distances
proteins_aligned/*.txt       Kabsch superimposed backbones
proteins_aligned/pairs_aligned.csv
code/packedness.py           the measurement, both definitions, with its own tests
code/frechet.py              discrete and continuous Frechet, with its own tests
code/build_*.py              the pipeline that produced each family
code/finalize.py             alignment, summary and LaTeX generation
figures/hurricane_c_witness.png   the square that gives c_cube, for three tracks
```

Curve files are whitespace separated coordinates, one vertex per line, in the
same format as the protein backbone curves.

## 5. Sources

NOAA National Hurricane Center, HURDAT2, files hurdat2-1851-2015-070616.txt and
hurdat2-nepac-1949-2015-050916.txt (https://www.nhc.noaa.gov/data/), cited as
Landsea and Franklin, "Atlantic Hurricane Database Uncertainty and Presentation
of a New Database Format", Monthly Weather Review 141(10):3576-3592, 2013,
doi:10.1175/MWR-D-12-00254.1. Used through the CSV conversion by Aleksey Bilogur,
`github.com/ResidentMario/hurdat2` (MIT License); rebuilding from its current
`atlantic_storms.csv` and `pacific_storms.csv` reproduces all 1295 tracks
(checked 15 September 2026). UCI Character Trajectories, via
`github.com/ShuningZhao/Character-Trajectories`. Protein backbones from the RCSB
PDB. Ball definition from Driemel, Har-Peled and Wenk, "Approximating the Frechet
Distance for Realistic Curves in Near Linear Time", Discrete & Computational
Geometry 48(1):94-127, 2012. Cube definition and the 2-approximation from Joachim
Gudmundsson, Yuan Sha and Sampson Wong, "Approximating the Packedness of
Polygonal Curves", ISAAC 2020, LIPIcs 181, 9:1-9:15,
doi:10.4230/LIPIcs.ISAAC.2020.9; journal version in Computational Geometry 108
(2023), article 101920, doi:10.1016/j.comgeo.2022.101920.
