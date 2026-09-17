# Coverage extension: rule fixed before optimization

Alpha=0.5 only. Preserve the original ten pairs and their observations. Add the
nine requested paper protein identities; preparation eligibility is determined
from sources, not solver outcomes. Reuse the existing aligned w=16 arrays for
4eca.b and 4eca.c, including their documented omissions. Do not reconstruct them.
For new unequal-length inputs reuse small_pair_data.structural_fit (US-align,
rigid sequential structural superposition) on entire observed chains, then reuse
decimation_algorithms.decimate independently with w=16, retaining both endpoints.
No equal-length cropping, padding, guessed index correspondences or trimming to
the paper's length is allowed. Internal breaks in a newly prepared chain block
the pair: the existing residue preparation splits those chains, and using a
longest fragment or bridging gaps would introduce a different input policy.
Existing prepared gapped geometries remain explicitly identified as such.

Hurricanes: consider only unordered-unique pairs in the existing pairs.csv,
with finite 2D tracks, >=2 vertices, matching declared lengths, and matching
basin classifications (AL=Atlantic; EP/CP=Pacific). Retain every supplied row and the common basin frame.
Exclude the five existing pairs. Prefer unused storms. To cover ordinary sizes,
use the existing cohort's range of maximum track lengths when it contains enough
eligible pairs; split that range's candidate size ranks into five quantile bands.
Choose one pair from each band, in increasing band order. Within a band first
minimize reused storms, then maximize minimum standardized feature distance to
the original and already selected pairs, then use pair ID as deterministic tie
break. Features are log(max length), log(length ratio), mean chord/path ratio,
and centroid separation divided by mean path length. Existing feature() is reused.
No fidelity results, compression, runtime, or solver success enter selection.

Thresholds: delta1=0.5 mean_edge(A), delta2=0.5 mean_edge(B),
delta3=discrete_Frechet(A,B), full precision. All adopted budgets, thread settings,
fresh output-local caches, warm-up, validation and certificate=False are retained.
Run each eligible new slot once, sequentially. A resource exit is an observation.
Blocked inputs are neither observations nor infeasible instances. Readiness is
checked before each new configuration, outside measured solver phases.
