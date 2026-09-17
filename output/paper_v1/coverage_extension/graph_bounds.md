# Graph counting and numerical bounds

## Sources and assumptions
Gozoltzani thesis: Chapter 3, Observation 4, printed p.18 (PDF p.26);
Chapter 5, section 5.1, printed p.26 (PDF p.34). Manuscript:
sections 3.1, 4.1 and 4.4, PDF pp.7-8, 11-12, 15-16.
The experimental chapters are provisional and are not evidence for these results.
Both constructions use original output-vertex indices and matching positions on
the extended input curves. For a curve with n original vertices, each of n
spheres intersects each of n-1 nondegenerate straight edges in at most two
locations. Thus U(n)=2n(n-1) is a finite upper bound on added locations and
n*=n+a <= n+U(n). Assumption: no edge is contained in a sphere boundary.
The implementation skips degenerate edges and preserves their original endpoints,
counts distinct interior curve parameters (not distinct spatial coordinates),
and merges parameters within 1e-10. Endpoint intersections are not added.
The numerical tolerances do not create more than two candidates per edge/sphere.
No packedness assumption is made or tested here.

## Mapping to the implementation
For one curve, dense indices (i,p) cover S_A=n(n+aA) locations: owner i in
the extended curve and selectable original vertex p. CPS-3F has aA=0 and
S_A=n^2. These dense slots are NOT all valid single-curve graph vertices.
The code removes vertices not on an anchored source-to-terminal path; the
reported active count H_A is the remaining number. Source and terminal are
included. Every active vertex has one zero-cost self/wait transition.
E_A counts outgoing entries after this pruning, including these waits.
Auxiliary locations cannot be selected as simplified output vertices.

The finite input-only configuration-state bound is
B_input=n*m*(n+2n(n-1))*(m+2m(m-1)) for CPS-2F;
for CPS-3F it is n^2*m^2 (one bound, no duplicate auxiliary-based bound).
When completed single-curve graphs supply auxiliary counts,
B_extended=S_A*S_B is a tighter finite combinatorial bound. When active
counts are saved, B_active=H_A*H_B is tighter still. None enumerates all
coupling-valid states: the output coupling condition can only remove states.
The root and terminal require no additional vertices beyond these products.
The theory's O(n^3*m^3) counts potential CPS-2F configurations, not edges,
auxiliary locations or DP labels. It is displayed separately, never evaluated
with a hidden constant of one. CPS-3F's finite n^2*m^2 construction bound has
no auxiliary locations; auxiliary utilization 0/0 is not defined.

Discovered states are len(values), including the source. Dequeued states
count heap pops, including a terminal pop. completed_states_lower_bound counts
previously expanded nonterminal states; the last interrupted expansion may be
partial. Queue size is pending heap entries, not a progress percentage.
The product graph is implicit; all valid or reachable states are not enumerated
in advance. No claim that discovered states equal all valid states is made.

Candidate transition checks are counted only after excluding simultaneous waits
and testing destination coupling. A budget-triggering check is included even
when not relaxed. They are not reported as stored unique graph edges. A finite
bound on these checks, where both single graphs completed, is
E_A*E_B-H_A*H_B: each pair of outgoing single-graph entries defines one product
source and destination, the H_A*H_B simultaneous waits are excluded, and each
source is popped at most once because the queued set is never cleared.
Each single graph has unique (source,destination) entries by the construction's
index enumeration. Coupling and reachability further reduce checked candidates.
No inner-loop instrumentation was added to count stored product edges.
Single-graph edges including waits obey H*(H+1)/2 because every non-wait edge
strictly advances the dense topological index and occurs at most once.

DP labels X_v[r] are distinct from states: r counts advances in A's selected
index; at most n labels per state, each storing a best B count. Labels are not
instrumented. The objective reported here counts vertices including the first
one: max(r+1,z+1), not max(r,z). Unequal lengths and repeated coupling indices
are allowed. Feasibility validation alone does not prove optimality.

All percentages mean 100*observed/finite_bound, NOT search progress or a
prediction of remaining work. Interrupted counts are observations before
termination. Historical missing counters stay missing; no runs were repeated
to fill them. A missing percentage is explained by missing evidence, no graph
execution, or a zero/non-applicable denominator, never replaced by zero.
