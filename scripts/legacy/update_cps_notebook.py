"""Update the notebook builder's presentation and input contract."""

# Historical helper: run from the repository root.
import sys as _sys
from pathlib import Path as _Path
_sys.path[:0] = [str(_Path(__file__).resolve().parents[2] / 'src'), str(_Path(__file__).resolve().parents[2])]

from pathlib import Path

p=Path('src/build_cps_notebook.py')
s=p.read_text(encoding='utf-8')
start=s.index('**Scope:**')
end=s.index('\n\nSources:',start)
s=s[:start]+'''**Inputs:** by default, use every vertex of the supplied **aligned** decimated CSV curves, including their original endpoints and any gap-crossing edges. Tables 1–2 supply seven pairs. The seven Table 3 pairs are supplemented from the full cached audited chains, decimated separately with endpoints retained; unequal lengths are supported. Their cached rigid US-align transforms are held fixed. No historical result counts are treated as new measurements. Set `INPUT_SCOPE="qualified"` to use the previous continuous-interval preparation instead.'''+s[end:]
start=s.index('Use the project environment in')
end=s.index("''')",start)
s=s[:start]+'''Use the project environment in `requirements-analysis.txt`. The default loads the supplied aligned CSV coordinates without cropping or another decimation. Table 3 supplements use full cached residue coordinates and a previously computed rigid alignment. The coordinate-frame provenance is recorded per pair. The curves are mathematical polylines: edges crossing missing residues remain edges in this default experiment.'''+s[end:]
s=s.replace('INCLUDE_TABLE3 = True','INCLUDE_TABLE3 = True\nINPUT_SCOPE = "supplied"  # or "qualified": previously audited common continuous intervals')
s=s.replace('load_inputs(ROOT, WS, INCLUDE_TABLE3)','load_inputs(ROOT, WS, INCLUDE_TABLE3, scope=INPUT_SCOPE)')
s=s.replace("'source_first','source_last','full_chain_endpoints_preserved','TM_score','regime'", "'source_first','source_last','full_chain_endpoints_preserved','fragments_A','fragments_B'")
start=s.index('The existing preparation uses the first model;')
end=s.index("''')",start)
s=s[:start]+'''The residue maps preserve first-model Cα selections and document missing residues, alternate atom choices and fragment IDs. In `supplied` mode, gaps are reported but not removed; the exact input polyline is simplified. Source row 0 and the final source row are retained. Table 3 curves are sampled independently, so A and B can have different lengths; discrete Fréchet coupling does not require an index-wise correspondence. In `qualified` mode, the earlier common continuous intervals and their exclusions are used instead.'''+s[end:]
start=s.index("code('''from curve_algorithms import shortest_independent")
end=s.index("\nmd('''## 5.",start)
s=s[:start]+'''md("### 4.1 Global continuous Fréchet: free-space geometry and GCS\\n\\nThese executable definitions implement the global cost-layer sweep. Auxiliary functions include a continuous decision predicate and a discrete distance evaluator.")
geometry = (root/'src/curve_algorithms.py').read_text(encoding='utf-8')
geometry = geometry[:geometry.index('def cps_discrete_reference')]
code(geometry.replace('cache=True', 'cache=False'))
md("### 4.2 Joint solvers: discrete and continuous configuration graphs\\n\\nThe two graph builders share the two-count dynamic program. `certificate=False` forces graph execution; the default also accepts a feasible pair attaining the independent lower bound.")
joint = (root/'src/cps_paper_algorithms.py').read_text(encoding='utf-8')
joint = joint.replace('from curve_algorithms import (ball_segment_interval, segment_frechet_decision,\\n    continuous_decision, discrete_frechet, shortest_independent)', '# Geometry functions are defined in the preceding cell.')
code(joint.replace('cache=True', 'cache=False'))
'''+s[end:]
s=s.replace("from cps_notebook_run import run_comparison\nresults", "import cps_notebook_run as experiment\n# Run the solver definitions above; edits to those cells affect this experiment.\nexperiment.shortest_independent = shortest_independent\nexperiment.solve_cps = solve_cps\nexperiment.continuous_distance = continuous_distance\nexperiment.discrete_frechet = discrete_frechet\nexperiment.discrete_coupling = discrete_coupling\nrun_comparison = experiment.run_comparison\nresults")
s=s.replace('Common continuous interval qualification is a documented scope choice. It can remove full-chain endpoints and substantial regions; Table 3 supplements do not recover the unavailable historical extraction. A future full-chain study should analyze explicitly chosen fragment pairs or adopt a stated break-aware curve model.', 'The default retains all supplied vertices and treats gap-crossing edges as straight segments. The optional qualified mode removes regions outside common continuous intervals. Neither mode reconstructs missing residues. Table 3 supplements retain unequal full-chain lengths, but do not recover the unavailable historical extraction.')
s=s.replace('US-align is a structural alignment heuristic. Fixed input transforms, cropped intervals, present-day coordinates and decimation mean that the paper\'s historical δ₃ and counts need not match. The supplied Kabsch numbers are not silently reused.', 'The supplied pairs use their provided Kabsch frame; Table 3 supplements use cached US-align transforms fitted on qualified intervals and applied to full curves. These pair-specific frames, present-day coordinates and decimation mean that historical distances and counts need not match.')
s=s.replace("units='angstrom', endpoints=", "input_scope=INPUT_SCOPE, units='angstrom', endpoints=")
p.write_text(s,encoding='utf-8')
