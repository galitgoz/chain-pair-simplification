"""Build the reviewable notebook; analysis implementations are embedded as cells."""
from pathlib import Path
import nbformat as nbf
import hashlib,json

ROOT=Path(__file__).resolve().parents[1]
cells=[]
def md(text): cells.append(nbf.v4.new_markdown_cell(text.strip()))
def code(text): cells.append(nbf.v4.new_code_cell(text.strip()))

md(r"""
# Protein backbone curves: independent global simplification

**Purpose.** Load the supplied backbone coordinates, reconcile the protein pairs in Fan et al., Section 7, Tables 1–3, prepare traceable continuous fragments, and compute a shortest **vertex-restricted, endpoint-preserving** simplification of each fragment under **global continuous Fréchet** fidelity.

The independent solver implements the free-space recurrence of **Algorithm 1 in the supplied full version of _Global Curve Simplification_, Section 3.3, printed page 13**. It stores integer cost layers rather than explicitly materializing elementary intervals. The conference version numbers its algorithms differently. A local shortcut/subcurve test is used only for an upper bound, never as the global feasibility criterion.

The notebook also provides an exact small-instance **discrete CPS-3F** reference and prepares all 21 paper experiments for later joint computation. **CPS-2F is reserved and unimplemented.** Protein-scale joint CPS results and reproduction of the paper's reported counts are not claimed.

**Units:** Å. **Indices:** zero-based curve and fragment indices; original PDB author and mmCIF label residue identifiers are retained separately. Run from this project directory with its Python environment. Dependencies are recorded in `requirements-analysis.txt`. The supplied notebook includes executed outputs; cached records allow subsequent offline use.

## 1. Load the data
""")
code("""
from pathlib import Path
import json, hashlib, time, itertools, inspect
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from IPython.display import display
from protein_data import prepare, read_raw, paper_experiments, SELECTIONS, atom_table

ROOT = Path.cwd()
assert (ROOT / 'data' / 'raw protein backbones').exists(), 'Run from the project root.'
OUTPUT = ROOT / 'output'
PREPARED = ROOT / 'data' / 'prepared'
FIGURES = OUTPUT / 'figures'
for folder in (OUTPUT, PREPARED, FIGURES): folder.mkdir(parents=True, exist_ok=True)
raw_files = sorted((ROOT / 'data' / 'raw protein backbones').glob('*.raw'))
raw_curves = {p.name.removesuffix('.raw'): read_raw(p) for p in raw_files}
sources = pd.DataFrame(json.loads((ROOT / 'data' / 'pdb_cache' / 'sources.json').read_text()))
experiments = paper_experiments()
POLICY = dict(model='lowest deposited model number per author chain',
    local_raw='preserve coordinates and exactly match deposited C-alpha atoms, including alternate locations',
    other_chains='C-alpha highest occupancy; ties prefer blank, then A, then lexical altloc',
    ca_break_A=4.5, peptide_break_A=2.0,
    breaks='split at absent selected label_seq_id, long C-alpha/C-N separation, or uncheckable C-N connectivity',
    endpoints='fixed within each independent fragment', alignment='none; deposited coordinates')
THRESHOLDS = (2.0, 4.0, 8.0)  # Exploratory sweep, distinct from the paper's CPS thresholds.
plt.rcParams.update({'font.size': 10, 'figure.dpi': 110, 'axes.spines.top': False, 'axes.spines.right': False})
""")
md('## 2. General properties\nThree brief cells describe the inputs before organizing the curves.')
code("""
raw_summary = pd.DataFrame([
    dict(chain=key, points=len(P), finite=bool(np.isfinite(P).all()),
         median_adjacent_A=np.median(np.linalg.norm(np.diff(P, axis=0), axis=1)),
         max_adjacent_A=np.max(np.linalg.norm(np.diff(P, axis=0), axis=1)))
    for key, P in raw_curves.items()])
display(raw_summary.round(3))
print('Coordinate-only files contain no residue, model, alternate-location, or chain-break annotations.')
""")
code("""
display(sources[['pdb_id', 'status', 'retrieved_utc']])
print('PDB source URLs and SHA-256 hashes are retained in data/pdb_cache/sources.json.')
print('107J is unavailable (HTTP 404); 4CEA exists and is checked separately from 4ECA.')
""")
code("""
display(experiments.groupby('table').agg(experiments=('B_printed','size'),
    min_delta1_A=('delta1','min'), max_delta1_A=('delta1','max'),
    min_delta3_A=('delta3','min'), max_delta3_A=('delta3','max')))
print(f'{len(raw_curves)} supplied raw chains; 15 selected chains; 21 experiments in Tables 1–3.')
print('Large raw coordinate jumps are diagnostics, not permission to interpolate.')
""")
md(r"""
## 3. Organize traceable curves and verify identifiers

`chains[key]` holds the ordered coordinates, a residue mapping table, the current deposited C-alpha coordinates, and a list of **separate** continuous fragments. `fragments[id]` is the direct input interface for the algorithms. `experiments` preserves the paper's printed labels and parameters.

**Selection policy.** The lowest numbered deposited model is selected; other available models are recorded and never concatenated. Raw coordinates are preserved and matched uniquely to deposited C-alpha atoms within 0.002 Å; every supplied point matched with zero error. For chains without raw files, choose one C-alpha per polymer residue by occupancy, with deterministic alternate-location tie breaking. Raw alternate choices are preserved even if different from that policy, and flagged. Nonstandard polymer residues are retained if a C-alpha is present; calcium ions and nonpolymer atoms are excluded using `label_seq_id`.

**Connectivity policy.** Follow polymer `label_seq_id`, preserving author numbering and insertion codes. Split at missing selected residues, Cα–Cα > 4.5 Å, C–N > 2.0 Å, or absent compatible C/N atoms. The latter is a conservative rule. Blank or unknown insertion-code sentinels become an empty insertion code in the curve mapping. Unresolved termini are recorded, not synthesized. No line, simplification edge, or distance calculation joins fragments. Current files are mmCIF, so author/label chain boundaries and peptide geometry are used rather than assuming PDB `TER` records exist.

The raw `4eca` files omit residue 12; this is a **raw-data omission**, not an unresolved deposited residue. The complete deposited representation remains available separately. No chain is trimmed to make its length agree with the paper.
""")
code("""
chains, audit, breaks, missing = prepare(ROOT, ca_break=POLICY['ca_break_A'], peptide_break=POLICY['peptide_break_A'])
fragments = {f['id']: f for chain in chains.values() for f in chain['fragments']}
display(audit[['printed','verified','reported_length','polymer_sequence_length','deposited_ca_length',
               'prepared_ca_length','label_chain','fragments']])
print(f'{len(fragments)} fragments; {len(breaks)} explicit interior breaks.')
# Verify the literal transposed identifier as well as its proposed interpretation.
literal_cif, literal_atoms = atom_table(ROOT/'data'/'pdb_cache'/'4cea.cif')
literal_ca = literal_atoms[(literal_atoms.label_atom_id=='CA') & (~literal_atoms.label_seq_id.isin(['.','?']))]
literal_counts = literal_ca.drop_duplicates(['auth_asym_id','label_seq_id']).groupby('auth_asym_id').size()
assert literal_counts['B'] == 150 and 'D' not in literal_counts.index
print('Literal 4CEA title:', ' '.join(literal_cif['_struct.title']))
print('Literal 4CEA polymer C-alpha counts by author chain:', literal_counts.to_dict())
""")
md(r"""
**Interpretations, not silent corrections.** `107j.a` is interpreted as **1O7J author chain A**: the literal identifier gives HTTP 404 and the local coordinates match 1O7J A. `4cea.b` and `4cea.d` are interpreted as **4ECA B/D**: the local coordinates match those chains, whereas 4CEA is HIV integrase, with 150 C-alpha residues in B and no protein D. These are evidence-supported typo interpretations, not confirmed author errata. `1toh` omits a chain; A is its only deposited protein chain and matches the raw file. `1hfj.c` means **author C, label B**.

The reported 297 for `1d9q.d` conflicts with 325 local and current deposited C-alpha atoms. Table 3 differs for 3NTX A (322 / 331), 1WLS A (316 / 328), 2ZSK A (219 / 226), 1ZQ1 A (418 / 437), and 3JQ0 A (457 / 478), where each pair is **reported / prepared**. The historical trimming and preprocessing remain unresolved. Sequence length, observed C-alpha count, and selected curve length are different quantities.
""")
code("""
display(audit[['verified','model','models_available','alternate_ca_residues',
               'selected_noncanonical_altloc','selected_zero_occupancy_ca','raw_max_coordinate_error_A']])
assert audit.selected_zero_occupancy_ca.sum() == 0, 'Zero-occupancy C-alpha requires an explicit new policy.'
display(breaks[['chain','before_index','after_index','before_residue','after_residue','missing_label_seq_ids','reason']])
display(missing.groupby(['chain','reason']).size().rename('residue_count').reset_index())
""")
code("""
sample = chains['1qd1.b']['mapping']
display(sample[['curve_index','fragment_index','author_chain','label_chain','label_seq_id',
                'author_seq_id','insertion_code','residue_name','model','altloc','occupancy']].head(8))
fragment_inventory = pd.DataFrame([dict(fragment=f['id'], n=len(f['xyz']),
    first_author_residue=f['residues'].iloc[0].author_seq_id,
    last_author_residue=f['residues'].iloc[-1].author_seq_id) for f in fragments.values()])
display(fragment_inventory)
for key, chain in chains.items():
    chain['mapping'].to_csv(PREPARED / f'{key}.residues.csv', index=False)
    for f in chain['fragments']:
        assert np.all(np.diff(f['residues'].label_seq_id) == 1)
        np.savez_compressed(PREPARED / (f['id'].replace('/', '.') + '.npz'),
            xyz=f['xyz'], chain_indices=f['chain_indices'])
audit.to_csv(OUTPUT / 'identifier_audit.csv', index=False)
breaks.to_csv(OUTPUT / 'chain_breaks.csv', index=False)
missing.to_csv(OUTPUT / 'missing_residues.csv', index=False)
experiments.to_csv(OUTPUT / 'paper_experiments.csv', index=False)
(PREPARED / 'policy.json').write_text(json.dumps(POLICY, indent=2))
""")
md(r"""
### Prepared fragments in 3D
All panels have the same two viewpoints, equal x/y/z aspect, and the **same numerical axis limits**. Each chain is translated to its own centroid **for display only**. No rotation, fitting, scaling, or coordinate change is applied to the analysis inputs. Fragment colors and endpoint markers expose breaks. These panels compare shape; they do not imply structural registration between proteins.
""")
code("""
from notebook_plots import plot_fragments, plot_simplified, plot_metrics
plot_fragments(chains, FIGURES)
""")
md(r"""
## 4. Continuous global Fréchet simplification

For each fragment $P=(p_0,\ldots,p_{n-1})$, minimize $|P'|$ over increasing vertex subsequences that include both endpoints, subject to $F(P,P')\leq\delta$. This is continuous, strong Fréchet distance in 3D. The curves are never resampled.

**Algorithm 1 correspondence.** All $u<v$ are shortcut-DAG edges. A spine stores the continuous free-space intervals between $P$ and a fixed retained vertex $p_v$; each shortcut strip compares **all of $P$** against segment $p_up_v$. The paper propagates the integer cost function $\varphi$ over elementary intervals and takes lower envelopes. Here, `reach[v,k,i]` stores the earliest reachable parameter on original edge $i$ at spine $v$ with **exactly $k$ links**. At a fixed cost and within one original edge, the reachable set is a suffix of its closed free-space interval: convexity lets the traversal pause at $p_v$ and advance to the interval's end. A minimum of lower endpoints therefore represents the union exactly. Taking the smallest feasible $k$ recovers the lower envelope's optimum. Parent states recover retained vertices.

Inside a cell, arrival from the bottom permits the whole right free interval; arrival only from the left restricts its lower endpoint. The analogous rule advances the top boundary. This is the strip reachability step of the paper's recurrence, with integer cost layers replacing explicit elementary-interval subdivision. **It is not a literal line-by-line implementation of the paper's interval data structure or its later improved algorithm.** With a feasible local bound $K\leq n-1$, time is $O(n^3K)$ and memory $O(n^2K)$; worst case $O(n^4)$ time and $O(n^3)$ memory.

The local shortcut solver supplies **only $K$**. It cannot exclude globally valid shortcuts. The decision routines use closed intervals, handle zero-length edges, and use a numerical predicate tolerance of `1e-10`. Mathematical optimality assumes exact predicates; these experiments use double precision. Reported continuous distances are bisection brackets of width at most $10^{-5}$ Å, subject to that predicate tolerance.
""")
algorithm_source=(ROOT/'src/curve_algorithms.py').read_text()
# Embed actual executable implementations in pedagogical groups.
chunks=[
 ('Geometric primitive and segment decision',algorithm_source[algorithm_source.index('import itertools'):algorithm_source.index('@njit(cache=True)\ndef continuous_decision')]),
 ('Full continuous and discrete distance routines',algorithm_source[algorithm_source.index('@njit(cache=True)\ndef continuous_decision'):algorithm_source.index('@njit(cache=True)\ndef local_shortcut_path')]),
 ('Feasible upper bound and global strip propagation',algorithm_source[algorithm_source.index('@njit(cache=True)\ndef local_shortcut_path'):algorithm_source.index('def shortest_independent')]),
 ('Shortest global simplification interface',algorithm_source[algorithm_source.index('def shortest_independent'):algorithm_source.index('@njit(cache=True)\ndef discrete_frechet')]),
 ('Discrete Fréchet dynamic program',algorithm_source[algorithm_source.index('@njit(cache=True)\ndef discrete_frechet'):algorithm_source.index('def cps_discrete_reference')])]
for title,body in chunks:
    md('### '+title)
    code(body.replace('cache=True','cache=False'))
md('### Verify global behavior before using protein results\nThe saved test report records 811 checks, including exhaustive vertex-subsequence enumeration on 160 random small curves. The cell below also exercises the notebook-defined solver against an exhaustive oracle and a global-versus-local counterexample.')
code("""
checks = json.loads((OUTPUT / 'algorithm_checks.json').read_text())
rng = np.random.default_rng(749)
for trial in range(20):
    P = np.cumsum(rng.normal(size=(7,3)),axis=0)
    delta = float(rng.uniform(.2,2.5))
    sizes = [len(ids) for k in range(2,len(P)+1)
        for middle in itertools.combinations(range(1,len(P)-1),k-2)
        for ids in [(0,)+middle+(len(P)-1,)] if continuous_decision(P,P[list(ids)],delta)]
    assert len(shortest_independent(P,delta)) == min(sizes)
example = checks['local_global_counterexample']
P = np.array(example['P']); d = example['delta']
assert len(shortest_independent(P,d)) == 3
assert len(local_shortcut_path(P,d)) == 4
line = np.array([[0.,0,0],[1.,0,0],[2.,0,0]])
assert continuous_decision(line,line[[0,2]],0)
assert discrete_frechet(line,line[[0,2]]) == 1
print(f"Saved suite: {checks['checks']} checks passed. Notebook: 20 exhaustive optimality checks and metric regressions passed.")
print('Global counterexample: 3 retained vertices; local-shortcut restriction: 4.')
""")
md(r"""
## 5. Run on the prepared fragments

The first cell runs the complete 325-point 1O7J A curve live at 4 Å. The following cell consumes the already executed 69-run sweep (23 fragments × 3 thresholds), verifies its data/implementation signature, and independently checks every saved simplification's fidelity. Set `RECOMPUTE_ALL = True` to rerun the sweep with the notebook-defined solver. This can take several minutes and uses Numba compilation on the first call. Cached timings are the original measured solver timings; they exclude data preparation and distance evaluation.
""")
code("""
live_fragment = fragments['1o7j.a/f00']
started = time.perf_counter()
live_ids = shortest_independent(live_fragment['xyz'], 4.0)
print(f'Live GCS: {len(live_fragment["xyz"])} → {len(live_ids)} vertices; {time.perf_counter()-started:.2f} s')
display(live_fragment['residues'].iloc[live_ids][['fragment_index','curve_index','author_seq_id','insertion_code']].head(10))
""")
code("""
RECOMPUTE_ALL = False
RESULT_DIR = OUTPUT / 'independent'
RESULT_DIR.mkdir(exist_ok=True)
def input_signature():
    inputs = sorted((ROOT/'data'/'pdb_cache').glob('*.cif')) + raw_files
    inputs += [ROOT/'src/protein_data.py', ROOT/'src/curve_algorithms.py']
    payload = {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    payload['thresholds'] = list(THRESHOLDS)
    payload['policy'] = POLICY
    return hashlib.sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest()

signature = input_signature()
manifest_path = RESULT_DIR/'manifest.json'
manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
cache_valid = (manifest.get('signature')==signature and (RESULT_DIR/'results.csv').exists()
               and (RESULT_DIR/'paths.json').exists())
if RECOMPUTE_ALL or not cache_valid:
    result_rows=[]; paths={}
    for f in fragments.values():
        P=f['xyz']
        for delta in THRESHOLDS:
            started=time.perf_counter(); ids=shortest_independent(P,delta); seconds=time.perf_counter()-started
            lower,upper=continuous_distance(P,P[ids])
            result_rows.append(dict(chain=f['chain'],fragment=f['id'],delta_A=delta,n=len(P),k=len(ids),
                retained_fraction=len(ids)/len(P),continuous_F_lower_A=lower,continuous_F_upper_A=upper,
                discrete_F_A=float(discrete_frechet(P,P[ids])),runtime_s=seconds,local_k=len(local_shortcut_path(P,delta))))
            paths[f'{f["id"]}@{delta:g}']=dict(fragment_indices=ids.tolist(),chain_indices=f['chain_indices'][ids].tolist())
    results=pd.DataFrame(result_rows)
    results.to_csv(RESULT_DIR/'results.csv',index=False)
    (RESULT_DIR/'paths.json').write_text(json.dumps(paths,indent=2))
    manifest_path.write_text(json.dumps(dict(signature=signature,mode='notebook recomputed'),indent=2))
else:
    results=pd.read_csv(RESULT_DIR/'results.csv')
    paths=json.loads((RESULT_DIR/'paths.json').read_text())
    print('Using the verified, previously executed full sweep. Set RECOMPUTE_ALL=True after editing solver cells.')

assert len(results)==len(fragments)*len(THRESHOLDS)
retained_residues=[]
for row in results.itertuples():
    f=fragments[row.fragment]; P=f['xyz']; record=paths[f'{row.fragment}@{row.delta_A:g}']
    ids=np.array(record['fragment_indices'])
    assert ids[0]==0 and ids[-1]==len(P)-1 and np.all(np.diff(ids)>0)
    assert len(ids)==row.k and np.array_equal(f['chain_indices'][ids],record['chain_indices'])
    assert continuous_decision(P,P[ids],row.delta_A)
    assert row.continuous_F_upper_A<=row.delta_A+1e-5
    assert row.k<=row.local_k
    chosen=f['residues'].iloc[ids].copy(); chosen['delta_A']=row.delta_A
    chosen['fragment']=row.fragment; chosen['simplification_index']=np.arange(len(ids))
    retained_residues.append(chosen)
for _, group in results.groupby('fragment'):
    assert np.all(np.diff(group.sort_values('delta_A').k)<=0)
pd.concat(retained_residues).to_csv(RESULT_DIR/'retained_residue_mapping.csv',index=False)
assert len(live_ids)==int(results[(results.fragment=='1o7j.a/f00') & (results.delta_A==4)].k.iloc[0])
display(results.groupby('delta_A').agg(fragments=('fragment','size'),input_vertices=('n','sum'),
    retained_vertices=('k','sum'),solver_seconds=('runtime_s','sum')).round(2))
""")
md(r"""
## 6. Reveal simplification, fidelity, and cost
The overlays retain the original sampling as a thin gray trace and the shortest simplification as an orange polyline. Each fragment is simplified separately. Summed vertex counts below describe the **fragment collection**, not a single artificial chain. The diagonal in the fidelity panel is equality; discrete fidelity can exceed the continuous threshold even for a valid continuous simplification.
""")
code("""
plot_simplified(chains, paths, delta=4.0, directory=FIGURES)
""")
code("""
chain_summary = plot_metrics(results, FIGURES)
chain_summary.to_csv(RESULT_DIR/'chain_totals.csv', index=False)
display(chain_summary[chain_summary.delta_A==4][['chain','n','k']].reset_index(drop=True))
""")
md(r"""
## 7. Prepare the discrete CPS comparison

The paper's discrete CPS-3F objective is to minimize $\max\{|A'|,|B'|\}$ with ordered vertex selections, subject to $d_{dF}(A,A')\leq\delta_1$, $d_{dF}(B,B')\leq\delta_2$, and $d_{dF}(A',B')\leq\delta_3$. It does **not** require selected endpoints to equal original endpoints; its initial/final configurations permit other selected vertices. The independent GCS formulation here **does** preserve endpoints. Fair comparisons must declare that difference or run an additional endpoint-matched variant.

The exact exhaustive function below is a **small-instance reference**, not the scalable dynamic program (Algorithm 2) in Fan et al. It enumerates all feasible subsequences and finds the smallest maximum size. A guard prevents accidental exponential runs on full proteins. The demo uses contiguous eight-residue windows from a paper pair and the Table 1 tolerances; it is a software demonstration, not reproduction of a paper row.
""")
code(algorithm_source[algorithm_source.index('def cps_discrete_reference'):])
code("""
window_A=fragments['1o7j.a/f00']['xyz'][:8]
window_B=fragments['1hfj.c/f00']['xyz'][:8]
cps_demo=cps_discrete_reference(window_A,window_B,4,4,1,fix_endpoints=False)
cps_fixed=cps_discrete_reference(window_A,window_B,4,4,1,fix_endpoints=True)
demo_rows=[]
for name,solution in [('Discrete CPS, paper endpoints',cps_demo),('Discrete CPS, fixed endpoints',cps_fixed)]:
    if solution is None:
        demo_rows.append(dict(method=name,status='infeasible in deposited frame'))
    else:
        ia,ib=solution['indices_A'],solution['indices_B']
        demo_rows.append(dict(method=name,status='optimal on eight-residue windows',k=solution['k'],
            fidelity_A=discrete_frechet(window_A,window_A[ia]),fidelity_B=discrete_frechet(window_B,window_B[ib]),
            pair_discrete_F=discrete_frechet(window_A[ia],window_B[ib])))
ia=shortest_independent(window_A,4); ib=shortest_independent(window_B,4)
demo_rows.append(dict(method='Independent continuous GCS, fixed endpoints',status='different fidelity constraint',
    k=max(len(ia),len(ib)),fidelity_A=continuous_distance(window_A,window_A[ia])[1],
    fidelity_B=continuous_distance(window_B,window_B[ib])[1],pair_discrete_F=discrete_frechet(window_A[ia],window_B[ib])))
display(pd.DataFrame(demo_rows).round(4))
print('A window: author residues',fragments['1o7j.a/f00']['residues'].iloc[:8].author_seq_id.tolist())
print('B window: author residues',fragments['1hfj.c/f00']['residues'].iloc[:8].author_seq_id.tolist())
""")
code("""
by_printed=audit.set_index('printed')
pair_jobs=[]
for row in experiments.to_dict('records'):
    a=by_printed.loc[row['A_printed'],'verified']; b=by_printed.loc[row['B_printed'],'verified']
    pair_jobs.append(dict(**row,A_verified=a,B_verified=b,
        A_fragments=[f['id'] for f in chains[a]['fragments']],B_fragments=[f['id'] for f in chains[b]['fragments']],
        coordinate_frame='deposited; no fitted transform',
        status='pending registration choice and fragment pairing; full CPS solver not implemented'))
(OUTPUT/'pair_jobs.json').write_text(json.dumps(pair_jobs,indent=2))
display(pd.DataFrame(pair_jobs)[['table','A_printed','A_verified','B_printed','B_verified',
    'delta1','delta2','delta3','paper_CPS3F']])
# Whole-pair distances are computed only when BOTH chains contain one fragment.
pair_diagnostics=[]
for b in audit.verified:
    if b=='1o7j.a' or len(chains[b]['fragments'])!=1: continue
    distance=float(discrete_frechet(chains['1o7j.a']['fragments'][0]['xyz'],chains[b]['fragments'][0]['xyz']))
    pair_diagnostics.append(dict(A='1o7j.a',B=b,deposited_frame_discrete_F_A=distance))
display(pd.DataFrame(pair_diagnostics).round(3))
pd.DataFrame(pair_diagnostics).to_csv(OUTPUT/'pair_frame_diagnostics.csv',index=False)
""")
md(r"""
## 8. Planned comparisons and unresolved choices

| Track | Fidelity / pair constraint | Objective | Current status |
|---|---|---|---|
| Independent GCS | Continuous $F(P,P')\leq\delta$ separately | Minimum vertices per fragment; endpoints fixed | Implemented, checked, and run on every prepared fragment |
| Joint discrete CPS-3F | Three discrete Fréchet bounds | Minimum $\max(|A'|,|B'|)$ | Exact small-window reference; 21 full-pair job specifications prepared; scalable Algorithm 2 pending |
| CPS-2F | To be specified with additional data | To be specified | Named interface deliberately raises `NotImplementedError` |

**Frame diagnostic:** the supplied 1O7J A and 1HFJ C coordinates have a discrete Fréchet distance of about **54.87 Å** in their deposited frames, whereas Table 1 uses a pair threshold of 1 Å. Their present placement therefore cannot be treated as the paper's registered experiment. The window reference can correctly return infeasible; that is not a reason to silently align or relax the bounds.

**Comparison plan.** Start with identical prepared fragment pairs and fixed, recorded registration transforms. Sweep the archived paper thresholds separately from the exploratory 2/4/8 Å sweep. For each method report retained counts for both chains, their maximum and sum, both continuous and discrete fidelity measures, simplified-pair discrete Fréchet distance, feasibility, time, and memory. Report the vertex-placement and endpoint rules. Add an independent **discrete-fidelity** baseline to isolate the benefit of joint optimization from the change in metric. Compare against paper counts only after historical selection and registration are resolved. Preserve explicit infeasibility rather than increasing tolerances silently.

**Unresolved choices, with current behavior stated:**

- **Historical data and typo confirmation:** the printed identifiers are preserved. Their interpretations are supported by PDB records and exact local-coordinate matches, but author confirmation and the historical subsets are unavailable. Current Table 3 C-alpha selections are used without trimming.
- **Registration and biological correspondence:** deposited coordinates are used; no rigid alignment is inferred. Choose a correspondence and fitting method, or explicitly retain these frames. Freeze the same transform before comparing methods; a visualization's centroid shift is never reused for pair distances.
- **Fragment pairing:** a broken chain is a collection, not one curve. Decide homologous residue windows and how to pair fragments, including unmatched termini and fragments. Full-pair jobs remain pending until that choice is made.
- **Endpoint fairness:** independent GCS fixes each fragment's endpoints. The paper's CPS does not. Report both policies for CPS or formulate an additional independent baseline with compatible endpoint rules.
- **Structural variants and sensitivity:** the selected model, raw alternate choices, and break thresholds are recorded. Decide whether to restore the omitted 4ECA residue 12 in an additional dataset, evaluate other models/conformers, and vary the conservative break rules. Do not change the primary raw dataset silently.
- **CPS-2F definition:** specify exactly which two Fréchet constraints, whether each is continuous or discrete, the coupling and vertex domains, objective, endpoints, and treatment of gaps. Its name alone does not determine an algorithm.
- **Scale and numerical limits:** full CPS needs a scalable implementation and a runtime/memory budget. The GCS solver uses floating point; threshold-adjacent results should be checked with tighter/robust predicates if they determine a scientific conclusion. Editing inline solver cells requires `RECOMPUTE_ALL=True`; cache signatures cover the shipped modules and input data, not arbitrary cell edits.

### Sources and reproducibility

- Supplied full GCS paper: `papers/global curve simplification.pdf`, Section 3.3, Algorithm 1, printed page 13. This notebook implements its global recurrence with a different state representation; it does not use the differently numbered conference Algorithm 1.
- Fan, Filtser, Katz, Wylie, and Zhu, [_On the Chain Pair Simplification Problem_, full version](https://arxiv.org/pdf/1409.2457), Sections 4–5 and 7, Tables 1–3. Cached as `papers/chain_pair_simplification.pdf`.
- Primary coordinate records: [RCSB PDB downloads](https://www.rcsb.org/downloads); individual entry links are in `output/identifier_audit.csv`, and source URLs, retrieval times, and SHA-256 hashes are in `data/pdb_cache/sources.json`. Author and label identifiers are read directly from the mmCIF records.

Prepared coordinates and residue maps are under `data/prepared/`; all audits, paper parameters, pair jobs, algorithm checks, selected residue mappings, and figures are under `output/`. Dependencies and execution instructions are in `README-analysis.md`.
""")
nb=nbf.v4.new_notebook(cells=cells,metadata=dict(kernelspec=dict(display_name='Python 3 (protein analysis)',language='python',name='python3'),language_info=dict(name='python',version='3.14')))
nbf.validate(nb)
nbf.write(nb,ROOT/'notebooks/legacy/0.analyze_independent.ipynb')
print(f'Wrote {len(cells)} cells.')
