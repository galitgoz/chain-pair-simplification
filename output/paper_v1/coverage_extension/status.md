# Coverage extension: saved results

## Old timing study
120/120 scheduled repetitions recorded: 117 validated, 3 transition-budget exits. The four previously pending repetitions completed successfully. The three CPS-2F AL072006__AL011900 exits are preserved; no extra timing rounds or replacement attempts were launched.

## Protein-paper coverage
Tables 1–3 of the local paper were extracted and checked (see provenance/paper_tables.json for PDF pages and hash). All 14 requested B identities are inventoried; 12 pairs have eligible inputs. Two pairs, 3ntx.a and 3jq0.a, are preparation-blocked because the existing policy splits their internally broken chains. No cropping, invented residue matching, or new gap bridging was introduced.

The original five protein inputs remain unchanged. 4eca.b/c reuse supplied aligned w=16 arrays with documented residue omissions and 0.001-angstrom coordinate storage. The five eligible varying-length chains reuse existing US-align structural_fit on all selected observed residues before independent endpoint-preserving w=16 sampling. Their B input lengths vary from 15 to 29; A retains 22 vertices. This preparation differs from historical index-Kabsch inputs, so cross-pair changes cannot be attributed to size alone.

Identifier corrections: printed 107j.a maps to verified 1o7j.a; 4cea.b/d to 4eca.b/d; 1toh to author chain A. Paper 1d9q.d length 297 differs from local selected/deposited C-alpha length 325. Other length discrepancies, author versus label chains, alternate conformers, model selection and gaps are recorded in protein_paper_coverage.csv and protein_audit/. Identity coverage is not exact reproduction of the paper's coordinates, thresholds, or results.

## New hurricane selection
Earlier additions used size-quantile bands within the prior maximum-track-size range, unused-storm preference, and geometry diversity. The final20 cohort instead uses seeded, preference-constrained stratified random sampling from the entire eligible supplied-pair pool, without a size cap; see its separate report and selection_pool.json. No solver outcomes entered either selection procedure. All supplied vertices and the shared basin frame are retained; temporal synchronization is not inferred. Cohort identity is retained in the combined input and result tables. Unweighted final20 percentages describe that stratified cohort, not the archive.

- AL182011__AL061926: 31/27 vertices.
- EP141997__EP092000: 34/25 vertices.
- EP211990__EP152008: 26/38 vertices.
- EP041981__EP211992: 25/40 vertices.
- AL181887__AL091888: 26/64 vertices.
- EP172000__EP061970: 32/29 vertices.
- EP111995__EP032014: 41/41 vertices.
- EP061986__EP131991: 38/70 vertices.
- EP092002__EP171994: 26/31 vertices.
- AL091944__AL031917: 29/35 vertices.
- EP101985__EP151985: 37/38 vertices.
- EP031965__EP141991: 25/42 vertices.
- AL061921__AL101926: 43/57 vertices.
- EP041976__EP041984: 30/48 vertices.
- EP081967__EP101980: 25/26 vertices.
- EP031994__EP031978: 29/35 vertices.
- EP021985__EP181984: 41/36 vertices.
- AL101886__AL051950: 25/38 vertices.
- AL091947__AL011953: 36/53 vertices.
- EP072011__EP092011: 26/50 vertices.
- AL171980__AL112009: 25/31 vertices.
- AL102008__AL061989: 33/36 vertices.
- EP102002__EP042014: 31/42 vertices.
- EP032000__EP041975: 27/41 vertices.
- EP081993__EP121989: 56/45 vertices.
- EP032006__EP201992: 27/69 vertices.
- AL022011__AL071977: 29/27 vertices.
- EP101972__EP141967: 29/37 vertices.

## Completion
- New observations: 140/140 executable slots; 139 validated; 1 resource exits; 8 preparation-blocked slots; 0 pending executable slots.
- Combined: 180 observed of 188 intended slots (180 eligible); 178 validated; 2 resource exits. Blocked/unexecuted slots are not observations.
- Original 40 single observations, new extension observations and 120 timing repetitions retain separate provenance. No selection of favorable attempts.

