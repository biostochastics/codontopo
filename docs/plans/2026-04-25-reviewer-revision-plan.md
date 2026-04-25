# 2026-04-25 — Reviewer-driven revision plan (3 reviews + methodological nuance)

This plan addresses three converging review documents (Review-1 "Overall verdict", Review-2 "RIGORA Deep", Review-3 "Diff v1→v2") plus the methodological nuance about per-table validation potentially testing standard-code proximity rather than independent variant optimality.

## Triage summary

| # | Issue | Source | Action class | Touches code? |
|---|-------|--------|--------------|---------------|
| 1 | K_4^3 / K_3^4 → H(3,4) = K_4^□3 | R1.1 | Manuscript-only notation sweep | No |
| 2 | Candidate-space arithmetic 64×20=1280 with explicit y∈A_20∪{Stop} | R1.2 | Manuscript Methods rewrite (+ 1-line code docstring fix) | docstring only |
| 3 | Topology-breaking definition harmonization | R1.3 | Manuscript text — code already implements "creates new disconnection in previously-connected family" matching JSON 6/931 | docstring/method label |
| 4 | n=27→n=28 in Methods §2.3.4 | R1.4, R3.1 | Manuscript text fix | No |
| 5 | ρ-sweep n=1,000→n=10,000 with exact p-values | R1.5 | Manuscript text fix; values already in JSON | No |
| 6 | Bonferroni-across-38-tests false claim | R1.6 | Manuscript Limitations rewrite (family-wise statement) | No |
| 7 | Abstract conditional-logit wording (tRNA proxy not independent) | R1.7, R3.4 | Manuscript abstract rewrite + dynamic ΔAICc binding | No |
| 8 | tRNA enrichment Option A vs B reframing | R1.8 | Manuscript section reframing, optional code subset analysis | optional |
| 9 | Add table 22 to Leu disconnection catalogue | R1.9 | catalogue.py + claim_hierarchy.py + manuscript prose | YES |
| 10 | Dual-function NCBI codons handling sentence | R1.10 | Manuscript Methods sentence | No |
| A | Topology-definition sensitivity table (Δβ₀>0 vs new-disconnection) | R1.A | New code: topology_avoidance.py — add Δβ₀>0 variant | YES |
| B | All-24-encoding sensitivity for Q_6 topology avoidance | R1.B | New code: encoding-sweep over topology | YES |
| C | Leave-one-clade conditional-logit refits | R1.C, R2.M1, R3.M2 | New code: evolutionary_simulation.py — clade-exclusion regime | YES |
| D | Denominator sensitivity table | R1.D | New code: report 4 candidate universes | YES |
| E | Direct tRNA subset analysis for topology-breaking cases | R1.E | New code: trna_evidence.py — topology-restricted subset | YES |
| F | Posterior-predictive validation surfacing | R1.F | Manuscript text only — already in JSON | No |
| 2A | Syn57 Δlocal Mann-Whitney degenerate test | R2.M2 | Manuscript Section 3.9.1.2 reframing | No |
| 2B | IIA assumption surfacing in main text | R2.M3 | Manuscript Section 4.2 sentence | No |
| 2C | Tables 1/11 sense-codon-identity disclosure | R2.5 | Manuscript per-table-results paragraph | No |
| 2D | Per-table standard-code-proximity caveat | R2.6, "methodological nuance" | Manuscript paragraph + optional encoding-shuffled control | YES (control) |
| 2E | Syn61 7-corrections characterization | R2.7 | Manuscript paragraph | No |
| 2F | Ding 2024 framing in Table 8 | R2.8, R3.* | Manuscript Table 8 caption + cross-kingdom-scope paragraph | No |
| 2G | catalogue.json WS2-DIR mismatch | R2.9 | catalogue.py — demote WS2-DIR to exploratory + p=0.075 | YES |
| 2H | ΔtRNA proxy-vs-effect interpretation | R2.10 | Manuscript Section 3.5 sentence | No |
| 2I | Mechanistic discriminant Miyata-vs-Grantham | R2.11 | Manuscript exploratory section | No |
| 2J | Yeast-mito triple-convergence note | R2.12 | Manuscript Section 4.2 paragraph | No |
| 3A | Stale ΔAICc=100/67 in Discussion 4.2 | R3.2 | Already-dynamic — verify | partial |
| 3B | Figure 4 caption ΔAICc≥67 | R3.3 | Already-dynamic — verify | partial |
| 3C | Section 3.5 inline figure values | R3.5 | Re-render Fig 4 from current JSON | YES (R re-render) |
| 3D | catalogue.json WS2-DIR | R3 (overlap with 2G) | (see 2G) | YES |
| Misc | Reference formatting (doubled periods, "et al.") | R1 secondary | references.bib cleanup | bib only |
| Misc | codon-topo version 0.3.1 in supplement S8 | R1 secondary | supplement.typ | No |
| Misc | Syn57 60,005 vs 60,240 reconciliation | R1 secondary | Manuscript footnote | No |
| Misc | "evolution favors" → softer | R1 secondary | Manuscript prose | No |
| Misc | Table 7 organism display | user request | Verify R figure / supplement table | YES (verify) |

## Phases

**Phase 0 — Pre-flight** (informational)
- Verify exact current state of `topology_avoidance.json` Q6 numerator/denominator semantics → confirmed: it implements "creates new disconnection from a previously connected family" (6/931, 6/846). The Methods text says Δβ₀>0. **Decision:** keep code semantics (matches existing numbers + manuscript Table 4); update Methods text to match the implemented semantics OR add Δβ₀>0 as a sensitivity column.
  
  Per reviewer R1's preference for Δβ₀>0 (matches conditional-logit Δ_topo): **add Δβ₀>0 as a parallel column in topology_avoidance.json + Table 4** so both definitions are present. The "creates new disconnection" definition stays primary; the Δβ₀>0 column is the sensitivity.

**Phase 1 — Code (analysis modules)**
1. Catalogue/claim-hierarchy reconciliation:
   - `src/codon_topo/reports/catalogue.py`: WS2-DIR → status="exploratory", evidence_strength="weak", p_value=0.075 (codon-permutation null gives p=1.0; bare positional p=0.006 was misleading)
   - `src/codon_topo/reports/claim_hierarchy.py`: confirm Leu disconnection catalogue entry references both table 16 AND table 22 (algal mito clade)
   - Add catalogue entry note for table 22 if missing
2. Topology-avoidance Δβ₀>0 sensitivity:
   - `src/codon_topo/analysis/synbio_feasibility.py` (or wherever topology landscape lives) — add `_count_disconnection_increases` returning (possible_breaks_β0, observed_breaks_β0) under both Q_6 and K_4^3
   - Output to `topology_avoidance.json` as `Q6_beta0` and `K43_beta0` keys
3. Q_6 topology avoidance encoding-sweep:
   - For each of 24 base→bit bijections, recompute (rate_observed, rate_possible, depletion, hypergeom_p)
   - Output to `topology_avoidance.json` under `Q6_encoding_sweep` key
4. Leave-one-clade conditional-logit refits:
   - `src/codon_topo/analysis/evolutionary_simulation.py` — add `fit_clade_excluded_models()` that loops the 7 clade-exclusion regimes from `phylogenetic_sensitivity.json` and refits M1-M4 on each
   - Output to new `condlogit_clade_sensitivity.json`
5. tRNA topology-breaking subset:
   - `src/codon_topo/analysis/trna_evidence.py` — add `topology_breaking_subset()` that filters to events that disconnect a previously-connected family under K_4^3 (the 6 events: Thr/T3, Leu/T16, Leu/T22, Ala/T26, Ser/T12, plus one TBD)
   - Output `trna_topology_subset` block in `trna_evidence.json`
6. Denominator-sensitivity table: in `topology_avoidance.json`, add `denominator_sensitivity` block reporting (1280 21-label-no-id, 1219 AA-only-no-id, 1280 AA-with-noop, 1344 stop-incl-with-noop) counts and p-values
7. Standard-code-proximity control (R2.6 / methodological nuance):
   - `src/codon_topo/analysis/coloring_optimality.py` — add `per_table_proximity_null()` that fits the same per-table block-preserving null but **excluding** standard-code-identical permutations (or reporting Hamming distance from standard for null sample)
   - This addresses the "near-standard-code" concern — if the variant code's quantile remains low after conditioning out proximity-to-standard-code, the per-table optimality is independent.
   - Output `per_table_proximity_audit` block in `manuscript_stats.json`
8. Re-render `manuscript_stats.json` by re-running `codon-topo all`
9. Re-run R figure scripts (Figures 1-5, plus any sensitivity panels)

**Phase 2 — Manuscript (output/manuscript.typ)**
Mechanical/notation:
1. K_4^3 / K_3^4 → H(3,4) = K_4^□3 (sweep all instances; define once at first use; keep `K_4^3` as a parenthetical alias the first time only)
2. Methods §2.3.4: n=27→n=28; K=931 → K=931 (Q_6) and K=846 (K_4^3); explicit candidate-space `M(C) = {(x,y): x∈codons, y∈A_20∪{Stop}, y≠C(x)}, |M(C)|=1280`
3. Methods §2.3.3: ρ-sweep n=1,000 → n=10,000; cite exact p-values
4. Limitations §4.6: replace "Bonferroni across 38 tests = 1.3e-3" paragraph with family-wise statement (R1.6 wording)
5. Abstract: ΔAICc≥89 phrasing (already done) + clarify tRNA proxy as "does not improve" not "adds independent power"; verify
6. Section 3.4 / Table 4: add Δβ₀>0 sensitivity column (or note as supplementary)
7. Section 3.4: add denominator-sensitivity reference + clade-exclusion summary
8. Section 3.5: surface IIA caveat (R2.M3); surface posterior-predictive p=0.60 (R1.F); add condlogit clade-exclusion sensitivity reference
9. Section 3.9.1.2 Syn57 Δlocal: reframe as descriptive — "by construction the alanine arm has Δlocal=+19.0 ± 0.0 (zero variance), so the Mann-Whitney p-value is a property of the design rather than of the underlying signal; the inferential weight comes from the Napolitano arginine analysis (Spearman ρ=−0.33 vs RBS deviation, n=12,888)"
10. Section 4.2: add IIA caveat sentence + yeast-mito triple-convergence paragraph
11. Section 3.x (Results — per-table): add Tables 1/11 and 27/28 sense-codon-identity disclosure + standard-code-proximity caveat with reference to proximity-audit supplement
12. Disconnection catalogue prose: "Leu in chlorophycean mitochondrial code (translation tables 16 and 22)" rather than "16" alone; explicitly state the lineage-collapsed event count
13. Methods: add dual-function NCBI codons handling sentence (R1.10)
14. Section 3.10 (KRAS or wherever exploratory): one-sentence Miyata-vs-Grantham asymmetry interpretation
15. Section 3.9.1.1 Syn61: characterize the 7 corrections (or remove the "99.96% viable" rhetorical use)
16. Table 8: rewrite Ding row note ("included qualitatively for cross-kingdom scope")
17. Section 4.x: soften "evolution favors" → "observed reassignment histories are enriched for..."
18. Reference formatting: clean references.bib doubled-period-before-DOI bug; replace "others" → "et al."
19. Supplement S8: align codon-topo version to 0.3.1
20. Syn57 60,005 vs 60,240: add footnote — "remainder = stop-codon recodings (n=235)"
21. Table 7 (or wherever): verify all 18 organisms display correctly per user request

**Phase 3 — Verification**
1. `python3.11 -m pytest` — all tests pass
2. `codon-topo all --output-dir=./output` — full pipeline regen
3. Re-render figures (R scripts)
4. `typst compile output/manuscript.typ output/manuscript.pdf` — clean compile
5. PyMuPDF text extraction — grep rendered PDF for stale numbers
6. `typst compile output/supplement.typ output/supplement.pdf` — clean compile

**Phase 4 — xen consensus** (user-requested)
Consult mimo, glm-5.1, kimi-k2.5, minimax-m2.7, gpt-5.2-pro on:
- Methodological nuance: is the per-table block-preserving null a defensible test of variant-code independent optimality, or does standard-code proximity dominate? What is the minimal valid framing?
- Overall verdict (Review-1): are the proposed mechanical fixes sufficient, or are there deeper structural issues we have missed?

## Open questions
- Confirmed: manuscript already uses dynamic ΔAICc bindings (line 114 abstract, line 504 Fig 4 caption, line 720 Discussion). Stale-prose risk is low if JSON is re-rendered.
- Confirmed: depth_calibration.py already includes table 22 (Scenedesmus obliquus mito Leu). Catalogue text (manuscript line 576) needs updating.
- Confirmed: phylogenetic_sensitivity.json already has yeast-mito clade-exclusion drop (21.4% → 8.3%). Just needs surfacing.

## Status tracker

See associated TaskCreate entries.
