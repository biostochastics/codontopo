# Clayworth Codon-Topo Explainer — Master Build Prompt v2

**Purpose**: Drop this entire block into a fresh Claude Code session (or v0.dev, or Cursor) to build an interactive web explainer that walks a scientific reader through the REFINED findings of the Clayworth codon-topo analysis.

**Stack**: TanStack Start v1.x • Mantine UI v8 • Recharts v3 • react-three-fiber • d3-force • Zustand v5 • Motion • TypeScript strict • KaTeX

---

## MASTER PROMPT (copy-paste into new session)

```
I want to build a polished, interactive explainer that walks a scientific
audience through the findings of the Clayworth Codon-Topo Framework analysis.
This framework proposed that the genetic code has algebraic structure
visible in GF(2)^6 encoding. Multi-model adversarial review (10 LLMs, 12
passes) identified which claims are real, which are tautological, which are
falsified by counterexample, and which were mistested. Several refinements
were then implemented that changed the story significantly.

The reader — likely a mathematical biologist or computational biology
grad student — should come away understanding:

  1. What claim was originally made
  2. What specifically broke
  3. What was fixed
  4. What genuinely survives as a publishable finding

This is an honest scientific post-mortem with live, interactive math.
Tone: rigorous, even-keeled, distill.pub style.

## STACK (strict)

- TanStack Start v1.x (file-based routing, SSR, Vercel deploy target)
- Mantine UI v8 (all components, typography, theme, dark mode default)
- Recharts v3 for 2D plots (histograms, bar, line)
- react-three-fiber + drei for 3D Q_6 hypercube visualization
- d3-force for graph layout computation (render with React)
- Zustand v5 for client state (encoding, filters)
- Motion (Framer Motion successor) for transitions
- KaTeX via better-react-mathjax for formulas
- TypeScript strict mode
- Monochrome + single accent (Mantine teal) palette, dark mode default

READ DOCS BEFORE CODING — don't rely on memorized APIs:
  https://tanstack.com/start/latest/docs/framework/react/overview
  https://mantine.dev/guides/typescript/
  https://docs.pmnd.rs/react-three-fiber/getting-started/introduction
  https://recharts.org/en-US/guide

## STRUCTURE

Single-page, scrollable, with sticky left navigation. AppShell from Mantine,
7 main claim sections plus intro + survivors. Each section has a status
badge (RED=broken, ORANGE=false-in-generality, YELLOW=mistest, GREEN=survives,
BLUE=deadly negative).

### Navigation structure

  0. Preamble — The encoding (interactive picker)
  1. M1: Two-fold bit-5 filtration                [RED — tautology]
  2. M2: Four-fold prefix filtration              [RED — trivial]
  3. M3: Serine min-distance invariant             [ORANGE — counterexample]
  4. M4: Fano lines / KRAS clinical                [BLUE — falsified]
  5. M5: Holomorphic embedding                    [RED — not a character]
  6. M6: Null Model A labeling                    [YELLOW — restate]
  7. M7: PSL(2,7) symmetry                        [RED — pre-rejected]
  -------------------------------------------------
  8. Bit-position bias (refined)                  [GREEN — p<0.05 after fix]
  9. tRNA duplication pattern                     [GREEN — trend 4/4]
  10. Hypercube coloring optimality (central)     [GREEN — p=0.006 proper null]
  11. What actually survives + paper outline

## THE KEY VISUALIZATIONS

### Vis 0 — Encoding Picker (preamble, persistent state)

Mantine SegmentedControl picks one of 24 permutations of {C,U,A,G} →
{(0,0),(0,1),(1,0),(1,1)}. Default: C=(0,0), U=(0,1), A=(1,0), G=(1,1).

Shows:
- 4×2 grid: base → 2-bit pair
- Parsed semantics:
  - Bit 1: "pyrimidine/purine" or "?" depending on encoding
  - Bit 2: "amino/keto" (for default) or "weak/strong" or "?"
  - H-bonding: "Bit 1 XOR Bit 2" parity (when default)

Zustand persists the selected encoding for all subsequent visualizations.

Add two preset buttons:
- [Default] — C=00, U=01, A=10, G=11
- [Kimi counterexample] — U=00, C=11, A=01, G=10

### Vis 1 — Q_6 Hypercube (used in M1, M2, M3)

3D rendering via react-three-fiber:
- 64 vertices = 64 codons, labeled as RNA triplets
- 192 edges = Hamming-1 adjacencies
- Color by amino acid (colorblind-safe: Tol bright palette)
- Hover reveals codon + AA
- Buttons: "Highlight only 2-fold AAs", "Highlight only 4-fold AAs",
  "Highlight only Serine (show disconnection)"
- Epsilon slider (1→6): edges appear/disappear as threshold changes
- Connected-component count readout updates live

Layout: use PCA of adjacency matrix OR a hand-tuned 3D projection
preserving the 6-dim hypercube structure (4×4×4 grid with connections).

### Vis 2 — Encoding vs Invariance (M1, M3)

Split panel:
- LEFT: Default encoding, synonymous pairs labeled "differ at bit X"
- RIGHT: User-selected (or kimi counterexample) encoding, same display

Below: a table "Under how many of the 24 encodings does 'bit 5 only' hold?"
showing 16/24. And "In how many encodings is Serine min distance = 4?"
showing 8/24 (the remaining 16 give min distance = 2).

Highlight the counterexample preset: clicking Kimi preset makes the Serine
panel flash red, showing min distance = 2 with AGU-UCU highlighted.

### Vis 3 — Linear Dependency Count (M4)

Show that XOR-zero triples in GF(2)^6 are NOT Fano lines — they are
2-dimensional subspaces of PG(5,2).

Display a scrollable list of 651 triples. Highlight GGU-GUU-CAC as "the
framework's example." A shuffle button picks a random triple and verifies
(arithmetic panel) that it XORs to zero. Nothing special about any one.

Math panel with KaTeX:
  Number of 2-dim subspaces in GF(2)^6 =
  $\binom{6}{2}_2 = \frac{(2^6-1)(2^6-2)}{(2^2-1)(2^2-2)} = \frac{63·62}{6} = 651$

### Vis 4 — KRAS Clinical Failure (M4)

Embed the cBioPortal MSK-IMPACT 2017 test result. Bar chart of 6 Fisher's
exact test p-values (one per G12 variant). All bars at p=1.0. Big red
"PREDICTION FAILED" label.

Below: "No biological mechanism propagates GF(2)^6 XOR through DNA
polymerase." Quote the Antoneli-Forger paper on why this was predictable.

### Vis 5 — Character Failure (M5)

Interactive diagram:
- (GF(2)², +) on the left (4 elements, order ≤ 2)
- ℂ* on the right with {1, i, -1, -i} highlighted
- Draw an arrow from (0,1) to i
- Compute live:
  χ((0,1) + (0,1)) = χ(0,0) = 1
  χ((0,1))² = i² = -1
  "Homomorphism condition fails: 1 ≠ -1"

Color the mismatch in red.

Subpanel: "The 4 actual characters of (GF(2)², +) take values in {+1, -1}."
Show all 4.

### Vis 6 — Null Model A Labeling (M6)

Code side-by-side:

```python
count_serine_unique = 0    # <-- variable name
if score["exactly_one_disconnected"]:    # <-- actual condition
    count_serine_unique += 1
```

Narration: "Under null permutations, Serine's label is randomized away.
This counts 'exactly one AA disconnected' — NOT 'Serine specifically.'"

### Vis 7 — Bit-Position Bias (refined, p<0.05 survives)

Two histograms stacked:
- TOP: Observed bit-flip distribution [5, 4, 0, 8, 13, 5]
- BOTTOM: Expected under various nulls:
  - Uniform null: [5.83, 5.83, 5.83, 5.83, 5.83, 5.83] → χ²=16.26, p=0.006
  - Nuclear Ts/Tv (2.5, 2.8, 4.5): weighted → χ²=11.13, p=0.049
  - Mitochondrial Ts/Tv (8, 6, 15): weighted → χ²=12.60, p=0.027

Annotation: "Signal weakens from p=0.006 to p=0.049/0.027, but STAYS
significant at p<0.05 under realistic nulls."

Segmented control for null-type, p-value updates live.

### Vis 8 — tRNA Duplication Trend (4/4)

Horizontal bar chart, 4 rows:
  S. cerevisiae (mito)   | Thr | 2 vs 1 (Yarrowia mito)
  S. obliquus (mito)     | Leu | 2 vs 1 (C. reinhardtii mito)
  Pachysolen             | Ala | 14 vs 11 (Lachancea)
  Candida albicans       | Ser | 16 vs 13 (Lachancea)

All bars show elevation in the disconnection organism. Binomial
p = 0.0625 (one-sided). "4/4 cases — trend level — one confirmed
mechanistically (Su et al. 2011 PMC3113583)."

Callout: "Yarrowia lipolytica replaced Chlamydomonas as the yeast mito
control after reviewer feedback — same phylogenetic subphylum."

### Vis 9 — Hypercube Coloring Optimality (CENTRAL)

Full-width chart. Show the Monte Carlo null distribution (histogram of
F values from 10,000 random codes under Freeland-Hurst block-preserving
null). Red vertical line = standard code F = 13,477. Label the quantile
(0.60%) and p-value (0.006).

Side panel: toggle null type
- Freeland-Hurst (block-preserving, stops fixed): p=0.006
- Class-size-preserving (weaker): p<0.001 (bounded by n=10000)

Second side panel: toggle include_stops
- With stops: p=0.006, F=13477
- Without stops: p=0.006, F=10467
- Stop contribution is CONSTANT (confirming fixed-stop null is correct)

Third side panel: Grantham matrix heatmap (20×20, amino acid physicochemical
distances). Hover shows the value. Leu-Ile=5 (closest), Trp-Cys=215 (farthest).

Text below: "Under the proper Freeland-Hurst null that preserves synonymous
block structure, the standard genetic code sits in the top 0.6% of random
colorings for edge-mismatch minimization (p=0.006). This is the proper
refinement of Freeland & Hurst (1998)'s 'one in a million' finding in
GF(2)^6 coordinates — stricter null, therefore smaller apparent p."

### Vis 10 — Cross-Table Optimality

Scatter plot: each of 25 NCBI tables as a point. X = table ID, Y = F score
(or F/F_standard). Standard is at 100%. Most variant codes are 97-108%.

Annotation: "Variant codes with reassignments don't dramatically improve
or degrade optimality. The structure is conserved."

### Vis 11 — The Survivors (summary)

Paper outline card with three columns:
  SUPPORTED            | SUGGESTIVE          | REJECTED
  ===================  | ===================  | =========
  Hypercube coloring   | tRNA duplication     | Serine ε=4 invariant
  optimality (p=0.006) | pattern (4/4, p=0.063)| Fano clinical predictions
                        | Bit-position bias    | Holomorphic embedding
                        | (p<0.05 weighted)    | PSL(2,7) symmetry
                                                | Depth calibration (n=6)

Journal target: Journal of Theoretical Biology (primary), BMC Bioinformatics
(backup). Not PLoS Comp Bio in current form (per gpt-5.4-pro review).

## DESIGN LANGUAGE

- Inter or system UI for body; SF Mono / JetBrains Mono for codons
- Dark mode by default, Mantine colorScheme switcher available
- Codons inline in mono: `UCU`, `AGY`, `GGU`
- KaTeX for all math expressions
- Status badges with Mantine Badge, colored by category:
  RED (broken-tautology), ORANGE (false-in-generality), BLUE (deadly-negative),
  YELLOW (mistest-but-fixable), GREEN (survives)
- Mantine Paper with radius="md", shadow="sm" for each section
- At the end of each section: "Model Consensus" strip showing the 10 LLM
  evaluations that examined this claim (pill-style badges with model name
  and verdict emoji)

## STATE SHAPE (Zustand)

```ts
interface StoreState {
  encoding: { C: [number, number]; U: [number, number];
              A: [number, number]; G: [number, number] };
  setEncoding: (e: Encoding) => void;
  activeClaimId: string;
  epsilon: number;
  setEpsilon: (e: number) => void;
  nullType: "freeland_hurst" | "class_size";
  includeStops: boolean;
}
```

Persist to localStorage.

## SERVER ROUTES (TanStack Server Functions)

- POST /api/compute/hamming — takes encoding + codon set, returns distances
- POST /api/compute/monte-carlo — takes null_type, n, returns histogram
- GET /api/data/grantham — returns the 20×20 matrix
- GET /api/data/ncbi25 — returns all 25 NCBI translation tables
- GET /api/data/trna — returns curated tRNA repertoires

## DATA FILES (static JSON in /app/data)

- ncbi_25_tables.json (all 25 NCBI translation tables)
- grantham.json (verified against Grantham 1974 Table 1)
- reassignment_events.json (61 events with bit-flip data)
- trna_repertoires.json (6 curated repertoires: 4 disconnection + 2 controls)
- kras_msk_impact.json (cBioPortal pre-extracted KRAS data)
- precomputed_null_distributions.json (for fast UI replay without
  client-side Monte Carlo)

## MILESTONES (iterative)

1. Bootstrap TanStack Start + Mantine + routing + dark mode theme (no
   business logic yet)
2. Vis 0 encoding picker + Zustand store
3. Vis 9 hypercube coloring (centerpiece — prove the central result)
4. Vis 1 Q_6 hypercube 3D
5. Vis 2 encoding vs invariance (shows the counterexample)
6. Vis 7-8 bit-bias + tRNA
7. Remaining claims (M1-M7)
8. Polish: transitions, keyboard nav, a11y fallback for 3D
9. Deploy to Vercel with TanStack Start preset

## WHAT I WILL PROVIDE YOU

- The 4 refined Python modules (coloring_optimality.py, trna_evidence.py,
  reassignment_db.py, null_models.py)
- Precomputed null distributions as JSON
- The 25 NCBI tables as JSON
- KRAS data subset from cBioPortal
- All 12 model-consensus findings for the status strips
- Freeland-Hurst 1998 paper and citations

## WHAT NOT TO DO

- Don't use D3 selection API; use d3-force just for layout math, render with React
- Don't use MUI, Radix, or Tailwind — Mantine v8 covers everything
- Don't use Redux; Zustand only
- Don't use SWR; TanStack Query ships with TanStack Start
- Don't style inline — use Mantine style props + CSS modules
- Don't over-emoji the UI — this is a serious scientific post-mortem
- Don't memorize APIs — read docs for each library first
- Don't misreport p-values as exactly 0; always show conservative (k+1)/(n+1)
- Don't claim "one in a million" for Freeland-Hurst unless running 1M+ samples

## CRITICAL COPY CORRECTIONS TO AVOID

Avoid language like:
  - "algebraic law" or "discovery of algebraic structure"
  - "universal topological invariant" (unless qualified to ε=1 only)
  - "Fano lines" (use "2-dim subspaces of PG(5,2)" or "3-term linear
    relations")
  - "holomorphic embedding" (use "coordinate-wise root-of-unity map")
  - "PSL(2,7) symmetry" (drop — pre-rejected)
  - "one in a million" (Freeland-Hurst magnitude requires n≥10^6)

Prefer language like:
  - "The code is significantly error-minimizing under a block-preserving null"
  - "Suggestive evidence for tRNA duplication accompanying reassignment"
  - "Minimum distance = 4 holds in 8/24 encodings, = 2 in 16/24"
  - "The ε=1 disconnection is the only invariant surviving encoding
    permutation"

Start by building Milestone 1 and show me the result. We'll iterate.
```

---

## Notes for when you run this

### On running the prompt

1. **Seed data first**: After the prompt, drop in the four JSON files:
   `ncbi_25_tables.json`, `grantham.json`, `reassignment_events.json`,
   `trna_repertoires.json`. Let the new Claude discover these in context.

2. **Vis 9 is the centerpiece**: Prioritize the hypercube coloring
   optimality visualization. If nothing else lands, that one tells the
   whole "corrected" story. Get it working with live Monte Carlo replay
   from precomputed histogram.

3. **Vis 2 tells the counterexample story**: Clicking the Kimi preset
   button should make Serine's panel flash red. This is the
   single-clearest demonstration of why the "invariance" claim broke.

4. **The honesty framing is the sell**: This isn't a paper that
   discovered something; it's a paper that shows careful refinement under
   adversarial review yielded two publishable findings (hypercube
   optimality + tRNA pattern) while falsifying several others. That
   honesty framing is what makes the post-mortem interesting to read.

5. **Performance**: Hypercube 3D with 64 vertices / 192 edges is trivial
   for three-fiber. Animations on epsilon slider should be smooth.
   Precomputed null distributions keep the Monte Carlo viz instant.

6. **Accessibility**: 3D Q_6 viz needs a "show me the table" toggle for
   screen readers (renders the same info as a Mantine Table).

7. **Deployment**: TanStack Start + Vercel via `@tanstack/react-start`
   framework preset. Should work out of the box with zero config.

### On the "proper" replication framing

gpt-5.4-pro's key insight: our p=0.006 under block-preserving null is a
stricter test than Freeland-Hurst 1998's looser null that gave p≈10^-6.
The paper should NOT claim "we reproduced one-in-a-million" — it should
claim "we reproduced the phenomenon under a more conservative null."

That's the right framing for the explainer too. The narrative is:
"Stricter null, still strongly optimal. The code is error-minimizing
under *realistic* constraints, not just toy randomizations."

### On target journal

Per gpt-5.4-pro: **Journal of Theoretical Biology** is the right fit.
The explainer's tone should be consistent with JTB audience — math-tolerant,
computational biology background, appreciates rigor.
