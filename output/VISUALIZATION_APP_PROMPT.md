# Prompt for Building the Clayworth Claims Explainer App

**Intended target**: A fresh session of Claude Code or v0.dev where you want to build a polished, interactive web explainer that walks a reader through each claim, showing *exactly* what's broken and why — with live visualizations.

**Stack**: TanStack Start (React meta-framework with SSR) + Mantine UI v8 + Recharts (or D3 for the hypercube) + Zustand (for state) + Motion (Framer Motion successor) for animations.

Copy the block below as the opening prompt to your new session.

---

## MASTER PROMPT (copy-paste into a fresh session)

```
I want to build a polished, interactive web explainer that walks a scientific
audience through the mathematical and biological claims of the "Clayworth
Codon-Topo Framework" — a framework that claims to find algebraic structure in
the genetic code using GF(2)^6 encoding. Multi-model adversarial review has
shown that most of these claims are tautological, false-in-generality, mistested,
or redundant. This app presents each claim, the critique, and live interactive
visualizations that let the reader verify the math themselves.

The tone is rigorous but friendly — like a really good Distill.pub post.

## STACK

- TanStack Start v1.x (React with file-based routing, SSR)
- Mantine UI v8 (all layout, forms, typography, theme)
- Recharts v3 for 2D plots
- react-three-fiber + drei for 3D hypercube visualization (Q_6 projected)
- d3-force for network/graph layouts
- Zustand v5 for client state (encoding selection, filters)
- Motion (Framer Motion successor) for page transitions and claim reveals
- TypeScript strict mode
- shadcn monochrome + one accent color (mantine purple or teal)

IMPORTANT: I'll provide package versions and API docs as we go. Do NOT rely on
memorized API shapes. Read docs before writing code for each library.

## STRUCTURE

Single-page app with sections that scroll, using Mantine's AppShell with a
persistent left navigation listing each claim. As the user scrolls, the nav
highlights the current claim.

### Top hero
- Title: "What Actually Holds Up: The Clayworth Codon-Topo Framework Under 10-Model Adversarial Review"
- Subtitle: "A claim-by-claim verification using interactive mathematics"
- Badge row: "12 model evaluations", "280 tests pass", "1 counterexample verified"

### Navigation (sticky left)
Sections:
  1. Preamble: The Encoding
  2. Claim M1: Two-fold bit-5 filtration (TAUTOLOGICAL)
  3. Claim M2: Four-fold prefix filtration (TRIVIAL)
  4. Claim M3: Serine disconnection min-distance-4 (FALSE IN GENERALITY)
  5. Claim M4: Fano lines (TERMINOLOGY WRONG + BIOLOGY DEAD)
  6. Claim M5: Holomorphic embedding (NOT A CHARACTER)
  7. Claim M6: Null Model A labeling (STATISTIC MISMATCH)
  8. Claim E1: Bit-2 biochemistry (AMINO/KETO, NOT WEAK/STRONG)
  9. Biology claims B1-B7 (summarized)
  10. What survives: the hypercube coloring theorem
  11. Open hypotheses to test

Each numbered section is a `Paper` component with:
  - Section label (e.g., "CLAIM M3")
  - Status badge (colored: red = broken, yellow = salvageable, green = survives)
  - Original claim (blockquote)
  - "What's broken" header
  - Inline reasoning
  - THE INTERACTIVE VISUALIZATION (the centerpiece)
  - Model consensus block (small pill-row showing which models agree)
  - "Verify it yourself" expandable code snippet

## THE KEY VISUALIZATIONS

### Vis 1 — The Encoding Picker (preamble)
Mantine SegmentedControl lets the user select one of 24 base->bit-pair encodings.
Shows a 4x2 grid: base (C, U, A, G) -> 2-bit pair. Shows parsed meaning:
  - "Bit 1 groups: [pyrimidine/purine]"
  - "Bit 2 groups: [amino/keto | weak/strong | something else]"
  - "H-bonding as parity: [yes/no]"
This encoding choice persists across subsequent visualizations via Zustand store.

### Vis 2 — Q_6 Hypercube (centerpiece of M1, M2, M3)
Interactive 3D rendering using react-three-fiber. Project Q_6 using a
principal-component layout where the first three singular vectors of an
adjacency matrix give a 3D embedding.
  - 64 vertices (codons), labeled as RNA triplets
  - 192 edges (Hamming-1 connections)
  - Colored by amino acid (use color-blind-safe palette: Tol bright)
  - Hover reveals codon + amino acid
  - Button row to toggle: "Show only 2-fold AAs" / "Show only 4-fold AAs" /
    "Show only Serine (highlight disconnection)"
  - Slider for epsilon (Hamming threshold) — edges show/hide as you drag
  - Readout: "Number of connected components: X" updates live
  - Encoding selector above affects positions/colors

Critical animation: when user drags epsilon from 1 → 4 with Serine highlighted,
they see the two blocks visibly join via a bridge edge.

When user switches to kimi's counterexample encoding (preset button), Serine's
bridge forms at epsilon=2, demonstrating the counterexample lives.

### Vis 3 — Bit-Position Sensitivity (M1)
Two parallel panels:
  LEFT: Default encoding, 2-fold families connected by an edge with label
    "differ at bit X". All should say "bit 5".
  RIGHT: User-selected encoding, same display, possibly different bit labels.
Below: "Under how many encodings does 'differ only at bit 5' hold?" X/24.

### Vis 4 — Count Visualization (M4)
Show the 2-dim subspaces of GF(2)^6 as a scrollable list of 651 triples.
Highlight GGU/GUU/CAC as "the framework's example". Show there's nothing
special about it — hit "Shuffle" to pick a random one, verify it also XORs
to zero.

Subplot: Gaussian binomial coefficient calculation laid out mathematically
with MathJax/KaTeX, showing 651 is exact.

### Vis 5 — Character Failure (M5)
Interactive diagram of (GF(2)^2, +) mapping to C*. User can pick one of
the framework's claimed bijections (e.g., (0,1) -> i). Show:
  "χ((0,1) + (0,1)) = χ(0,0) = 1"
  "χ((0,1))² = i² = -1"
  "1 ≠ -1 → NOT a homomorphism"

Color the clash in red.

Side panel: "If we instead require values in {+1, -1} (valid characters of
(GF(2)^2, +)), here are the 4 genuine characters of the additive group."

### Vis 6 — Serine Counterexample (M3)
This is THE most important viz. Side-by-side:
  LEFT: Default encoding (C=00,U=01,A=10,G=11)
    - Show UCN block vertices in blue
    - Show AGY block vertices in green
    - Min distance between blocks: 4
    - Epsilon slider, showing blocks merge at epsilon=4
  
  RIGHT: Counterexample encoding (U=00,C=11,A=01,G=10)
    - Same blocks colored same way
    - Min distance between blocks: 2
    - Epsilon slider, showing blocks merge at epsilon=2

Big label: "Same genetic code. Different encoding. Different 'invariant'.
Therefore: min distance is NOT encoding-invariant."

Button: "Run live computation" pulls from an API route that actually computes
the Hamming distance on the browser using the selected encoding. Prove it's
not just static SVG.

### Vis 7 — Null Model A Labeling Bug (M6)
Two lines of code side by side:
  count_serine_unique  # <-- variable name
  if score["exactly_one_disconnected"]:  # <-- actual condition
      count_serine_unique += 1

Below: "Under null permutations, no codon is called 'Serine' — labels are
randomized. So this counter is measuring 'exactly one AA class has
Hamming-1-disconnection', not 'Serine specifically is the unique one'."

Then: Monte Carlo visualization — press "Run 10,000 permutations" button.
Shows histogram of "number of disconnected AA classes" for random codes with
matching degeneracy profile. Red line shows the observed value (1). Computes
p-value honestly.

### Vis 8 — KRAS/Fano Clinical Failure (B7, M4)
Embedded real data from MSK-IMPACT 2017 (extracted subset).
  - Bar chart of Fisher's exact test p-values per G12 variant
  - All bars at p=1.0
  - Annotation: "Zero enrichment of Fano-predicted amino acids. The clinical
    prediction track is dead."

### Vis 9 — The Bit-Position Bias (B4)
Two histograms stacked:
  TOP: "Observed bit-flip distribution across 61 reassignment events"
    - chi-sq p=0.006 under uniform null
  BOTTOM: "With position-specific transition/transversion weighting"
    - User can toggle weights (mitochondrial vs nuclear, Ts/Tv ratio 2 or 10)
    - Chi-sq p value updates live
    - If p rises above 0.05 under realistic weights, bar changes to "not
      significant"

Conclusion text: "The signal may be purifying selection at position 2
rediscovered in bit-coordinates."

### Vis 10 — The Salvageable Direction: Hypercube Coloring Optimality
Title: "What Would Actually Be Novel"
Interactive Monte Carlo sandbox:
  - Set degeneracy profile (locked to standard code by default)
  - Set physicochemical distance matrix (Grantham vs Miyata)
  - Press "Run 1,000 random colorings"
  - Histogram of edge-mismatch scores; red line = standard code's score
  - Quantile readout: "Standard code is at the X-th percentile"

This teaches the reader: THIS is what a legitimate optimality claim looks
like. The framework doesn't currently do this, but it could.

## DESIGN LANGUAGE

- Typography: Mantine's default (Inter or system UI) with mono for codons
- Color: Monochrome + one accent. Use Mantine's colorscheme switcher.
- Dark mode default
- Codons in mono typeface throughout: `UCU` not UCU
- Use MathJax/KaTeX (via better-react-mathjax or @matejmazur/react-katex) for:
  - GF(2)^6 notation
  - Gaussian binomial formula
  - Character homomorphism equations
- Badges: Red for BROKEN, yellow for SALVAGEABLE, green for SURVIVES
- Card/Paper containers for each claim (Mantine Paper with radius="md"
  shadow="sm")
- Collapsible "Model Consensus" strip under each claim with 12 avatar-style
  pills showing model names and their verdict emoji

## STATE MANAGEMENT

Zustand store shape:

interface StoreState {
  encoding: {
    C: [number, number];
    U: [number, number];
    A: [number, number];
    G: [number, number];
  };
  setEncoding: (newEncoding: Encoding) => void;
  activeClaimId: string;
  setActiveClaim: (id: string) => void;
  epsilonThreshold: number;
  setEpsilon: (e: number) => void;
}

Persist encoding to localStorage so the reader's encoding choice survives
section scrolling.

## ROUTING (TanStack Start)

File layout:
  /app/routes/__root.tsx      -> AppShell + navigation
  /app/routes/index.tsx       -> Landing + all claim sections
  /app/routes/api/compute.ts  -> Server route for Hamming computations
                                  (offload Monte Carlo if heavy)
  /app/routes/methods.tsx     -> How the analysis was done
  /app/routes/citations.tsx   -> Bibliography with clickable DOIs

TanStack server functions for:
  - POST /api/compute/hamming { encoding, codon_set_A, codon_set_B }
    returns min_distance, component_count
  - POST /api/compute/monte-carlo { n_samples, degeneracy_profile }
    returns histogram data
  - GET /api/data/msk-impact  returns parsed KRAS mutation data

## DATA

Include as static JSON (in /app/data/):
  - ncbi_25_tables.json (all 25 translation tables)
  - kras_msk_impact_2017.json (co-mutations)
  - reassignment_events.json (61 events with bit-flip data)
  - grantham_distance.json (20x20 AA distance matrix)

## WHAT I CAN PROVIDE YOU

- The Python code from my repository (codon_topo package)
- The 280 test cases and their expected values
- All NCBI translation tables as JSON
- The full 12-model consensus findings for "Model Consensus" strips
- The Kimi-k2.5 counterexample verified outputs

## MILESTONES

1. Bootstrap TanStack Start + Mantine + routing skeleton (don't write business
   logic, just get "Hello world" app running with nav)
2. Build Vis 1 (encoding picker) + Zustand store — this is the foundation
3. Build Vis 6 (Serine counterexample — this is THE centerpiece)
4. Build Vis 2 (Q_6 hypercube in react-three-fiber)
5. Build remaining Vis in order of importance
6. Polish: animations, dark mode, accessibility (keyboard nav, aria labels
   for each visualization)
7. Deploy to Vercel (TanStack Start has first-class Vercel support)

## WHAT NOT TO DO

- Don't use D3 selection/enter/exit APIs; use d3-force just for layout
  computation and render with React
- Don't use MUI or Radix; Mantine v8 has everything
- Don't use Redux; Zustand only
- Don't use SWR/React Query; TanStack Query v5 is included with TanStack Start
- Don't style inline; use Mantine's style props and CSS modules
- Don't use emoji everywhere; the design should feel serious/academic
- Don't memorize APIs — read docs for each library before writing

## REFERENCES TO READ BEFORE CODING

- https://tanstack.com/start/latest/docs/framework/react/overview
- https://mantine.dev/guides/typescript/
- https://zustand.docs.pmnd.rs/getting-started/introduction
- https://recharts.org/en-US/guide
- https://docs.pmnd.rs/react-three-fiber/getting-started/introduction

Start by building the skeleton (Milestone 1), show me the result, and we'll
iterate.
```

---

## Notes for when you run this

1. **Seed the new session with the data files**: After the prompt above, drop in the JSON data exports from `output/` (the NCBI tables, KRAS data, reassignment events). Let the new Claude pick them up in context.

2. **The key visual is Vis 6** (Serine counterexample side-by-side). If nothing else lands, that one slide tells the whole story. Get that working first after the skeleton.

3. **Performance note for Vis 2**: 64 vertices with 192 edges is trivial for three-fiber. Animations on epsilon slider should be butter-smooth. If they aren't, profile the re-renders — likely culprit is not memoizing the graph geometry.

4. **Accessibility**: The hypercube viz needs a "skip 3D; show me the data as a table" fallback for screen readers.

5. **Deployment**: TanStack Start deploys to Vercel as two output types (React and server). Configure `vercel.json` with the `tanstack-start` framework preset.
