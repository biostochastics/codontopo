# Clayworth Codon-Topo Framework: Claim-by-Claim Breakdown of What's Broken

**Author**: Sergey Kornilov
**Date**: April 13, 2026
**Purpose**: Detailed verbose analysis of every major claim, explaining specifically what is broken and how, based on multi-model adversarial review (10 distinct LLM models, 12 evaluation passes) and direct code verification.

---

## Summary Matrix

| Category | Claims Affected | Severity | Recoverable? |
|----------|-----------------|----------|--------------|
| Tautological (true by construction) | M1, M2, B1 | HIGH | Reframe only |
| False in generality (counterexamples exist) | M3, M5, M7, E1 | CRITICAL | Partial reframe |
| Mistest (wrong null or underpowered) | M6, B4, B6 | HIGH | Redo with proper methodology |
| Biologically dead | B7 | FATAL | No - must drop |
| Restatement (known biology renamed) | B2, B3 (partial), B5 | MODERATE | Reframe as survey |

---

## PART I: MATHEMATICAL CLAIMS

### M1. "Two-fold filtration = bit-5 difference"

**The original claim**

> For every amino acid with exactly 2 synonymous codons in the standard code, those codons differ at exactly bit position 5 (the last bit). Verified across all 25 NCBI translation tables with zero exceptions.

**What's broken: Tautological — true by construction, not a discovery**

The standard code's 2-synonymous codon families all end in either {U,C} (pyrimidine pair) or {A,G} (purine pair):

```
Phe: UUU, UUC        Tyr: UAU, UAC
His: CAU, CAC        Gln: CAA, CAG
Asn: AAU, AAC        Lys: AAA, AAG
Asp: GAU, GAC        Glu: GAA, GAG
Cys: UGU, UGC
```

Under the default encoding C=(0,0), U=(0,1), A=(1,0), G=(1,1):
- C and U differ only at the **second bit** (both pyrimidines → same first bit)
- A and G differ only at the **second bit** (both purines → same first bit)

**So by construction**: any two bases that are both pyrimidines or both purines differ only at their second bit. When these occupy the third codon position, the difference lands at global bit 5. The claim doesn't discover anything about the genetic code — it discovers a property of the encoding choice.

**Model consensus (3 pure-math models agree)**
- glm-5: "Forced by encoding choice, not discovered"
- kimi-k2.5: "Tautological — it is not a discovery about the code but a restatement of the encoding choice"
- minimax-m2.5: "Partially forced by construction"

**What would make this a discovery**

Finding a property that holds under the biological code but NOT under alternative encodings. Null Model C attempts this, but it shows only 16/24 encodings preserve the "bit-5" property — meaning the claim is encoding-dependent, confirming it's a property of the chosen encoding, not of the genetic code itself.

---

### M2. "Four-fold filtration = shared 4-bit prefix"

**The original claim**

> For every amino acid with exactly 4 synonymous codons, the codons share an identical 4-bit prefix and their 2-bit suffixes exhaust GF(2)² = {00, 01, 10, 11}.

**What's broken: Trivial under ANY bijection**

A 4-fold degenerate amino acid has codons {XYA, XYC, XYG, XYU} where X and Y are fixed. Under **any** bijection from {A,C,G,U} to GF(2)²:
- First two bases fixed → first 4 bits identical (shared prefix) **by definition**
- Third base runs through all 4 nucleotides → last 2 bits exhaust GF(2)² **by definition of bijection**

This isn't a claim about the standard encoding — it holds for any 2-bit encoding. The mathematical content is zero. It's literally the definition of "four-fold degenerate" restated in binary.

**Model consensus**
- glm-5: "Trivial by construction"
- kimi-k2.5: "Definitional — it is the definition of the set, not a derived property"
- minimax-m2.5: "Trivial identity, not a discovery"

**What would make this a discovery**

Identifying non-trivially-structured 4-fold blocks — e.g., if four-fold families violated this pattern. But they can't, because the biological definition of "4-fold degenerate" already implies this structure.

---

### M3. "Serine minimum Hamming distance = 4, invariant across all 24 encodings"

**The original claim**

> Serine has 6 codons: UCN (UCU, UCC, UCA, UCG) and AGY (AGU, AGC). Minimum inter-block Hamming distance = 4. Serine has 2 connected components at epsilon=1, reconnecting at epsilon=4. This holds across all 25 NCBI codes AND all 24 base-to-bit encodings.

**What's broken: The "minimum distance = 4" part is FALSE**

Kimi-k2.5 provided a counterexample which was independently verified by minimax-m2.7, mimo, gemini3, glm-5, and directly confirmed against the codebase.

**The counterexample encoding**: φ(U)=(0,0), φ(C)=(1,1), φ(A)=(0,1), φ(G)=(1,0)

Computed against the actual repository code:

```
Base pairwise distances:
  d(U, A) = 1   [was 2 in default encoding]
  d(C, G) = 1   [was 2 in default encoding]

UCN vectors: UCU=(0,0,1,1,0,0), UCC=(0,0,1,1,1,1),
             UCA=(0,0,1,1,0,1), UCG=(0,0,1,1,1,0)
AGY vectors: AGU=(0,1,1,0,0,0), AGC=(0,1,1,0,1,1)

Minimum UCN ↔ AGY Hamming distance = 2 (not 4)
```

**What actually survives**

- "Serine disconnected at epsilon=1 in all 24 encodings" → **TRUE** (min distance always ≥ 2 because UCN and AGY differ at 2 nucleotide positions, and any bijection preserves base-inequality — so at minimum each position contributes 1 to Hamming distance)
- "Minimum Hamming distance = 4" → **FALSE** in some encodings (can be 2)
- "Reconnects at epsilon=4" → **FALSE** for the counterexample encoding (reconnects at epsilon=2)

**Why this matters**

The "universal topological invariant at epsilon=4" language — the strongest single claim of the framework — must be weakened to "universal disconnection at epsilon=1 with encoding-dependent reconnection depth." That weaker statement is itself near-tautological (2 codon sets differing at 2 nucleotide positions trivially have Hamming ≥ 2 under any per-position encoding).

**Model consensus (independently verified by 5 models)**
- kimi-k2.5: Original counterexample
- minimax-m2.7: Verified with own calculations
- mimo: "Disproving the 'min distance 4' claim"
- gemini3: "I concur... Hamming distance of 2"
- glm-5: "I CONCUR... distance = 2"

---

### M4. "Fano lines: GGU XOR GUU XOR CAC = 0"

**The original claim**

> XOR triples in GF(2)⁶ are "Fano lines." GGU XOR GUU XOR CAC = 0 connects to KRAS G12V biology.

**What's broken: Three separate errors**

**Error 1 — Terminology**: The Fano plane is PG(2,2), the 7-point, 7-line projective plane over GF(2), living in a 3-dimensional ambient space with a specific incidence structure. GF(2)⁶ has dimension 6. XOR-zero triples in GF(2)⁶ are lines of the projective space PG(5,2), not Fano lines. These are completely different objects.

**Error 2 — Count**: The correct count of 2-dimensional subspaces in GF(2)⁶ (each giving one unordered XOR-zero triple) is the Gaussian binomial coefficient:

$$\binom{6}{2}_2 = \frac{(2^6 - 1)(2^6 - 2)}{(2^2 - 1)(2^2 - 2)} = \frac{63 \cdot 62}{3 \cdot 2} = 651$$

My original analysis said 2016 (wrong — that counts unordered pairs). 651 was independently computed by glm-5, kimi-k2.5, minimax-m2.5, minimax-m2.7, and mimo. With 651 such triples in a 64-element space, they aren't rare or special — they're generic linear dependencies.

**Error 3 — Biology**: Tested against cBioPortal MSK-IMPACT 2017 across all 6 KRAS G12 variants. Fisher's exact test with Bonferroni correction: all p = 1.0. Zero enrichment of "Fano partner" amino acids. No known biological mechanism would make XOR relationships in GF(2)⁶ constrain somatic mutation co-occurrence.

**Why this isn't fixable**

Renaming the triples (call them "affine lines" or "3-term linear relations") fixes Error 1. But:
- Error 2 (count) is mathematical fact
- Error 3 (biology) killed the clinical claim permanently

There's no reformulation that makes the clinical KRAS prediction work.

---

### M5. "Holomorphic embedding φ: GF(2)⁶ → ℂ³"

**The original claim**

> φ maps each 2-bit pair to a fourth root of unity: (0,0)→1, (0,1)→i, (1,0)→−1, (1,1)→−i. This is a holomorphic embedding extending a character of GF(8)*.

**What's broken: Three separate errors**

**Error 1 — Not holomorphic**: "Holomorphic" requires a map from an open subset of ℂⁿ to ℂᵐ that is complex-differentiable. GF(2)⁶ is a **finite discrete set** (64 points). No classical sense in which a map from a finite set can be holomorphic. Category error.

**Error 2 — Not a character**: Paul's defense reframes as "a character of GF(8)* extended to GF(2)⁶." This fails.

Proof (multiple models):

In (GF(2)², +), every element has order ≤ 2. Characters χ: (GF(2)², +) → ℂ* must satisfy χ(v+v) = χ(v)² = χ(0) = 1. Since the only square roots of 1 in ℂ* are ±1, characters can only take values in {+1, −1}.

But the framework's map sends (0,1) → i, where i has order 4. Check:
- χ((0,1) + (0,1)) = χ(0,0) = 1
- χ((0,1))² = i² = −1

These don't match. Not a homomorphism.

No natural group structure on GF(2)² (additive, multiplicative via GF(4)*, etc.) makes this map a homomorphism.

**Error 3 — GF(8) confusion**: GF(8)* is cyclic of order 7. Fourth roots of unity live in cyclic order 4. gcd(7, 4) = 1, so no non-trivial homomorphism. The claim "GF(2)⁶ = GF(8) × GF(8)" is wrong as a field equality (GF(8) × GF(8) isn't even a field — not a direct product of fields). As additive groups, (GF(2))⁶ ≅ (GF(2))³ × (GF(2))³ ≅ (GF(8), +) × (GF(8), +), but this doesn't rescue the character claim.

**What it actually is**

A coordinate-wise bijection from GF(2)² (as a set) to the fourth roots of unity (as a set). Call it an "embedding" in the set-theoretic sense, nothing more.

**Model consensus**

4 pure-math models independently verified the non-character argument with identical reasoning.

---

### M6. "Null Model A: Serine uniqueness p < 0.05"

**The original claim**

> Null Model A runs 100,000 random codon-to-amino-acid assignments preserving the degeneracy structure. In < 5% of these does Serine appear as the unique disconnected amino acid.

**What's broken: Statistic-claim mismatch**

Under the null, there IS no Serine — labels are scrambled. The code (`null_models.py` line 133) actually computes:

```python
if score["exactly_one_disconnected"]:
    count_serine_unique += 1
```

Variable name: `count_serine_unique`. Actual condition: "exactly one amino acid class (of any identity) has its codons disconnected at Hamming distance 1."

**The real statistic**: "probability that random codes with matching degeneracy profile have exactly one disconnected AA class"

**The claimed statistic**: "probability that Serine specifically is uniquely disconnected"

These are different statements. Under null permutations, the concept of "Serine" doesn't even exist — its identity has been shuffled into anonymous classes.

**Why it matters**

Any peer reviewer will catch this immediately. The paper's headline claim — "Serine's disconnection is uniquely improbable under random codes" — collapses to the weaker, less interesting statement.

**Model consensus**

gpt-5.2-pro flagged this directly from code inspection.

**How to fix**

Relabel the variable. Rewrite the claim: "The standard code exhibits exactly one disconnected-at-ε=1 amino acid class (Serine); under block-size-preserving random permutations, exactly-one-disconnected occurs in < X% of samples."

That's a legitimate observation, just much weaker than "Serine is special."

---

### M7. "PSL(2,7) as the fundamental symmetry group of the genetic code"

**The original claim** (mostly in Paul's broader Clayworth Algebra framing):

> The algebraic engine has a constraint kernel built on the Fano plane together with PSL(2,7).

**What's broken: Already tested and rejected in peer-reviewed literature**

Antoneli & Forger (2011, Mathematical and Computer Modelling 53:1469–1488) conducted an exhaustive search for finite simple groups whose irreducible representations could accommodate the 64-codon structure. PSL(2,7) (order 168) was examined and explicitly rejected — it has no 64-dimensional irreducible representation (its irreps have dimensions 1, 3, 6, 7, 8).

Paper's conclusion: "None of the groups can reproduce the multiplet distribution of the genetic code without freezing in the last step."

**Weaker reformulation**: "PSL(2,7) acts on the GF(2)³ substructure whose XOR triples generate lines." This is trivially true (PSL(2,7) ≅ GL(3,2) is the automorphism group of PG(2,2) by definition), but provides no biological content that survives Antoneli & Forger's rejection.

**Why it matters**

Any peer reviewer familiar with the algebraic genetic code literature will know Antoneli & Forger and flag this claim. Needs to be dropped or explicitly addressed with a new differentiating argument.

---

### E1. "Bit 2 of the encoding = weak/strong hydrogen bonding"

**The original claim** (from TN-2026-11 and my first analysis):

> Bit 1 distinguishes pyrimidines from purines. Bit 2 distinguishes weak from strong hydrogen bonding.

**What's broken: Bit 2 actually groups amino vs keto, NOT weak/strong**

Check the default encoding carefully:

| Base | Encoding (bit1, bit2) | Pyr/Pur | Amino/Keto | H-bonds |
|------|----------------------|---------|------------|---------|
| C | (0,0) | Pyrimidine | Amino | 3 (Strong) |
| U | (0,1) | Pyrimidine | Keto | 2 (Weak) |
| A | (1,0) | Purine | Amino | 2 (Weak) |
| G | (1,1) | Purine | Keto | 3 (Strong) |

Bit 2 groupings:
- Bit 2 = 0: {C, A} — both **amino** (have NH₂ groups at position 6/4)
- Bit 2 = 1: {U, G} — both **keto** (have =O groups at position 4/6)

Weak/strong H-bonding would require:
- Weak (2 bonds): {A, U}
- Strong (3 bonds): {C, G}

Under default encoding:
- {A, U}: bit 2 = 0 and 1 → not the same — NOT captured by bit 2
- {C, G}: bit 2 = 0 and 1 → not the same — NOT captured by bit 2

**Gemini3's elegant insight**: Weak/strong H-bonding IS recoverable as **bit 1 XOR bit 2** (a parity function):
- C=(0,0): 0⊕0 = 0 → Strong ✓
- U=(0,1): 0⊕1 = 1 → Weak ✓
- A=(1,0): 1⊕0 = 1 → Weak ✓
- G=(1,1): 1⊕1 = 0 → Strong ✓

So H-bonding is a **derived** linear functional, not primary.

**Why it matters**

Saying "we encode pyrimidine/purine and weak/strong hydrogen bonding in 2 bits" sounds more biologically meaningful than amino/keto. The weak/strong version has stronger thermodynamic implications (base-pair stability). The amino/keto version is a tautomer classification, less central to function. When writing the encoding motivation, you must state what bit 2 actually is, not what sounds best.

**Model agreement (4 vs 1)**

- gpt-5.2-pro: bit 2 = amino/keto
- minimax-m2.7: bit 2 = amino/keto (confirmed)
- mimo: bit 2 = amino/keto (confirmed — called it a "biochemical misattribution")
- gemini3: bit 2 = amino/keto; H-bonding = XOR (elegant insight)
- glm-5: INCORRECTLY said bit 2 = weak/strong

---

## PART II: BIOLOGICAL CLAIMS

### B1. "Degeneracy structure reflects algebraic order"

**The original claim**

> The 2-fold and 4-fold patterns of synonymous codons reveal a hidden algebraic structure non-randomly embedded in the code.

**What's broken: Redescription, not explanation**

Claims M1 and M2 are tautological (forced by encoding). So whatever "algebraic structure" they reveal is a property of the representation, not the genetic code. The 2-fold/4-fold patterns are fully explained by:

1. **Crick's wobble hypothesis (1966)**: Third codon position is least constrained by tRNA-mRNA pairing. Non-Watson-Crick base pairing (wobble) at position 34 of tRNA anticodon accommodates multiple codons.

2. **tRNA isoacceptor families**: Multiple tRNAs with different anticodons read the same amino acid. Example: E. coli has 6 different tRNA^Leu genes with different anticodons.

3. **Anticodon modifications**: Inosine (I) at position 34, introduced by ADAT enzymes, reads codons ending in U, C, or A via non-canonical pairing. 5-methylaminomethyl-2-thiouridine (mnm⁵s²U) restricts reading to A and G.

The framework adds no new mechanism. It notates existing patterns in a new coordinate system.

**Model consensus**: gpt-5.4-pro, gemini3, claude-opus all converge.

---

### B2. "Serine's disconnection is a universal invariant across 25 codes"

**The original claim**

> The fact that Serine's UCN/AGY split persists across every known genetic code variant is a deep topological invariant.

**What's broken: Phylogenetic inertia, not mathematical invariance**

All 25 NCBI translation tables derive from the standard code via a small number of codon reassignments. None of these reassignments happen to bridge the UCN/AGY gap (which would require tandem 2-nucleotide changes, not single-codon reassignments). So the "universal persistence" reflects **what evolution hasn't done**, not an algebraic law forbidding it.

Additionally, the biological fact (Serine has two disjoint codon families requiring tandem substitution) is textbook:
- Schwartz et al. (2019, Sci Rep 9:17238): "Serine is the only amino acid encoded by two disjoint codon sets (TCN and AGY) so that a tandem substitution of two nucleotides is required to switch between the two sets."
- Rogozin et al. (2016, PNAS 113:13109): "TCN↔AGY switches are driven by selection, requiring tandem double substitutions."
- Inouye & Inouye (2020, PNAS 117:28572): Evidence from E. coli for AGY as primordial.
- Bernhardt (2016, Life 6:10): Class II tRNA structure and split codons as evolutionary fossils.

This fact predates the framework by decades.

**Model consensus**

- gpt-5.4-pro: "Phylogenetic inertia, not invariance"
- gemini3: "Simple phylogenetic inertia — all 25 codes are minor variants of the standard code"
- claude-opus: "Tautological given the evolutionary derivation"

---

### B3. "Novel disconnections in variant codes (Thr, Leu, Ala, 3-component Ser)"

**The original claim**

> In variant codes (yeast mitochondrial Thr, chlorophycean Leu, Pachysolen Ala, Candida 3-component Ser), new topological disconnections emerge, revealing structural signatures of evolution.

**What's partly salvageable — THE MOST PROMISING ORIGINAL CLAIM**

The underlying biological events are well-documented:
- Yeast mito: CUN block reassigned Leu → Thr (Miranda et al. 2006, PMC3113583)
- Chlorophycean mito: UAG → Leu
- Pachysolen: CUG → Ala
- Candida (CUG clade): CUG → Ser

Calling these "topological disconnections" is geometric language for "the codon set after reassignment is split in Hamming space." Not a new biological discovery per se.

**What could make this publishable**: Claude-opus proposed the test — do these disconnections correlate with tRNA gene duplications? Tavily search (April 13, 2026) confirmed this is mechanistically real for at least one case: yeast mitochondrial Thr acquired a new tRNA^Thr **derived from tRNA^His** (Su, Ho, Wang, Chang 2011, PMC3113583). If this pattern holds across all four disconnected cases, THAT is a publishable finding.

**Status**: Salvageable if the tRNA-duplication correlation is demonstrated across the full disconnection catalogue.

---

### B4. "Bit-position bias in reassignments, chi-square p = 0.006"

**The original claim**

> Across 61 codon reassignment events, bit position 4 (first wobble bit) accumulates 13/35 bit-flips while position 2 has zero. Chi-square test under uniform null gives p = 0.006.

**What's broken: The null hypothesis is invalid**

The uniform-across-6-bit-positions null model is not what any biologist believes. The three codon positions have vastly different mutational and selective dynamics:

- **Position 1**: Second-strongest selection; biophysical constraint on amino acid properties
- **Position 2**: STRONGEST selection; determines amino acid class (hydrophobic/polar/charged)
- **Position 3 (wobble)**: Weakest selection; tolerates transitions freely; many synonymous changes

A proper null would weight bit-flip probabilities by:
1. Position-specific transition/transversion rates (observable from within-lineage substitutions)
2. Mitochondrial vs. nuclear substitution biases (mitos have strand-asymmetric mutation, e.g. APOBEC3 / 8-oxo-G)
3. Phylogenetic non-independence (reassignments aren't independent events)
4. GC content biases
5. Codon usage frequency in source genomes

Under the current uniform null, finding p = 0.006 is roughly equivalent to rediscovering "second codon positions are under strongest purifying selection" — textbook biology (Woese 1965, Knight et al. 2001).

**What could save it**: Redo the chi-square with position-specific substitution rates weighting expected values. If p < 0.01 survives, it represents algebraic channeling beyond the known purifying selection gradient. If not, the signal is mutational spectrum.

**Status**: Needs redoing with proper null.

**Model consensus**
- claude-opus flagged this explicitly with detailed methodology
- gemini3: "The math merely rediscovers basic purifying selection"
- gpt-5.4-pro: "Not trivial in the sense of 'false,' but not novel unless shown to survive in biologically native coordinates"

---

### B5. "83% of 1,280 single-codon reassignments preserve the algebraic structure"

**The original claim**

> 83% of possible single-codon reassignments preserve the filtration and Serine disconnection, indicating a large viable design space for synthetic biology.

**What's broken: Measures encoding artifact, not biological meaning**

Given that:
- Two-fold bit-5 pattern is forced by encoding (M1 tautology)
- Four-fold prefix pattern is trivial under any bijection (M2)
- Serine disconnection persists because UCN/AGY differ at 2 positions (M3 partial-tautology)

...when you compute how many single-codon reassignments "preserve" these properties, you're largely measuring which reassignments leave the definitional structure intact. The 83% number has no clear biological meaning.

**What would give it meaning**: Weighting by actual biological costs:
- tRNA gene duplication requirements
- Aminoacyl-tRNA synthetase (aaRS) retargeting costs
- Release factor (RF1/RF2) modifications
- Proteome mistranslation burden
- Fitness cost in engineered strain construction

A synthetic biologist at Jason Chin's (Cambridge) or Farren Isaacs's (Yale) lab wouldn't use this score because it doesn't correlate with the physical difficulty of implementing the reassignment.

**Model consensus**

- gemini3: "Not viable. Synthetic biologists need tRNA abundance, ribosome affinity, and mRNA folding — not algebraic preservation"
- gpt-5.4-pro: "Little mechanistic meaning because the 'structure' is partly encoding-imposed"
- claude-opus: "Likely reflects that most single-codon reassignments change a codon to a nearby AA, which is already explained by tRNA misreading mechanisms"

---

### B6. "Depth calibration: ε vs evolutionary age shows ρ = 0.0 (n=6)"

**The original claim**

> Spearman correlation between reconnection epsilon and evolutionary age is zero. Therefore, topology reflects structural position in GF(2)⁶, not time — the genuinely interesting finding.

**What's broken: n = 6 is statistically meaningless**

With 6 data points, Spearman's ρ has almost no power. The possible ρ range under the null is enormous. Getting ρ = 0.0 tells you essentially nothing — you couldn't have detected a positive result either.

Additionally, the underlying hypothesis — "older events produce deeper disconnections" — has no strong biological motivation. Codon reassignments are episodic events driven by:
- Genome reduction pressure (mitochondria, endosymbionts)
- tRNA loss / duplication
- Release factor innovations
- Ambiguous intermediate decoding states

No continuous molecular-clock-like process predicts older events produce deeper disconnections. The "correlation being zero" isn't surprising or informative.

**Status**: Cannot be salvaged. It's a non-result, not a negative result. Drop from the paper.

**Model consensus**
- gpt-5.4-pro: "n=6 is a statistical non-result"
- claude-opus: "Statistically meaningless"
- gemini3: "Statistically meaningless. Codon reassignment is episodic, not a continuous molecular clock"

---

### B7. "KRAS G12V Fano-line co-occurrence predicts clinical mutation patterns"

**The original claim**

> GGU XOR GUU XOR CAC = 0 in GF(2)⁶ predicts CAC/His should be enriched at KRAS G12V co-mutation sites, giving a liquid biopsy interpretation tool.

**What's broken: Tested against real data, result was zero enrichment**

Against MSK-IMPACT 2017 via cBioPortal (1,670 KRAS mutations), Fisher's exact test with Bonferroni correction returned p = 1.0 for every G12 variant (G12A, G12C, G12D, G12R, G12S, G12V). Zero enrichment.

**Why this was inevitable**

No biological mechanism makes XOR relationships in GF(2)⁶ constrain somatic mutation co-occurrence. Somatic mutations are driven by:
- DNA replication errors (low fidelity of specific polymerases)
- UV, ionizing radiation, chemical mutagens (cigarette smoke, aflatoxin, etc.)
- APOBEC3 cytidine deaminases (characteristic TCW → TKW signatures)
- DNA repair failures (MMR deficiency, BRCA deficiency)
- Positive/negative selection on protein function

None of these "compute" XOR or care about GF(2)⁶ structure. The framework asked DNA polymerase to do abstract algebra. It does not.

**Why the framework can't recover this**

The clinical application track is dead. There's no reformulation of XOR in GF(2)⁶ that would plausibly predict somatic co-mutation patterns. Any attempt to rescue (e.g., "the effect is subtle and masked by mutational signatures") would be unfalsifiable.

**Status**: Irrecoverable. Must be dropped from the paper entirely. The only useful thing is to report the clean negative so future researchers don't waste time on it.

---

## PART III: THE OVERALL PATTERN

```
┌────────────────────────────────────────────────────────────────┐
│  CATEGORY           │  CLAIMS          │  WHAT'S WRONG         │
├────────────────────────────────────────────────────────────────┤
│                     │                  │  Forced by encoding    │
│  Tautological       │  M1, M2, B1     │  definition, not by    │
│                     │                  │  biology               │
├────────────────────────────────────────────────────────────────┤
│                     │                  │  Counterexamples exist,│
│  False in           │  M3, M5, M7, E1 │  terminology wrong, or │
│  generality         │                  │  already refuted       │
├────────────────────────────────────────────────────────────────┤
│                     │                  │  Wrong null model,     │
│  Mistest            │  M6, B4, B6     │  wrong statistic, or   │
│                     │                  │  underpowered          │
├────────────────────────────────────────────────────────────────┤
│                     │                  │  Tested against real   │
│  Biologically       │  B7             │  data, failed cleanly  │
│  dead               │                  │                        │
├────────────────────────────────────────────────────────────────┤
│                     │                  │  Rediscovery of known  │
│  Restatement        │  B2, B3 (part), │  biology in new        │
│                     │  B5             │  notation              │
└────────────────────────────────────────────────────────────────┘
```

---

## What Actually Survives

| Item | Status | Path to Publication |
|------|--------|---------------------|
| Systematic NCBI 25-table survey | VALID | Computational biology methods paper |
| Disconnection catalogue (B3) | SALVAGEABLE | Link to tRNA duplication (confirmed for yeast Thr) |
| Bit-position bias (B4) | SALVAGEABLE | Redo with position-weighted null |
| KRAS/Fano clean negative (B7) | VALID AS NEGATIVE | Report to prevent future false positives |
| Open-source pipeline | VALID | Tool / software paper possible |
| Hypercube coloring optimality theorem | NEW DIRECTION | Gemini3/glm-5 proposal, bridges to Freeland/Koonin |

Everything else is either tautological, false, mistested, dead, or redundant.

---

*This report is based on 12 independent evaluations across 10 distinct LLM models (gpt-5.2-pro, gpt-5.4-pro, glm-5, kimi-k2.5, minimax-m2.5, minimax-m2.7, mimo, deepseek-r1, gemini3, claude-opus), direct code verification of the kimi-k2.5 counterexample against this repository's `hamming_distance()` function, and confirmed against peer-reviewed literature via Tavily searches (April 13, 2026).*
