# Testable Predictions: What Can Actually Be Done

**Purpose**: Enumerate concrete, testable hypotheses that could advance the Clayworth framework from descriptive redescription to a genuine scientific contribution. Organized by difficulty, data requirements, and publication potential.

---

## Tier 1: Can be done this week (before the April 14-16 meeting)

### Prediction 1.1 (CRITICAL — affects meeting agenda)
**Hypothesis**: Serine's UCN/AGY minimum Hamming distance varies across the 24 possible base-to-bit encodings. Specifically, there exists at least one encoding (kimi-k2.5's counterexample) where the minimum distance is 2, not 4.

**Method**: Extend `analysis/null_models.py` to emit, for each of the 24 encodings, the minimum inter-block Hamming distance for Serine (and for all other disconnected AAs).

**Expected**: Min distance takes values in {2, 4} across encodings (not invariantly 4).

**Outcome if verified**: The paper's strongest claim ("universal topological invariant at ε=4") must be weakened. This is already done at the counterexample level, but publishing the full distribution strengthens transparency.

**Time**: 1 hour.

**Status**: Counterexample verified. Full distribution pending.

---

### Prediction 1.2
**Hypothesis**: The bit-position bias observation (current chi-square p = 0.006) does NOT survive a position-weighted null model.

**Rationale**: Codon positions have very different transition/transversion rates and selection pressures. The uniform null is not biologically meaningful.

**Method**: 
1. Extract per-position substitution rates from within-lineage whole-genome comparisons for mitochondrial and nuclear genomes (e.g., Pagel-Meade Ts/Tv estimates).
2. Reweight the expected bit-flip distribution by these rates.
3. Redo chi-square test.

**Expected outcomes**:
- If weighted p < 0.01: Real signal beyond purifying selection. Publishable.
- If weighted p > 0.05: Signal is purifying selection rediscovered. Drop from paper.
- If 0.01 < weighted p < 0.05: Marginal — include with caveats.

**Time**: 3-4 hours.

---

### Prediction 1.3
**Hypothesis**: Organisms with codon-graph disconnections in variant codes have extra tRNA gene copies or tRNA innovations compared to their nearest non-variant relatives.

**Rationale**: Already mechanistically confirmed for yeast mitochondrial Thr (tRNA^Thr derived from tRNA^His, Su et al. 2011, PMC3113583). Test whether this pattern is general.

**Method**:
1. For each of 4 organisms with disconnections:
   - Saccharomyces cerevisiae (Table 3, Thr disconnect)
   - Scenedesmus obliquus (Table 22, Leu disconnect)
   - Pachysolen tannophilus (Table 26, Ala disconnect)
   - Candida albicans (Table 12, 3-component Ser)
2. Query GtRNAdb and mitotRNAdb for tRNA gene copies by anticodon.
3. Compare to nearest sequenced relative using standard code.

**Test**: Is there a statistically significant increase in tRNA gene copies for the reassigned amino acid in disconnection organisms vs. controls?

**Time**: 4-6 hours (data retrieval + analysis).

**Publication impact**: If positive, this is THE primary novel contribution.

---

## Tier 2: For Paper 1 (next 3 months)

### Prediction 2.1
**Hypothesis**: The presence of specific tRNA modifying enzymes (ADATs, Tit1, TrmA, TrmL, MnmE) predicts which codon reassignments are accessible in a lineage.

**Specific subhypotheses**:
- ADAT presence correlates with reassignments of codons decoded via inosine (ANN wobble)
- MnmE/MnmG presence correlates with stable reassignments in U34-modified codons
- TrmA absence correlates with mcm⁵U-dependent reassignments being unavailable

**Method**: 
1. For each NCBI variant code, identify the host organism(s).
2. Query InterPro/Pfam for presence/absence of modifying enzymes.
3. Build a contingency table: reassignment type × enzyme presence.
4. Logistic regression + Bayesian hierarchical model (adjust for phylogenetic non-independence).

**Data**: InterPro for enzyme families, KEGG for pathway annotations, NCBI for genome annotations.

**Publication impact**: High — this generates a predictive model for synthetic biology.

---

### Prediction 2.2
**Hypothesis**: The CUG clade shows predictable patterns of CUG reassignment stability based on decoding-network topology.

**Background**: The CUG clade (Candida, Pachysolen, related yeasts) has independently reassigned CUG to different amino acids (Ser, Ala, Leu). 

**Specific test**: Compare in-clade species with differing CUG decoding:
1. How many species show pure-CUG reassignment vs. ambiguous decoding?
2. Does the algebraic "neighborhood" of the new AA predict stability (reversion frequency)?
3. Do species with the "harder" reassignment require more compensatory tRNA changes?

**Method**: Meta-analysis of CUG clade genomes (40+ sequenced species as of 2026) for CUG decoding state, tRNA repertoire, compensatory modifications.

**Data**: Comparative genomics from Butler et al. (2009) and subsequent sequencing projects.

**Publication impact**: Medium — specific to one clade, but provides a clean test system.

---

### Prediction 2.3 (NEW DIRECTION — the hypercube coloring theorem)
**Hypothesis**: The standard genetic code represents a locally optimal coloring of the Q_6 hypercube that minimizes physicochemical mismatch across 1-Hamming edges, within the class-size constraint.

**Formal statement**: Let C be the standard genetic code. Let 𝒞 be the set of all colorings of Q_6 by 22 classes with class-size distribution matching C. For an objective function:

$$F(C) = \sum_{\substack{v,w \in Q_6 \\ d(v,w)=1}} \mathbb{1}[\text{color}(v) \neq \text{color}(w)] \cdot \Delta(\text{color}(v), \text{color}(w))$$

where Δ is Grantham (1974) physicochemical distance. The claim: F(C) is in the top X% (e.g., top 0.1%) of the F distribution over 𝒞.

**Method**:
1. Implement Grantham distance matrix (20×20 table, already tabulated).
2. Generate 10^6 random colorings with matching class sizes.
3. Compute F(C) for the standard code and F(c) for each random c.
4. Report the quantile.
5. Extend: is the code ALSO optimal conditioned on affine-subspace synonymous blocks?

**Expected**: Standard code falls in top 1-5% (consistent with Freeland 2000 findings).

**Publication impact**: HIGH — this is the gemini3/glm-5 proposed theorem that could anchor the paper.

**Time**: 2-4 weeks.

---

### Prediction 2.4
**Hypothesis**: Forward simulation of genetic code evolution under different models reveals which factors best reproduce the observed 25 NCBI tables.

**Models to compare**:
1. Random reassignment (baseline)
2. Wobble-weighted (Crick hypothesis only)
3. Algebra-preserving (Clayworth framework only)
4. Wobble + algebra combined
5. tRNA-network-constrained (requires tRNA data per genome)

**Method**: Agent-based simulation starting from the standard code. Apply reassignment events weighted by each model. Run 10^4 replicates per model. Compare distribution of simulated codes to observed NCBI tables using likelihood-ratio tests.

**Expected**: Model 5 (tRNA-network-constrained) wins. Model 3 (algebra alone) is insufficient. Model 4 may add marginal value.

**Publication impact**: Medium — decisively settles whether algebra adds predictive power beyond wobble.

---

### Prediction 2.5
**Hypothesis**: Codon reassignment events cluster in "ambiguous intermediate" states that are detectable in ribosome profiling data.

**Background**: Schultz & Yarus (1994) proposed the ambiguous intermediate theory. Modern ribosome profiling can detect real-time mistranslation.

**Method**: Mine published ribosome profiling data (GWIPS-viz, Trips-Viz) for organisms in active codon reassignment transitions. Look for:
- Dual-tRNA readout at specific codons
- Increased near-cognate pairing
- Frame-shifted proteins

**Data**: Public ribosome profiling datasets from Candida, Scenedesmus, and mitochondrial polysome profiles.

**Publication impact**: Medium — connects the framework to live biology.

---

## Tier 3: Paper 2 / grant horizon

### Prediction 3.1 (EXPERIMENTAL)
**Hypothesis**: In orthogonal translation systems, reassignments that violate the GF(2)^6 filtration structure show lower translation fidelity than structure-preserving reassignments, CONTROLLING for tRNA abundance.

**Critical test**: Distinguishes "algebra = wobble biology" from "algebra = independent predictor."

**Method**: Collaboration with Jason Chin (Cambridge) or Farren Isaacs (Yale):
1. Design 10 matched-pair reassignments:
   - 5 that preserve the filtration
   - 5 that violate it
   - Matched for tRNA abundance, mRNA folding, codon usage
2. Engineer into E. coli orthogonal translation system
3. Measure via:
   - Ribosome profiling (mistranslation rates)
   - Mass spectrometry (proteome-wide amino acid misincorporation)
   - Growth rate (fitness proxy)
4. Statistical test: paired t-test on fidelity metrics

**Prediction**: If algebra is epiphenomenal (just reflects wobble), no significant difference. If algebra has independent predictive power, filtration-preserving reassignments show higher fidelity.

**Outcome interpretations**:
- Significant difference: algebraic score has independent predictive power for synbio design
- No difference: algebra is purely descriptive of wobble biology

**Timeline**: 18-24 months with secured collaboration.

---

### Prediction 3.2
**Hypothesis**: Specific tRNA bridging (an engineered tRNA that reads both UCN and AGY) will be deleterious in E. coli or yeast due to mistranslation cascade.

**Rationale**: Paul's original synbio suggestion — test if the Serine disconnection is "functionally maintained."

**Method**:
1. Design a chimeric tRNA^Ser with dual-specificity anticodon modification
2. Express in E. coli under inducible promoter
3. Measure:
   - Growth rate (fitness)
   - Mistranslation rate (mass spec)
   - Proteotoxic stress markers (heat shock response)
4. Control: same amount of native tRNA^Ser expressed

**Prediction**: Fitness cost proportional to bridging efficiency and mistranslation rate.

**Key distinction**: Does the cost come from mistranslation alone, or from "breaking" the topological invariant independently?

**Timeline**: 12-18 months with a standard microbiology lab.

---

### Prediction 3.3
**Hypothesis**: The genetic code's hypercube coloring is evolutionarily "adaptive" in the sense that perturbations to nearby codings decrease fitness in simulated evolution.

**Method**: In-silico evolutionary rescue experiment:
1. Start with a perturbed code (swap one AA assignment)
2. Simulate forward evolution with selection on protein stability
3. Measure reversion rate and time-to-recovery
4. Compare to swap between distant codings

**Expected**: Swaps that violate the filtration structure require more compensatory mutations.

**Publication impact**: Medium — connects to evolutionary theory, but simulation-only.

---

### Prediction 3.4
**Hypothesis**: Codon reassignment events can be predicted from pre-existing tRNA innovations in ancestral lineages.

**Method**: Reconstruct ancestral tRNA repertoires for lineages preceding known reassignment events (e.g., pre-CUG clade ancestors). Test whether:
- tRNA gene duplications PRECEDE reassignment events
- Modification enzyme presence correlates with subsequent reassignment accessibility
- Reassignment direction is predictable from ancestral tRNA state

**Data**: Large-scale phylogenomics (400+ yeast genomes) + molecular clock calibration.

**Publication impact**: High if successful — connects framework to ancestral state reconstruction.

---

## Tier 4: Speculative / exploratory

### Prediction 4.1
**Hypothesis**: The affine-subspace structure of synonymous codon blocks reflects a specific tRNA-decoding optimization that can be formalized as a linear-algebra problem.

**Method**: Express wobble decoding as a linear map from anticodons to codons, find its kernel and image. Check if the 21 AA classes = orbits of this map.

**Publication impact**: If true — a deep mathematical characterization. If false — the affine structure is less meaningful than claimed.

---

### Prediction 4.2
**Hypothesis**: Ribosomes act as approximate syndrome decoders in the error-correcting-code sense.

**Method**: Formalize a parity-check matrix for the genetic code. Compare theoretical decoding thresholds to observed mistranslation rates in ribosome profiling.

**Publication impact**: Very high if demonstrable, but likely false — ribosomes are kinetic proofreaders, not algebraic decoders.

---

### Prediction 4.3
**Hypothesis**: The bit-position bias, when decomposed, distinguishes mutational spectrum from algebraic channeling.

**Method**: Build a linear model:
  bit-flip rate at position i = 
    β_1 * (position's Ts/Tv bias) + 
    β_2 * (position's purifying selection) + 
    β_3 * (position's "algebraic depth") + 
    ε

If β_3 is significant after controlling for β_1 and β_2, the algebraic channeling is real and independent.

**Publication impact**: High. Refines B4 from descriptive to potentially causal.

---

### Prediction 4.4 (NEGATIVE PREDICTION, worth stating)
**Hypothesis**: The framework predicts NOTHING meaningful about somatic cancer mutations.

**Test**: Already done (KRAS/Fano, p=1.0 across 6 G12 variants).

**Implication**: Include as a clean negative in the paper. Prevents future wasted research.

---

## Priority Matrix

```
                    HIGH IMPACT                    LOWER IMPACT

EASY (days-weeks)  P1.1 (counterexample)         P4.4 (KRAS negative - DONE)
                   P1.3 (tRNA duplications)        
                   P2.3 (hypercube coloring)     
                   P1.2 (weighted null)          

HARDER (months)    P2.1 (enzyme prediction)      P2.2 (CUG clade)
                   P2.4 (evolution sim)          P2.5 (ribosome profiling)

HARDEST (yrs)      P3.1 (experimental fidelity)  P3.2 (tRNA bridging)
                   P3.4 (ancestral tRNAs)        P3.3 (rescue sim)
                                                  P4.1-4.3 (speculative)
```

---

## Recommended First 4 Actions

```
╔═══════════════════════════════════════════════════════════════╗
║  BEFORE APRIL 14-16 MEETING (TIER 1):                         ║
║                                                                ║
║  1. Extend null_models.py to emit per-encoding min distance   ║
║     → Supports the honest revision of "M3 invariance" claim   ║
║                                                                ║
║  2. Implement position-weighted chi-square for bit-position   ║
║     bias                                                       ║
║     → Test survival of B4 under realistic null                ║
║                                                                ║
║  3. Query GtRNAdb for tRNA gene copies in the 4 disconnection ║
║     organisms                                                  ║
║     → First evidence for or against B3-upgraded tRNA hypothesis║
║                                                                ║
║  4. Implement Grantham distance + Monte Carlo skeleton for    ║
║     hypercube coloring optimality                             ║
║     → Seed the new direction if Paul accepts the pivot        ║
╚═══════════════════════════════════════════════════════════════╝
```

---

*This document aligns with the consensus of 10 LLM models (12 evaluation passes). Predictions are ranked by publication potential and difficulty, with direct references to data sources and collaborators where applicable.*
