// Cover letter for the revised submission to BioSystems.
// Manuscript: BIOSYS-D-26-00689
// Render with: typst compile output/cover_letter.typ output/cover_letter.pdf

#set page(
  paper: "us-letter",
  margin: (x: 1in, y: 1in),
)
#set text(font: "New Computer Modern", size: 11pt, lang: "en")
#set par(justify: true, leading: 0.65em, first-line-indent: 0pt, spacing: 1.1em)

// ---------- Sender block ----------
#align(right)[
  #text(weight: "bold")[Sergey Kornilov, PhD] \
  Biostochastics, LLC \
  Seattle, WA, USA \
  #link("mailto:sergey@biostochastics.com")[sergey\@biostochastics.com]
]

#v(1em)

// ---------- Date ----------
July 8, 2026

#v(0.5em)

// ---------- Recipient ----------
Dr. Robert Prinz \
Editor \
#emph[BioSystems]

#v(0.8em)

// ---------- Re-line ----------
*Re: Revised manuscript BIOSYS-D-26-00689*  \
*"Robust error-minimization in the genetic code across physicochemical metrics and variant codes: a graph-theoretic analysis in GF(2)#super[6]"*

#v(0.4em)

// ---------- Salutation ----------
Dear Dr. Prinz,

We are pleased to submit a revised version of the above manuscript by Paul Clayworth and Sergey Kornilov for further consideration in #emph[BioSystems]. This letter summarises the changes made in response to all three reviewer reports; the accompanying #emph[Response to Reviewers] provides the detailed point-by-point reply. All new analyses are reproducible from the public `codon-topo` pipeline at the commit tagged in the manuscript reproducibility stamp.

*Reviewer 1* asked us to engage with recent BioSystems work by Yurova Axelsson and Khrennikov on 2-adic and dynamical-system representations of the genetic code. We have added a Discussion subsection (§4.5, "Correspondence with 2-adic codon algebra") that situates our finite Hamming-graph framework alongside the 2-adic attractor picture, together with related p-adic work by Khrennikov and Kozyrev and by Dragovich and Mišić. As a concrete finite bridge between the two perspectives, we have added a Walsh--Hadamard / 2-adic spectral probe (Supplement §S16) with two null models and two encoding-invariant Walsh signatures of the wobble structure --- spectral depth 544 and label-spectrum fraction $S = 0.7514$ --- both interpreted explicitly as algebraic invariants of wobble-box organisation, not as evidence of beyond-wobble optimisation.

*Reviewer 2* raised three methodological points, all now addressed. First, the reviewer correctly noted that our original four physicochemical measures share substantial polarity/volume content; we have removed all "disjoint conceptual bases" and "independent replication" language and reframe the panel as a #emph[convergent sensitivity envelope] of established but partially overlapping code-independent parameterisations, with pairwise and partial Spearman correlations quantified in Results §3.1 and Supplement §S15. Second, we added ProtSub (Jia & Jernigan 2021) as a structure-aware fifth measure and report a full Freeland--Hurst null Monte Carlo ($n = 10000$, $p = 0.0004$; Table 1); the five-metric panel passes Bonferroni at $alpha = 0.05$. Third, we now cite Di Giulio (2001) in Methods §2.2 as the explicit rationale for restricting our primary panel to code-independent scales, and treat ProtSub as a sensitivity check rather than primary evidence. Kosiol et al. (2004) is now cited in the metric-overlap discussion.

We have also added a short Discussion paragraph (§4.6) placing our topology-avoidance result alongside the very recent Tsour, Slavov et al. (2026, #emph[Nature]) alternate-decoding dataset. At the event level, 65.0 % of their 5,873 high-confidence sense-codon recoding events involve source--target codon pairs at nucleotide-Hamming distance 1 versus a 39.5 % baseline (binomial $p < 10^(-100)$), providing convergent event-level context for the biological relevance of the local codon neighbourhoods that our null models test.

*Reviewer 3* recommended acceptance and asked us to discuss --- and consider as an alternative --- the non-power model of the genetic code introduced by Gonzalez, Giannerini and Rosa (2016, #emph[Philos Trans R Soc A] 374(2063):20150062), which uses a redundant signed-binary numeration with Fibonacci-like positional weights $[1,1,2,4,7,8]$ so that codon-to-integer maps become intentionally many-to-one, mirroring code degeneracy. We have added a paragraph to Discussion §4.5 ("Broader connections") that cites Gonzalez et al., describes the non-power construction, and makes the point that our primary topology-avoidance result is defined on the encoding-independent Hamming graph $H(3,4)$ --- which depends only on which nucleotides differ between two codons, not on any codon-level bijection at all, whether power or non-power --- so the paper's principal topology claim is not exposed to this concern. A non-power extension of the $Q_6$ subgraph analysis (which is already encoding-dependent within the power family, §3.4 and Supplement §S4) is explicitly flagged as a natural continuation, but not attempted here because it requires reformulating the underlying mutation graph when codon-to-integer maps become many-to-one. The bibliography now includes the Gonzalez, Giannerini and Rosa (2016) entry.

The revised package accordingly comprises: the revised manuscript with the new §4.5 and §4.6 subsections and a five-metric Results §3.1; a substantially expanded supplement (new §S15 metric correlations, §S16 Walsh--Hadamard 2-adic probe, §S17 Slavov SAAP cross-analysis); an updated bibliography with nine new references; and the accompanying point-by-point Response to Reviewers.

All analyses remain fully reproducible through the open-source `codon-topo` pipeline (`codon-topo all --output-dir=./output --seed=135325`). This revision is original, has not been published previously, and is not under consideration elsewhere. All authors have approved the submission.

We thank the editor and the reviewers for their careful reading and constructive comments, and we are happy to address any further points raised in re-review.

#v(1em)

Sincerely,

#v(2.5em)

Sergey Kornilov, PhD \
Corresponding author \
Biostochastics, LLC \
Seattle, WA, USA \
#link("mailto:sergey@biostochastics.com")[sergey\@biostochastics.com]

#v(0.4em)

#emph[On behalf of all authors.]
