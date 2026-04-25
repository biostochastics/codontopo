"""Synthetic biology feasibility scoring for alternative genetic codes.

Scores each possible single-codon reassignment from the standard code
by how well it preserves the algebraic structure discovered in WS1:
two-fold filtration, four-fold filtration, Serine disconnection, and
overall disconnection profile.

The feasibility score is a weighted composite:
- Two-fold filtration intact: 0.25
- Four-fold filtration intact: 0.25
- Serine disconnected: 0.30
- Single disconnected AA (like standard): 0.20
"""

from codon_topo.core.encoding import ALL_CODONS
from codon_topo.core.genetic_codes import STANDARD
from codon_topo.core.filtration import check_twofold, check_fourfold
from codon_topo.core.homology import disconnection_catalogue


def score_variant_code(code: dict[str, str]) -> dict:
    """Score a genetic code variant for algebraic structure preservation.

    Returns dict with boolean flags and a composite feasibility score [0, 1].
    """
    # Two-fold check
    tw = check_twofold(code)
    tw_intact = all(ok for _, ok, _ in tw) if tw else True

    # Four-fold check
    ff = check_fourfold(code)
    ff_intact = all(ok for _, ok in ff) if ff else True

    # Disconnection analysis
    cat = disconnection_catalogue(code)
    disconnected_aas = [e["aa"] for e in cat]
    ser_disc = "Ser" in disconnected_aas
    n_disc = len(disconnected_aas)

    # Composite score
    score = 0.0
    if tw_intact:
        score += 0.25
    if ff_intact:
        score += 0.25
    if ser_disc:
        score += 0.30
    if n_disc == 1:  # Exactly one disconnected (like standard)
        score += 0.20
    elif n_disc == 0:
        score += 0.10  # Partial credit

    return {
        "twofold_intact": tw_intact,
        "fourfold_intact": ff_intact,
        "serine_disconnected": ser_disc,
        "n_disconnected_aas": n_disc,
        "disconnected_aas": disconnected_aas,
        "feasibility_score": round(score, 4),
    }


def single_reassignment_landscape(
    base_code: dict[str, str] | None = None,
) -> list[dict]:
    """Score every possible single-codon reassignment from the base code.

    For each of the 64 codons, try reassigning it to every other AA
    (including Stop) and score the result.

    Excludes identity reassignments (same AA).
    """
    ref = base_code or STANDARD
    all_aas = sorted(set(ref.values()))  # Includes 'Stop'

    landscape: list[dict] = []
    for codon in ALL_CODONS:
        original_aa = ref[codon]
        for new_aa in all_aas:
            if new_aa == original_aa:
                continue
            variant = dict(ref)
            variant[codon] = new_aa
            sc = score_variant_code(variant)
            landscape.append(
                {
                    "codon": codon,
                    "original_aa": original_aa,
                    "new_aa": new_aa,
                    "feasibility_score": sc["feasibility_score"],
                    "twofold_intact": sc["twofold_intact"],
                    "fourfold_intact": sc["fourfold_intact"],
                    "serine_disconnected": sc["serine_disconnected"],
                    "n_disconnected_aas": sc["n_disconnected_aas"],
                }
            )
    return landscape


def feasibility_summary(
    base_code: dict[str, str] | None = None,
) -> dict:
    """Summary statistics for the single-reassignment landscape.

    Classifies variants into high (>=0.8), medium (0.5-0.8), low (<0.5)
    feasibility tiers.
    """
    landscape = single_reassignment_landscape(base_code)
    high = [e for e in landscape if e["feasibility_score"] >= 0.8]
    medium = [e for e in landscape if 0.5 <= e["feasibility_score"] < 0.8]
    low = [e for e in landscape if e["feasibility_score"] < 0.5]

    best = sorted(landscape, key=lambda x: x["feasibility_score"], reverse=True)[:10]

    return {
        "total_variants": len(landscape),
        "high_feasibility": len(high),
        "medium_feasibility": len(medium),
        "low_feasibility": len(low),
        "best_variants": best,
    }


def topology_avoidance_test() -> dict:
    """Test whether natural reassignments avoid creating new disconnections.

    Compares the rate of disconnection-creating changes among OBSERVED
    natural reassignments vs ALL POSSIBLE single-codon reassignments.

    Prediction: observed events are depleted for topology-breaking changes
    (i.e., observed rate < possible rate).
    """
    from scipy.stats import fisher_exact, hypergeom

    from codon_topo.analysis.reassignment_db import build_reassignment_db

    # 1. What fraction of POSSIBLE reassignments create new disconnections?
    # (Use the standard code as base, matching the original landscape)
    landscape = single_reassignment_landscape(STANDARD)
    standard_cat = disconnection_catalogue(STANDARD)
    standard_disc_aas = {e["aa"] for e in standard_cat}

    possible_creates_disc = 0
    possible_no_disc = 0
    for entry in landscape:
        # Does this variant have MORE disconnected AAs than standard?
        new_disc_aas = set()
        variant = dict(STANDARD)
        variant[entry["codon"]] = entry["new_aa"]
        new_cat = disconnection_catalogue(variant)
        new_disc_aas = {e["aa"] for e in new_cat}
        novel_discs = new_disc_aas - standard_disc_aas
        if novel_discs:
            possible_creates_disc += 1
        else:
            possible_no_disc += 1

    # 2. What fraction of OBSERVED natural reassignments create disconnections?
    db = build_reassignment_db()
    observed_creates_disc = 0
    observed_no_disc = 0
    seen: set[tuple[str, str]] = set()  # de-duplicate
    for e in db:
        key = (e.codon, e.target_aa)
        if key in seen:
            continue
        seen.add(key)
        variant = dict(STANDARD)
        variant[e.codon] = e.target_aa
        new_cat = disconnection_catalogue(variant)
        new_disc_aas = {entry["aa"] for entry in new_cat}
        novel_discs = new_disc_aas - standard_disc_aas
        if novel_discs:
            observed_creates_disc += 1
        else:
            observed_no_disc += 1

    a = observed_creates_disc
    n_obs = observed_creates_disc + observed_no_disc
    N = possible_creates_disc + possible_no_disc
    K = possible_creates_disc

    rate_observed = a / max(n_obs, 1)
    rate_possible = K / max(N, 1)

    # Hypergeometric test: n_obs draws from landscape of N items,
    # K of which are "disc-creating". P(X <= a) tests depletion.
    hypergeom_p = float(hypergeom.cdf(a, N, K, n_obs))

    # Also report Fisher's exact for comparison.
    # Observed is a subset of possible, so the second row is the complement.
    b = observed_no_disc
    c_complement = K - a  # possible disc-creating NOT in observed
    d_complement = (N - K) - b  # possible non-disc NOT in observed
    fr = fisher_exact(
        [[a, b], [max(c_complement, 0), max(d_complement, 0)]],
        alternative="less",
    )
    odds_ratio = float(fr.statistic)  # type: ignore[attr-defined]
    fisher_p = float(fr.pvalue)  # type: ignore[attr-defined]

    # Table-preserving permutation null: within each table, keep which
    # codons change but permute their target AAs among that table's
    # observed targets. Recompute "creates novel disconnections."
    import random as _rng_mod

    from codon_topo.analysis.reassignment_db import build_reassignment_db as _build_db

    by_table: dict[int, list[tuple[str, str]]] = {}
    for e in _build_db():
        by_table.setdefault(e.table_id, []).append((e.codon, e.target_aa))

    rng = _rng_mod.Random(135325)
    n_perm = 10_000
    n_perm_extreme = 0
    for _ in range(n_perm):
        perm_creates = 0
        perm_seen: set[tuple[str, str]] = set()
        for _tid, events in by_table.items():
            codons = [c for c, _t in events]
            targets = [t for _c, t in events]
            rng.shuffle(targets)
            for c, t in zip(codons, targets):
                key = (c, t)
                if key in perm_seen:
                    continue
                perm_seen.add(key)
                variant = dict(STANDARD)
                variant[c] = t
                new_cat = disconnection_catalogue(variant)
                new_disc = {entry["aa"] for entry in new_cat}
                if new_disc - standard_disc_aas:
                    perm_creates += 1
        perm_rate = perm_creates / max(len(perm_seen), 1)
        if perm_rate <= rate_observed:
            n_perm_extreme += 1

    perm_p = (n_perm_extreme + 1) / (n_perm + 1)

    return {
        "observed_creates_disc": a,
        "observed_no_disc": b,
        "observed_total": n_obs,
        "possible_creates_disc": K,
        "possible_no_disc": N - K,
        "possible_total": N,
        "rate_observed": rate_observed,
        "rate_possible": rate_possible,
        "hypergeom_p": hypergeom_p,
        "permutation_p": perm_p,
        "n_permutations": n_perm,
        "odds_ratio": odds_ratio,
        "fisher_p": fisher_p,
        "hypothesis": (
            "Natural reassignments are depleted for topology-breaking "
            "changes relative to the space of possible reassignments"
        ),
        "caveat": (
            "Hypergeometric p assumes iid draws; permutation p preserves "
            "table structure. Both reported for transparency."
        ),
    }


# ====================================================================
# Phylogenetic lineage groupings per Sengupta, Yang & Higgs (2007)
# J Mol Evol 64:662-688, and Abascal et al. (2006) PLoS Biol 4:e127
# ====================================================================

# Maps each NCBI table to its major independent lineage.
# Tables sharing a lineage reflect the SAME ancestral reassignment event.
# ~15 independent lineages across 25 tables.
PHYLOGENETIC_LINEAGES: dict[int, str] = {
    1: "standard",  # no reassignment
    2: "vertebrate_mito",
    3: "yeast_mito",
    4: "mold_protozoan_mito",
    5: "invertebrate_mito",
    6: "ciliate_nuclear_UAR_Gln",  # multiple independent origins within
    9: "echinoderm_mito",
    10: "euplotid_nuclear_UGA_Cys",
    11: "bacterial_chloroplast",  # essentially standard
    12: "CUG_Ser_clade_1",
    13: "ascidian_mito",
    14: "flatworm_mito",
    15: "blepharisma_nuclear_UGA_Trp",
    16: "chlorophycean_mito",
    21: "trematode_mito",
    22: "scenedesmus_mito",
    23: "thraustochytrium_mito",
    24: "rhabdopleurid_mito",
    26: "CUG_Ala_pachysolen",
    27: "karyorelict_nuclear",
    28: "condylostoma_nuclear",
    29: "mesodinium_nuclear",
    30: "peritrich_nuclear",
    31: "blastocrithidia_nuclear",
    33: "cephalodiscus_mito",
}

# Major clade groups for exclusion sensitivity analysis
CLADE_GROUPS: dict[str, list[int]] = {
    "all_ciliates": [6, 10, 15, 27, 28, 29, 30],
    "all_yeast_mito": [3],
    "all_CUG_clade": [12, 26],
    "all_metazoan_mito": [2, 5, 9, 13, 14, 21],
    "all_algal_mito": [16, 22],
    "all_protist_mito": [4, 23],
    "all_hemichordate_mito": [24, 33],
}


def topology_avoidance_k43(seed: int = 135325) -> dict:
    """Topology avoidance test using K4^3 (nucleotide-level) adjacency.

    Hamming-1 in GF(2)^6 misses ~1/3 of single-nucleotide mutations. Here we
    define amino acid codon-family graphs using full single-nucleotide
    adjacency: two codons are neighbors iff they differ at exactly one
    nucleotide position (encoding-independent).

    A reassignment is topology-breaking if it increases the number of
    connected components in ANY amino acid's codon graph under K4^3.
    """
    import random as _rng_mod

    from scipy.stats import hypergeom

    from codon_topo.analysis.reassignment_db import build_reassignment_db
    from codon_topo.core.encoding import ALL_CODONS, nucleotide_distance

    def _k43_components(code: dict[str, str]) -> dict[str, int]:
        """Count connected components per AA under K4^3 adjacency."""
        from collections import defaultdict

        aa_codons: dict[str, list[str]] = defaultdict(list)
        for c, aa in code.items():
            if aa != "Stop":
                aa_codons[aa].append(c)

        result = {}
        for aa, codons in aa_codons.items():
            if len(codons) < 2:
                result[aa] = 1 if codons else 0
                continue
            # Union-find
            parent = list(range(len(codons)))

            def find(x: int) -> int:
                while parent[x] != x:
                    parent[x] = parent[parent[x]]
                    x = parent[x]
                return x

            def union(x: int, y: int) -> None:
                px, py = find(x), find(y)
                if px != py:
                    parent[px] = py

            for i in range(len(codons)):
                for j in range(i + 1, len(codons)):
                    if nucleotide_distance(codons[i], codons[j]) == 1:
                        union(i, j)
            result[aa] = len(set(find(i) for i in range(len(codons))))
        return result

    standard_comps = _k43_components(STANDARD)

    # All possible single-codon reassignments
    all_aas = sorted(set(STANDARD.values()))
    possible_breaks = 0
    possible_total = 0
    for codon in ALL_CODONS:
        orig = STANDARD[codon]
        for new_aa in all_aas:
            if new_aa == orig:
                continue
            variant = dict(STANDARD)
            variant[codon] = new_aa
            var_comps = _k43_components(variant)
            breaks = any(
                var_comps.get(aa, 0) > standard_comps.get(aa, 0)
                for aa in set(list(var_comps.keys()) + list(standard_comps.keys()))
            )
            if breaks:
                possible_breaks += 1
            possible_total += 1

    # Observed natural reassignments
    db = build_reassignment_db()
    observed_breaks = 0
    observed_total = 0
    seen: set[tuple[str, str]] = set()
    for e in db:
        key = (e.codon, e.target_aa)
        if key in seen:
            continue
        seen.add(key)
        variant = dict(STANDARD)
        variant[e.codon] = e.target_aa
        var_comps = _k43_components(variant)
        breaks = any(
            var_comps.get(aa, 0) > standard_comps.get(aa, 0)
            for aa in set(list(var_comps.keys()) + list(standard_comps.keys()))
        )
        if breaks:
            observed_breaks += 1
        observed_total += 1

    rate_obs = observed_breaks / max(observed_total, 1)
    rate_poss = possible_breaks / max(possible_total, 1)

    hyper_p = float(
        hypergeom.cdf(observed_breaks, possible_total, possible_breaks, observed_total)
    )

    # Permutation test (table-preserving)
    by_table: dict[int, list[tuple[str, str]]] = {}
    for e in db:
        by_table.setdefault(e.table_id, []).append((e.codon, e.target_aa))

    rng = _rng_mod.Random(seed)
    n_perm = 10_000
    n_extreme = 0
    for _ in range(n_perm):
        perm_breaks = 0
        perm_seen: set[tuple[str, str]] = set()
        for _tid, events in by_table.items():
            codons = [c for c, _t in events]
            targets = [t for _c, t in events]
            rng.shuffle(targets)
            for c, t in zip(codons, targets):
                key = (c, t)
                if key in perm_seen:
                    continue
                perm_seen.add(key)
                variant = dict(STANDARD)
                variant[c] = t
                var_comps = _k43_components(variant)
                if any(
                    var_comps.get(aa, 0) > standard_comps.get(aa, 0)
                    for aa in set(list(var_comps.keys()) + list(standard_comps.keys()))
                ):
                    perm_breaks += 1
        perm_rate = perm_breaks / max(len(perm_seen), 1)
        if perm_rate <= rate_obs:
            n_extreme += 1

    perm_p = (n_extreme + 1) / (n_perm + 1)

    # Risk ratio + CI (log-normal approximation)
    import math

    rr = rate_obs / max(rate_poss, 1e-10)
    # Standard error of log(RR) for two independent proportions
    if observed_breaks > 0 and observed_total > observed_breaks:
        se_log_rr = math.sqrt(
            (1 / max(observed_breaks, 1) - 1 / observed_total)
            + (1 / max(possible_breaks, 1) - 1 / possible_total)
        )
        rr_ci_lo = math.exp(math.log(rr) - 1.96 * se_log_rr)
        rr_ci_hi = math.exp(math.log(rr) + 1.96 * se_log_rr)
    else:
        se_log_rr = float("inf")
        rr_ci_lo = 0.0
        rr_ci_hi = float("inf")

    return {
        "adjacency": "K4^3 (nucleotide-level, encoding-independent)",
        "observed_breaks": observed_breaks,
        "observed_total": observed_total,
        "rate_observed": rate_obs,
        "possible_breaks": possible_breaks,
        "possible_total": possible_total,
        "rate_possible": rate_poss,
        "depletion_fold": rate_poss / max(rate_obs, 0.001),
        "risk_ratio": rr,
        "risk_ratio_ci_95": (rr_ci_lo, rr_ci_hi),
        "hypergeom_p": hyper_p,
        "permutation_p": perm_p,
        "n_permutations": n_perm,
    }


def topology_avoidance_phylogenetic_sensitivity() -> dict:
    """Re-run topology avoidance test with phylogenetic corrections.

    NCBI tables share evolutionary ancestry, so reassignment events across
    tables are not strictly independent. This routine quantifies the impact
    of that non-independence on the topology-avoidance result.

    Three analyses:
    1. Lineage-collapsed: one event per unique (lineage, codon, target_aa)
    2. Clade-exclusion: iteratively remove major clades and retest
    3. Conservative: one event per lineage (most extreme reduction)

    Reference: Sengupta, Yang & Higgs (2007) J Mol Evol 64:662-688
    """
    from scipy.stats import hypergeom

    from codon_topo.analysis.reassignment_db import build_reassignment_db

    db = build_reassignment_db()
    standard_cat = disconnection_catalogue(STANDARD)
    standard_disc_aas = {e["aa"] for e in standard_cat}

    # Compute full landscape once
    landscape = single_reassignment_landscape(STANDARD)
    possible_creates_disc = 0
    possible_total = 0
    for entry in landscape:
        variant = dict(STANDARD)
        variant[entry["codon"]] = entry["new_aa"]
        new_cat = disconnection_catalogue(variant)
        new_disc_aas = {e["aa"] for e in new_cat}
        if new_disc_aas - standard_disc_aas:
            possible_creates_disc += 1
        possible_total += 1
    rate_possible = possible_creates_disc / max(possible_total, 1)

    def _count_disc_events(
        events: list,
    ) -> tuple[int, int]:
        """Count disc-creating vs non-disc events."""
        creates = 0
        no_disc = 0
        seen: set[tuple[str, str]] = set()
        for e in events:
            key = (e.codon, e.target_aa)
            if key in seen:
                continue
            seen.add(key)
            variant = dict(STANDARD)
            variant[e.codon] = e.target_aa
            new_cat = disconnection_catalogue(variant)
            new_disc = {entry["aa"] for entry in new_cat}
            if new_disc - standard_disc_aas:
                creates += 1
            else:
                no_disc += 1
        return creates, no_disc

    # 1. Lineage-collapsed: de-dup by (lineage, codon, target_aa)
    seen_lineage: set[tuple[str, str, str]] = set()
    lineage_events = []
    for e in db:
        lineage = PHYLOGENETIC_LINEAGES.get(e.table_id, f"unknown_{e.table_id}")
        key = (lineage, e.codon, e.target_aa)
        if key not in seen_lineage:
            seen_lineage.add(key)
            lineage_events.append(e)

    lc_creates, lc_no_disc = _count_disc_events(lineage_events)
    lc_total = lc_creates + lc_no_disc
    lc_rate = lc_creates / max(lc_total, 1)
    lc_hyper_p = float(
        hypergeom.cdf(lc_creates, possible_total, possible_creates_disc, lc_total)
    )

    # 2. Clade-exclusion sensitivity
    clade_results = []
    for clade_name, table_ids in CLADE_GROUPS.items():
        excluded_tables = set(table_ids)
        filtered = [e for e in db if e.table_id not in excluded_tables]
        fc, fn = _count_disc_events(filtered)
        ft = fc + fn
        if ft == 0:
            continue
        fr = fc / ft
        fp = float(hypergeom.cdf(fc, possible_total, possible_creates_disc, ft))
        clade_results.append(
            {
                "excluded_clade": clade_name,
                "excluded_tables": sorted(table_ids),
                "n_events_remaining": ft,
                "creates_disc": fc,
                "rate_observed": fr,
                "hypergeom_p": fp,
                "significant_p05": fp < 0.05,
            }
        )

    return {
        "lineage_collapsed": {
            "n_events": lc_total,
            "creates_disc": lc_creates,
            "rate_observed": lc_rate,
            "rate_possible": rate_possible,
            "depletion_fold": rate_possible / max(lc_rate, 0.001),
            "hypergeom_p": lc_hyper_p,
        },
        "clade_exclusion": clade_results,
        "all_clade_exclusions_significant": all(
            r["significant_p05"] for r in clade_results
        ),
        "method": (
            "Phylogenetic sensitivity per Sengupta et al. 2007 (J Mol Evol "
            "64:662-688). Lineage-collapsed: one event per (lineage, codon, "
            "target_aa). Clade-exclusion: remove major taxonomic groups and "
            "retest."
        ),
    }


# ====================================================================
# Definitions audit: 2 adjacencies x 2 topology-breaking definitions
# ====================================================================


def _q6_components_per_aa(
    code: dict[str, str],
    encoding: dict[str, tuple[int, int]] | None = None,
) -> dict[str, int]:
    """Count Q_6 (Hamming-1, encoding-dependent) connected components per AA."""
    from collections import defaultdict

    from codon_topo.core.encoding import codon_to_vector, DEFAULT_ENCODING
    from codon_topo.core.homology import _partition

    enc = encoding or DEFAULT_ENCODING
    aa_codons: dict[str, list[str]] = defaultdict(list)
    for c, aa in code.items():
        if aa != "Stop":
            aa_codons[aa].append(c)

    out: dict[str, int] = {}
    for aa, codons in aa_codons.items():
        if not codons:
            out[aa] = 0
        elif len(codons) == 1:
            out[aa] = 1
        else:
            vectors = [codon_to_vector(c, enc) for c in codons]
            out[aa] = len(_partition(vectors, 1))
    return out


def _k43_components_per_aa(code: dict[str, str]) -> dict[str, int]:
    """Count K_4^3 (nucleotide-level, encoding-independent) components per AA."""
    from collections import defaultdict

    from codon_topo.core.encoding import nucleotide_distance

    aa_codons: dict[str, list[str]] = defaultdict(list)
    for c, aa in code.items():
        if aa != "Stop":
            aa_codons[aa].append(c)

    out: dict[str, int] = {}
    for aa, codons in aa_codons.items():
        if not codons:
            out[aa] = 0
            continue
        if len(codons) == 1:
            out[aa] = 1
            continue
        parent = list(range(len(codons)))

        def find(x: int) -> int:
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for i in range(len(codons)):
            for j in range(i + 1, len(codons)):
                if nucleotide_distance(codons[i], codons[j]) == 1:
                    pi, pj = find(i), find(j)
                    if pi != pj:
                        parent[pi] = pj
        out[aa] = len({find(i) for i in range(len(codons))})
    return out


def _classify_move(
    standard_q6: dict[str, int],
    standard_k43: dict[str, int],
    standard_q6_disc: set[str],
    standard_k43_disc: set[str],
    variant_q6: dict[str, int],
    variant_k43: dict[str, int],
) -> dict[str, bool]:
    """Classify a single candidate move under all four definitions.

    Returns four bool keys:
      - q6_new_disc:  variant has AAs newly disconnected (variant_q6[aa] > 1
                      AND aa not in standard_q6_disc)
      - q6_beta0:     variant has any AA whose Q_6 component count exceeds
                      its standard value (Δβ₀ > 0)
      - k43_new_disc: variant has AAs newly disconnected under K_4^3
      - k43_beta0:    variant has any AA whose K_4^3 component count exceeds
                      its standard value
    """
    variant_q6_disc = {aa for aa, n in variant_q6.items() if n > 1}
    variant_k43_disc = {aa for aa, n in variant_k43.items() if n > 1}
    novel_q6 = variant_q6_disc - standard_q6_disc
    novel_k43 = variant_k43_disc - standard_k43_disc

    q6_beta0 = any(
        variant_q6.get(aa, 0) > standard_q6.get(aa, 0)
        for aa in set(variant_q6) | set(standard_q6)
    )
    k43_beta0 = any(
        variant_k43.get(aa, 0) > standard_k43.get(aa, 0)
        for aa in set(variant_k43) | set(standard_k43)
    )

    return {
        "q6_new_disc": bool(novel_q6),
        "q6_beta0": q6_beta0,
        "k43_new_disc": bool(novel_k43),
        "k43_beta0": k43_beta0,
    }


def topology_definitions_audit() -> dict:
    """Compute the 2 (adjacency) x 2 (definition) sensitivity audit.

    Reviewer R1 noted that the manuscript Methods text claims the Δβ₀>0
    (any-component-increase) definition while Q_6 numbers (931 of 1280
    possible) actually correspond to the narrower "creates a new
    disconnection from a previously connected family" definition. This
    routine computes all four cells so the manuscript can pick one
    definition consistently and report the other as a sensitivity column.

    Cells (possible_breaks, observed_breaks, rate_obs, rate_poss,
    depletion, hypergeom_p, RR, CI):
      - Q_6 with new-disconnection (legacy primary)
      - Q_6 with Δβ₀>0
      - K_4^3 with new-disconnection
      - K_4^3 with Δβ₀>0 (legacy primary for K_4^3)
    """
    import math

    from scipy.stats import hypergeom

    from codon_topo.analysis.reassignment_db import build_reassignment_db
    from codon_topo.core.encoding import ALL_CODONS

    standard_q6 = _q6_components_per_aa(STANDARD)
    standard_k43 = _k43_components_per_aa(STANDARD)
    standard_q6_disc = {aa for aa, n in standard_q6.items() if n > 1}
    standard_k43_disc = {aa for aa, n in standard_k43.items() if n > 1}

    cells: dict[str, dict[str, int]] = {
        "q6_new_disc": {"possible": 0, "observed": 0},
        "q6_beta0": {"possible": 0, "observed": 0},
        "k43_new_disc": {"possible": 0, "observed": 0},
        "k43_beta0": {"possible": 0, "observed": 0},
    }
    possible_total = 0

    all_aas = sorted(set(STANDARD.values()))
    for codon in ALL_CODONS:
        orig = STANDARD[codon]
        for new_aa in all_aas:
            if new_aa == orig:
                continue
            possible_total += 1
            variant = dict(STANDARD)
            variant[codon] = new_aa
            v_q6 = _q6_components_per_aa(variant)
            v_k43 = _k43_components_per_aa(variant)
            flags = _classify_move(
                standard_q6,
                standard_k43,
                standard_q6_disc,
                standard_k43_disc,
                v_q6,
                v_k43,
            )
            for k, v in flags.items():
                if v:
                    cells[k]["possible"] += 1

    db = build_reassignment_db()
    seen: set[tuple[str, str]] = set()
    observed_total = 0
    for e in db:
        key = (e.codon, e.target_aa)
        if key in seen:
            continue
        seen.add(key)
        observed_total += 1
        variant = dict(STANDARD)
        variant[e.codon] = e.target_aa
        v_q6 = _q6_components_per_aa(variant)
        v_k43 = _k43_components_per_aa(variant)
        flags = _classify_move(
            standard_q6,
            standard_k43,
            standard_q6_disc,
            standard_k43_disc,
            v_q6,
            v_k43,
        )
        for k, v in flags.items():
            if v:
                cells[k]["observed"] += 1

    audit_rows = []
    for cell_name, c in cells.items():
        K = c["possible"]
        x = c["observed"]
        N = possible_total
        n = observed_total
        rate_obs = x / max(n, 1)
        rate_poss = K / max(N, 1)
        depletion = rate_poss / max(rate_obs, 1e-10)
        hp = float(hypergeom.cdf(x, N, K, n))
        if x > 0 and n > x and K > 0 and N > K:
            se_log_rr = math.sqrt((1 / max(x, 1) - 1 / n) + (1 / max(K, 1) - 1 / N))
            rr = rate_obs / max(rate_poss, 1e-10)
            rr_lo = math.exp(math.log(rr) - 1.96 * se_log_rr)
            rr_hi = math.exp(math.log(rr) + 1.96 * se_log_rr)
        else:
            rr = float("nan")
            rr_lo, rr_hi = float("nan"), float("nan")
        audit_rows.append(
            {
                "cell": cell_name,
                "adjacency": "Q_6" if cell_name.startswith("q6") else "K_4^3",
                "definition": (
                    "new_disconnection_in_previously_connected_family"
                    if cell_name.endswith("new_disc")
                    else "increase_in_components_(Δβ₀>0)"
                ),
                "possible_breaks": K,
                "possible_total": N,
                "rate_possible": rate_poss,
                "observed_breaks": x,
                "observed_total": n,
                "rate_observed": rate_obs,
                "depletion_fold": depletion,
                "risk_ratio": rr,
                "risk_ratio_ci_95": (rr_lo, rr_hi),
                "hypergeom_p": hp,
            }
        )

    return {
        "method": (
            "Two-by-two audit of topology-breaking under (Q_6, K_4^3) "
            "x (new-disconnection-in-previously-connected-family, Δβ₀>0). "
            "All four cells share the same denominators (1280 candidate "
            "moves; 28 de-duplicated observed events) but differ in how "
            "'topology-breaking' is defined. Reviewer R1 noted the "
            "manuscript Methods text described Δβ₀>0 while Q_6 reported "
            "counts matched the new-disconnection definition; both are "
            "now reported transparently."
        ),
        "audit_rows": audit_rows,
    }


def topology_avoidance_q6_encoding_sweep() -> dict:
    """Recompute Q_6 topology avoidance under all 24 base-to-bit encodings.

    K_4^3 adjacency is encoding-independent and does not need this sweep.
    The Q_6 result depends on which 192 of the 288 single-nucleotide edges
    are designated as Hamming-1 (vs the 96 Hamming-2 'diagonals'); reviewer
    R1.B asked for confirmation that the depletion holds across all 24
    bijections, not just the default C=00, U=01, A=10, G=11.

    Uses the new-disconnection definition (legacy primary for Q_6).
    """
    from scipy.stats import hypergeom

    from codon_topo.analysis.reassignment_db import build_reassignment_db
    from codon_topo.core.encoding import ALL_CODONS, all_encodings

    db = build_reassignment_db()

    rows = []
    for idx, enc in enumerate(all_encodings()):
        # encoding label as concatenated base->bits string
        label = "".join(f"{b}={enc[b][0]}{enc[b][1]}" for b in ("C", "U", "A", "G"))
        standard_q6 = _q6_components_per_aa(STANDARD, enc)
        standard_disc = {aa for aa, n in standard_q6.items() if n > 1}

        possible_breaks = 0
        possible_total = 0
        all_aas = sorted(set(STANDARD.values()))
        for codon in ALL_CODONS:
            orig = STANDARD[codon]
            for new_aa in all_aas:
                if new_aa == orig:
                    continue
                possible_total += 1
                variant = dict(STANDARD)
                variant[codon] = new_aa
                v_q6 = _q6_components_per_aa(variant, enc)
                variant_disc = {aa for aa, n in v_q6.items() if n > 1}
                if variant_disc - standard_disc:
                    possible_breaks += 1

        observed_breaks = 0
        observed_total = 0
        seen: set[tuple[str, str]] = set()
        for e in db:
            key = (e.codon, e.target_aa)
            if key in seen:
                continue
            seen.add(key)
            observed_total += 1
            variant = dict(STANDARD)
            variant[e.codon] = e.target_aa
            v_q6 = _q6_components_per_aa(variant, enc)
            variant_disc = {aa for aa, n in v_q6.items() if n > 1}
            if variant_disc - standard_disc:
                observed_breaks += 1

        rate_obs = observed_breaks / max(observed_total, 1)
        rate_poss = possible_breaks / max(possible_total, 1)
        depletion = rate_poss / max(rate_obs, 1e-10)
        hp = float(
            hypergeom.cdf(
                observed_breaks, possible_total, possible_breaks, observed_total
            )
        )
        rows.append(
            {
                "encoding_index": idx,
                "encoding_label": label,
                "rate_observed": rate_obs,
                "rate_possible": rate_poss,
                "depletion_fold": depletion,
                "hypergeom_p": hp,
                "observed_breaks": observed_breaks,
                "possible_breaks": possible_breaks,
            }
        )

    rates_obs = [r["rate_observed"] for r in rows]
    rates_poss = [r["rate_possible"] for r in rows]
    depletions = [r["depletion_fold"] for r in rows]
    pvals = [r["hypergeom_p"] for r in rows]
    return {
        "method": (
            "Q_6 topology-avoidance recomputed under all 24 base-to-bit "
            "encodings (the K_4^3 result is encoding-independent and "
            "needs no sweep). Uses the new-disconnection definition. "
            "All 24 encodings share the same denominators (1280 candidate "
            "moves; 28 de-duplicated observed events). The default "
            "encoding (C=00, U=01, A=10, G=11) is encoding_index 0."
        ),
        "n_encodings": len(rows),
        "rate_obs_min": min(rates_obs),
        "rate_obs_median": sorted(rates_obs)[len(rates_obs) // 2],
        "rate_obs_max": max(rates_obs),
        "rate_poss_min": min(rates_poss),
        "rate_poss_median": sorted(rates_poss)[len(rates_poss) // 2],
        "rate_poss_max": max(rates_poss),
        "depletion_min": min(depletions),
        "depletion_median": sorted(depletions)[len(depletions) // 2],
        "depletion_max": max(depletions),
        "p_max": max(pvals),
        "p_median": sorted(pvals)[len(pvals) // 2],
        "p_min": min(pvals),
        "all_significant_at_0p01": all(p < 0.01 for p in pvals),
        "all_significant_at_0p05": all(p < 0.05 for p in pvals),
        "rows": rows,
    }


def topology_denominator_sensitivity() -> dict:
    """Report the four candidate-universe definitions side-by-side.

    Reviewer R1.D / R1.2 asked for an explicit denominator-sensitivity
    table covering:
      1. 21-label targets, identity-excluded: |M| = 64 x 20 = 1280
         (each codon gets 20 alternative labels = 19 AAs + Stop, OR
         20 AAs if currently Stop)
      2. Amino-acid-only targets, identity-excluded: |M| = 1219
         (61 sense codons x 19 AA alternatives + 3 stop codons x 20)
      3. AA targets including no-op (current label allowed): |M| = 64 x 20
         = 1280, but where 64 of those are no-ops (codon -> current AA)
      4. Stop-inclusive with no-ops: |M| = 64 x 21 = 1344
    """
    from codon_topo.core.encoding import ALL_CODONS

    n_codons = 64
    n_sense = sum(1 for c in ALL_CODONS if STANDARD[c] != "Stop")
    n_stop = n_codons - n_sense

    # Define each universe
    rows = [
        {
            "universe": "U1_21label_no_identity",
            "description": (
                "y in {A_20 union Stop}, y != C(x); each codon has "
                "exactly 20 alternative labels"
            ),
            "size": n_codons * 20,  # 1280
            "is_primary": True,
        },
        {
            "universe": "U2_AA_only_no_identity",
            "description": (
                "y in A_20, y != C(x); 61 sense codons get 19 alts, "
                "3 stop codons get 20 alts"
            ),
            "size": n_sense * 19 + n_stop * 20,  # 1219 for standard code
        },
        {
            "universe": "U3_AA_with_noop",
            "description": (
                "y in A_20, no identity restriction; 64 candidates are "
                "no-ops (codon -> current AA)"
            ),
            "size": n_codons * 20,  # 1280; semantically different from U1
        },
        {
            "universe": "U4_stop_inclusive_with_noop",
            "description": (
                "y in {A_20 union Stop}, no identity restriction; 64 no-ops included"
            ),
            "size": n_codons * 21,  # 1344
        },
    ]
    return {
        "method": (
            "Reviewer R1.D / R1.2 denominator-sensitivity audit: four "
            "candidate-universe definitions for the topology-avoidance "
            "denominator. Primary for the manuscript is U1 "
            "(21-label, identity-excluded) = 1280 single-codon relabelings. "
            "U2 (AA-only) gives 1219 for the standard code; the difference "
            "matters when interpretating the hypergeometric parameter K. "
            "U3 and U4 are listed for completeness but are not used in "
            "the primary analysis."
        ),
        "rows": rows,
    }
