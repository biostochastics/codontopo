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
    odds_ratio, fisher_p = fisher_exact(
        [[a, b], [max(c_complement, 0), max(d_complement, 0)]],
        alternative="less",
    )

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
        "odds_ratio": float(odds_ratio),
        "fisher_p": float(fisher_p),
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

    Addresses the reviewer concern that Hamming-1 in GF(2)^6 misses ~1/3
    of single-nucleotide mutations. Here we define amino acid codon-family
    graphs using full single-nucleotide adjacency: two codons are neighbors
    iff they differ at exactly one nucleotide position (encoding-independent).

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
            idx = {c: i for i, c in enumerate(codons)}
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
                    for aa in set(
                        list(var_comps.keys()) + list(standard_comps.keys())
                    )
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

    Addresses the reviewer concern about non-independence of reassignment
    events across NCBI tables sharing evolutionary ancestry.

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
        fp = float(
            hypergeom.cdf(fc, possible_total, possible_creates_disc, ft)
        )
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
