"""cBioPortal API client for KRAS mutation co-occurrence analysis.

Queries cBioPortal (public, no API key required) for KRAS G12 mutations
and tests whether Fano-line-predicted co-mutations are enriched.

WS4 is the fast-fail gate: a null result kills the clinical prediction
track only, not the underlying mathematical framework.
"""

import logging
from collections import defaultdict
from dataclasses import dataclass

from scipy.stats import fisher_exact

from codon_topo.core.genetic_codes import STANDARD

logger = logging.getLogger(__name__)

# Single-letter amino acid codes for parsing protein change strings
AA_THREE_TO_ONE = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
}

AA_ONE_TO_THREE = {v: k for k, v in AA_THREE_TO_ONE.items()}

# KRAS wild-type codon 12 is GGT (DNA) = GGU (RNA) = Gly
KRAS_G12_VARIANTS: dict[str, str] = {
    "G12V": "GUU",  # Val
    "G12D": "GAU",  # Asp
    "G12A": "GCU",  # Ala
    "G12R": "CGU",  # Arg
    "G12C": "UGU",  # Cys
    "G12S": "AGU",  # Ser
}


@dataclass
class KRASMutation:
    """A parsed mutation record."""

    sample_id: str
    protein_change: str
    mutation_type: str
    ref_aa: str = ""
    alt_aa: str = ""
    position: int = 0


@dataclass
class CoOccurrence:
    """Co-occurring mutation with a KRAS G12 variant."""

    sample_id: str
    protein_change: str
    ref_aa: str


def _parse_protein_change(pc: str) -> tuple[str, int, str]:
    """Parse 'H47R' into (ref_aa='H', position=47, alt_aa='R').

    Returns ('', 0, '') if unparseable.
    """
    if not pc or len(pc) < 3:
        return ("", 0, "")
    ref = pc[0]
    alt = pc[-1]
    try:
        pos = int(pc[1:-1])
    except ValueError:
        return ("", 0, "")
    if ref.isalpha() and alt.isalpha():
        return (ref, pos, alt)
    return ("", 0, "")


def fano_predictions_for_kras() -> dict[str, dict]:
    """Compute Fano-line predictions for each KRAS G12 variant."""
    from codon_topo.core.fano import fano_partner

    predictions = {}
    for variant, mutant_codon in KRAS_G12_VARIANTS.items():
        partner = fano_partner("GGU", mutant_codon)
        partner_aa = STANDARD.get(partner, "Unknown")
        partner_aa_1 = AA_THREE_TO_ONE.get(partner_aa, "?")
        predictions[variant] = {
            "wt_codon": "GGU",
            "mutant_codon": mutant_codon,
            "fano_partner_codon": partner,
            "fano_partner_aa": partner_aa,
            "fano_partner_aa_1letter": partner_aa_1,
        }
    return predictions


def parse_mutation_data(raw_mutations: list[dict]) -> list[KRASMutation]:
    """Parse raw cBioPortal mutation dicts into KRASMutation objects."""
    results = []
    for m in raw_mutations:
        pc = m.get("proteinChange", "")
        ref, pos, alt = _parse_protein_change(pc)
        results.append(
            KRASMutation(
                sample_id=m.get("sampleId", ""),
                protein_change=pc,
                mutation_type=m.get("mutationType", ""),
                ref_aa=ref,
                alt_aa=alt,
                position=pos,
            )
        )
    return results


def compute_cooccurrence(raw_mutations: list[dict]) -> dict:
    """Compute co-occurrence of mutations with KRAS G12 variants."""
    by_sample: dict[str, list[dict]] = defaultdict(list)
    for m in raw_mutations:
        by_sample[m.get("sampleId", "")].append(m)

    results: dict[str, dict] = {}
    for variant in KRAS_G12_VARIANTS:
        variant_samples = []
        co_muts: list[dict] = []
        for sid, muts in by_sample.items():
            has_variant = any(m.get("proteinChange", "") == variant for m in muts)
            if has_variant:
                variant_samples.append(sid)
                for m in muts:
                    pc = m.get("proteinChange", "")
                    if pc != variant:
                        ref, pos, alt = _parse_protein_change(pc)
                        co_muts.append(
                            {
                                "sample_id": sid,
                                "protein_change": pc,
                                "ref_aa": ref,
                                "alt_aa": alt,
                                "position": pos,
                            }
                        )
        results[variant] = {
            "n_samples": len(variant_samples),
            "co_mutations": co_muts,
        }
    return results


class CBioPortalClient:
    """Client for the cBioPortal REST API."""

    def __init__(
        self,
        base_url: str = "https://www.cbioportal.org/api",
        offline: bool = False,
    ):
        self.base_url = base_url.rstrip("/")
        self.offline = offline

    def get_kras_mutations(
        self,
        study_id: str = "msk_impact_2017",
    ) -> list[dict]:
        """Fetch KRAS mutations from a cBioPortal study."""
        if self.offline:
            return []

        import requests  # type: ignore[import-untyped]

        url = f"{self.base_url}/molecular-profiles/{study_id}_mutations/mutations"
        params = {
            "entrezGeneId": 3845,
            "projection": "DETAILED",
        }
        try:
            resp = requests.get(url, params=params, timeout=30)
            resp.raise_for_status()
            return resp.json()
        except requests.exceptions.RequestException as exc:
            logger.warning("cBioPortal query failed for study %s: %s", study_id, exc)
            return []


def fetch_kras_mutations(
    study_id: str = "msk_impact_2017",
    mock_data: list[dict] | None = None,
    offline: bool = False,
) -> list[KRASMutation]:
    """Fetch and parse KRAS mutations."""
    if mock_data is not None:
        return parse_mutation_data(mock_data)
    client = CBioPortalClient(offline=offline)
    raw = client.get_kras_mutations(study_id=study_id)
    return parse_mutation_data(raw)


def fano_enrichment_test(
    raw_mutations: list[dict],
    variant: str = "G12V",
    background_variants: list[str] | None = None,
) -> dict:
    """Test whether the Fano-predicted AA is enriched in co-mutations.

    Compares the frequency of the Fano-predicted AA as a reference residue
    at co-mutation sites in samples with the given KRAS G12 variant vs.
    samples with other G12 variants (background).

    Uses Fisher's exact test on a 2x2 contingency table:
                     Fano AA   Other AAs
    Target variant     a           b
    Background         c           d
    """
    preds = fano_predictions_for_kras()
    if variant not in preds:
        return {"variant": variant, "error": f"Unknown variant {variant}"}

    predicted_aa_1 = preds[variant]["fano_partner_aa_1letter"]

    if background_variants is None:
        background_variants = [v for v in KRAS_G12_VARIANTS if v != variant]

    cooc = compute_cooccurrence(raw_mutations)

    # Count Fano AA in target variant co-mutations
    target_co = cooc.get(variant, {}).get("co_mutations", [])
    a = sum(1 for m in target_co if m["ref_aa"] == predicted_aa_1)
    b = len(target_co) - a

    # Count Fano AA in background co-mutations
    bg_co: list[dict] = []
    for bv in background_variants:
        bg_co.extend(cooc.get(bv, {}).get("co_mutations", []))
    c = sum(1 for m in bg_co if m["ref_aa"] == predicted_aa_1)
    d = len(bg_co) - c

    # Fisher's exact test
    if a + b == 0 or c + d == 0:
        odds_ratio, p_value = float("nan"), 1.0
    else:
        odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")

    return {
        "variant": variant,
        "fano_predicted_aa": predicted_aa_1,
        "fano_predicted_aa_3letter": preds[variant]["fano_partner_aa"],
        "fano_partner_codon": preds[variant]["fano_partner_codon"],
        "observed_count": a,
        "total_co_mutations": a + b,
        "background_count": c,
        "background_total": c + d,
        "odds_ratio": float(odds_ratio),
        "fishers_p": float(p_value),
    }


def ws4_gate_decision(
    raw_mutations: list[dict],
    p_threshold: float = 0.01,
) -> dict:
    """WS4 fast-fail gate: test all KRAS G12 variants for Fano enrichment.

    The gate PASSES if at least one variant shows significant enrichment
    after Bonferroni correction for multiple testing across all G12 variants.

    A null result (gate fails) kills the clinical prediction track only,
    not the underlying algebraic framework.
    """
    n_tests = len(KRAS_G12_VARIANTS)
    corrected_threshold = p_threshold / n_tests

    results = []
    significant = []
    for variant in KRAS_G12_VARIANTS:
        r = fano_enrichment_test(raw_mutations, variant=variant)
        results.append(r)
        if r.get("fishers_p", 1.0) < corrected_threshold:
            significant.append(variant)

    return {
        "pass": len(significant) > 0,
        "p_threshold": p_threshold,
        "corrected_threshold": corrected_threshold,
        "n_tests": n_tests,
        "variants_tested": list(KRAS_G12_VARIANTS.keys()),
        "significant_variants": significant,
        "details": results,
    }
