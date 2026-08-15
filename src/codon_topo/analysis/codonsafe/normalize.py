"""Codon normalization and validation utilities for CodonSafe.

Handles DNA->RNA conversion, codon validation, and study-specific
normalization requirements.
"""

from __future__ import annotations

_DNA_TO_RNA = str.maketrans("Tt", "Uu")


def dna_to_rna(codon: str) -> str:
    """Convert a DNA-alphabet codon to RNA alphabet (T->U)."""
    return codon.translate(_DNA_TO_RNA)


def normalize_codon(codon: str) -> str:
    """Normalize a codon to uppercase RNA alphabet.

    Raises ValueError if the result contains invalid characters.
    """
    result = codon.upper().translate(_DNA_TO_RNA)
    valid = set("ACGU")
    if not set(result).issubset(valid):
        raise ValueError(
            f"Codon {codon!r} -> {result!r} contains invalid characters. "
            f"Expected RNA alphabet {{A, C, G, U}}."
        )
    if len(result) != 3:
        raise ValueError(f"Codon must be 3 characters, got {len(result)}: {codon!r}")
    return result


def validate_codon_pair(source: str, target: str) -> tuple[str, str]:
    """Normalize and validate a source/target codon pair."""
    return normalize_codon(source), normalize_codon(target)
