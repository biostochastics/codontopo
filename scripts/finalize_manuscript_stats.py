#!/usr/bin/env python3.11
"""Splice new sensitivity blocks into manuscript_stats.json.

Used as a one-shot post-processor when `codon-topo all` was launched against
an older cli.py that did not yet include the new sensitivity keys
(topology_audit, topology_q6_encoding_sweep, condlogit.encoding_robustness,
topology_avoidance.k43_definitions_audit). Reads the per-analysis JSONs
and injects the new blocks into manuscript_stats.json without re-running
any heavy computation.

Idempotent: running multiple times produces the same result.
"""

from __future__ import annotations

import json
from pathlib import Path

OUTPUT = Path("/Users/biostochastics/clayworth/output")


def load_json(path: Path) -> dict:
    with path.open() as f:
        return json.load(f)


def write_json(path: Path, data: dict) -> None:
    with path.open("w") as f:
        json.dump(data, f, indent=2, default=str)


def main() -> None:
    stats = load_json(OUTPUT / "manuscript_stats.json")
    topo = load_json(OUTPUT / "topology_avoidance.json")
    evosim = load_json(OUTPUT / "evolutionary_simulation.json")

    # 1. Inject topology_audit (2x2 sensitivity)
    if "definitions_audit" in topo:
        stats["topology_audit"] = topo["definitions_audit"]
        print("  [+] topology_audit added")

    # 2. Inject topology_q6_encoding_sweep
    if "Q6_encoding_sweep" in topo:
        stats["topology_q6_encoding_sweep"] = topo["Q6_encoding_sweep"]
        print("  [+] topology_q6_encoding_sweep added")

    # 3. Inject denominator_sensitivity
    if "denominator_sensitivity" in topo:
        stats["topology_denominator_sensitivity"] = topo["denominator_sensitivity"]
        print("  [+] topology_denominator_sensitivity added")

    # 4. Inject encoding_robustness from evolutionary_simulation.json into
    #    the condlogit block in manuscript_stats.json
    if "encoding_robustness" in evosim:
        if "condlogit" not in stats:
            stats["condlogit"] = {}
        stats["condlogit"]["encoding_robustness"] = evosim["encoding_robustness"]
        print("  [+] condlogit.encoding_robustness added")

    # 5. Inject the K_4^3 model fits if missing (for completeness)
    if "model_fits" in evosim:
        if "condlogit" not in stats:
            stats["condlogit"] = {}
        existing_fits = stats["condlogit"].get("model_fits", {})
        for name, f in evosim["model_fits"].items():
            if name not in existing_fits:
                existing_fits[name] = f
                print(f"  [+] condlogit.model_fits.{name} added")
        stats["condlogit"]["model_fits"] = existing_fits

    # 6. Inject lr_tests (extended) if available
    if "likelihood_ratio_tests" in evosim:
        if "condlogit" not in stats:
            stats["condlogit"] = {}
        existing_lr = stats["condlogit"].get("lr_tests", {})
        for name, lr in evosim["likelihood_ratio_tests"].items():
            if name not in existing_lr:
                existing_lr[name] = lr
                print(f"  [+] condlogit.lr_tests.{name} added")
        stats["condlogit"]["lr_tests"] = existing_lr

    write_json(OUTPUT / "manuscript_stats.json", stats)
    print(f"\nUpdated {OUTPUT / 'manuscript_stats.json'}")


if __name__ == "__main__":
    main()
