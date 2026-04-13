# Codebase Cleanup & CLI Consolidation — Next Session Prompt

**Purpose**: Drop this prompt into a fresh Claude Code session to clean up, wire, document, and CLI-ify the entire `codon_topo` codebase after the recent adversarial-review-driven refactor.

---

## CONTEXT (paste this at the top of the new session)

```
I just completed a round of adversarial review + refinements on a
scientific Python package `codon_topo` (genetic code as GF(2)^6 analysis).
Multi-model review (10 LLMs, 12 passes) led to several fixes:

  - Added null_model_c_extended() for per-encoding min-distance emission
  - Added bit_position_bias_weighted() with Ts/Tv (actually position-weighted) priors
  - Added trna_evidence.py module (4 disconnection organisms + 2 controls)
  - Added coloring_optimality.py with Freeland-Hurst block-preserving null
  - Added data/grantham.json (verified against Grantham 1974 Table 1)
  - Added reports/claim_hierarchy.py with structured claim status
  - Renamed misleading function names
  - Made p-value reporting honest (bounded)
  - Excluded non-standard target AAs from hamming path stats

324/324 tests pass. But the codebase needs CLEANUP. Specifically:

1. There is NO CLI entry point yet — add one (with subcommands for each
   primary analysis).
2. Several modules grew large during refinement — check for dead code,
   stale comments, old parameter names, misleading docstrings, and
   inconsistencies left over from the adversarial rewrite.
3. Some `# TODO` and `# note:` comments may be stale.
4. Deprecated function aliases (e.g., `_generate_random_code_block_preserving`)
   need removal after a deprecation cycle.
5. Documentation (README.md, CLAUDE.md, ARCHITECTURE.md) needs refresh.

## KEY FILES TO REVIEW

src/codon_topo/
  __init__.py                    # Public API re-exports — may be stale
  core/
    encoding.py                  # GF(2)^6 basics
    filtration.py                # Two-fold/four-fold (now known tautological)
    homology.py                  # Connected components / disconnection catalogue
    fano.py                      # XOR triples (terminology has changed)
    genetic_codes.py             # 25 NCBI tables
    embedding.py                 # Coordinate-wise root-of-unity map
  analysis/
    cosmic_query.py              # KRAS (failed clinical test, keep as negative record)
    depth_calibration.py         # n=6 non-result (demote or delete)
    null_models.py               # A, B, C, C_extended — check consistency
    reassignment_db.py           # bit_position_bias + weighted version
    synbio_feasibility.py        # Partly encoding-artifact — add caveat
    trna_evidence.py             # NEW — curated data, needs GtRNAdb refresh note
    coloring_optimality.py       # NEW central result module
  reports/
    claim_hierarchy.py           # NEW — structured claim statuses
    catalogue.py                 # Prediction catalogue (likely stale)
  data/
    grantham.json                # NEW — Grantham 1974 matrix
  visualization/
    data_export.py               # CSV exports — check columns still current
tests/
  test_*.py                      # 324 passing — should stay passing
README.md
CLAUDE.md
pyproject.toml                   # No CLI entry currently

## YOUR TASKS

Do these in order, running `python3.11 -m pytest tests/ -q` after each
task to confirm 324/324 still passes.

### 1. Add a proper CLI entry point

Create `src/codon_topo/cli.py` with a click/typer-based CLI. Register it
via `[project.scripts]` in pyproject.toml so `codon-topo` becomes a
shell command.

Required subcommands:

  codon-topo filtration [--table=1] [--all-tables]
    Run two-fold / four-fold filtration checks. Output CSV or pretty table.

  codon-topo disconnections [--table=1] [--extended]
    Emit the disconnection catalogue. --extended runs null_model_c_extended
    across 24 encodings.

  codon-topo coloring [--null=freeland_hurst|class_size] [--n=10000]
                       [--seed=135325] [--no-stops]
    Run hypercube coloring Monte Carlo. Output JSON with score, quantile,
    p-values.

  codon-topo bit-bias [--compartment=uniform|nuclear|mitochondrial]
    Run bit-position bias test. Report observed counts and p-value under
    chosen null.

  codon-topo trna
    Run tRNA duplication correlation test. Report 4/4 pattern.

  codon-topo kras [--offline]
    Run KRAS-Fano test (will return the p=1.0 negative).

  codon-topo all [--output-dir=./output]
    Run everything, write CSV and JSON reports to output/.

  codon-topo claims
    Print the claim hierarchy table (SUPPORTED, SUGGESTIVE, etc.).

Use Mantine-compatible-friendly JSON output format (snake_case, fully
typed). Print Markdown tables for interactive use. Use rich for pretty
output if available; fall back to stdlib.

### 2. Purge dead code and stale comments

Systematically:
  - Remove deprecated function aliases after adding a deprecation warning
    step: `_generate_random_code_block_preserving` (in coloring_optimality.py)
    should be deleted or raise DeprecationWarning.
  - Find and remove commented-out code blocks left from refactoring.
  - Check every docstring references the CURRENT function signature.
  - Find TODO/FIXME/XXX comments — resolve or file as GitHub issues.
  - Check any file starting with "#" that references removed claims
    (e.g., PSL(2,7), holomorphic embedding) and either remove or
    re-caption.
  - Verify `src/codon_topo/__init__.py` exports match what currently exists.

### 3. Consolidate documentation

  - Update README.md:
    * Project status (324 tests passing, claim hierarchy summary)
    * CLI usage examples (once you've added the CLI)
    * Scientific caveats (per gpt-5.4-pro review):
       * "Replicates the PHENOMENON of Freeland-Hurst optimality under
         stricter null, not the 10^-6 magnitude."
       * "tRNA pattern is suggestive (n=4, p=0.0625 trend)."
       * "PSL(2,7) claim dropped, KRAS-Fano failed, holomorphic embedding
         misnamed."
    * Target journal: Journal of Theoretical Biology

  - Create ARCHITECTURE.md with:
    * Module dependency graph
    * Data flow (NCBI tables → encoding → scoring → reports)
    * Null-model taxonomy (A, B, C, C_extended, freeland_hurst,
      class_size)
    * How claim_hierarchy.py relates to reports

  - Update CLAUDE.md to point to claim_hierarchy.py as the source of
    truth for "what the paper claims," and list the CLI commands.

### 4. Naming and consistency

Go through every module and fix:
  - Functions named like "foo_weighted" should have a docstring explaining
    what they weight and why.
  - Boolean flags: prefer positive names (`include_stops`) over negative
    (`exclude_stops`); we're consistent now — verify this.
  - Seed constant: export `DEFAULT_SEED = 135325` from `__init__.py` and
    have all Monte Carlo functions default to it.
  - Paper-ready labels: have a constants module or reports module with
    canonical labels like "Freeland-Hurst 1998 block-preserving null"
    that matches what the paper will use.

### 5. Remove old/orphaned files

Check the root of the repository and `codon_topo_preliminary/` directory
for any legacy files that were superseded by the refinement:
  - codon_topo_preliminary/*.py — these were the initial scripts;
    are they still reachable from tests? If not, archive them in a
    `legacy/` folder or delete.
  - Any output/*.csv files that are regenerable via the CLI — commit
    one canonical run and document how to regenerate via `codon-topo all`.

### 6. Update pyproject.toml

```toml
[project]
name = "codon-topo"
version = "0.2.0"   # bump from 0.1.0 after major refactor
dependencies = [
  "numpy>=1.24",
  "scipy>=1.10",
  "requests>=2.28",
  "click>=8.1",    # or typer>=0.9
  "rich>=13.0",    # optional pretty output
]

[project.scripts]
codon-topo = "codon_topo.cli:main"
```

### 7. Final verification

  - `python3.11 -m pytest tests/ -q` → 324 passed (or more if you added
    CLI tests)
  - `python3.11 -m pytest --cov=codon_topo --cov-report=term-missing`
    → coverage should stay ≥ 96%
  - `codon-topo --help` → renders
  - `codon-topo all --output-dir=./output` → produces all CSVs/JSONs
    in under 30 seconds
  - `codon-topo claims` → prints claim hierarchy matching
    reports/claim_hierarchy.py

### 8. Commit strategy

Suggest one atomic commit per task. Do NOT rush to commit everything at
once — smaller commits are easier to review if regressions appear.

Don't create any commits until the entire cleanup passes tests. Then:

  git add -p   (review each change)
  git commit -m "refactor: consolidate CLI, purge dead code, update docs"

After each commit, re-run the test suite.

## DO NOT

- Do not change the numerical results. If a refactor changes a p-value
  or a count, that's a regression — roll back.
- Do not remove negative results (KRAS failure, depth calibration
  non-result). These are scientifically important to keep.
- Do not delete claim_hierarchy.py or its tests — those are the
  structural anchors for the paper's scope.
- Do not rename public API functions without an alias + deprecation
  warning.
- Do not change the default encoding.
- Do not change the default seed (135325).

## EXPECTED FINAL STATE

```
codon-topo v0.2.0
├── CLI: codon-topo {filtration,disconnections,coloring,bit-bias,
│         trna,kras,all,claims}
├── Tests: 324+ passing, coverage ≥ 96%
├── Docs: README.md, ARCHITECTURE.md, CLAUDE.md all current
├── Claim hierarchy: 1 supported, 1 suggestive, 2 exploratory, 6 rejected
└── Zero dead code / stale TODOs / deprecated aliases / misleading comments
```

When done, print a summary of every change made and the final test
count. I will review.
```

---

## For the user running this

Before pasting the prompt into the new session, do ONE thing:

1. Ensure you have a clean working tree (commit or stash current changes).
2. Optionally create a new branch: `git checkout -b cleanup-session`

That way if the cleanup session goes badly you can `git reset --hard` back.

The cleanup session will touch many files. Expect it to take ~30-60 minutes
of agent runtime.
