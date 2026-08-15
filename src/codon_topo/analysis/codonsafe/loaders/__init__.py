"""Dataset loaders for published genome recoding studies.

Each loader module provides a function that reads published supplementary
data and returns a list of CodonSwapEvent objects with RecodingOutcome.

Supplementary data files should be placed in:
  data/codonsafe/<study_name>/

Loader conventions:
  - All codons normalized to RNA alphabet (T->U) via normalize.py
  - Event IDs are globally unique: "{study}_{gene}_{pos}_{target}"
  - Covariates are stored in the event.covariates dict
  - Missing covariates are omitted (not set to NaN)
"""
