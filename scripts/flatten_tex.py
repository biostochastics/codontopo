"""Inline `\\input{...}` bodies into their wrapper to produce a single
self-contained `.tex` file for editorial systems that don't run a
multi-file `latexmk` build (Elsevier Editorial Manager, arXiv-flat, etc.).

Usage:
    python scripts/flatten_tex.py <wrapper.tex> <body_basename> <out.tex>

`<body_basename>` is the argument that appears inside `\\input{...}`. The
tool matches `\\input{<basename>}` or `\\input{<basename>.tex}`.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path


def flatten(wrapper_path: Path, body_basename: str, body_text: str) -> str:
    wrapper = wrapper_path.read_text()
    # Escape regex metacharacters in the basename (only `.` and `_` are
    # realistic in practice, but be robust).
    pat = re.compile(r"\\input\{" + re.escape(body_basename) + r"(?:\.tex)?\}")
    if not pat.search(wrapper):
        raise SystemExit(f"error: no \\input{{{body_basename}}} in {wrapper_path}")
    # str.replace-safe substitution -- avoids re.sub's backslash re-escaping.
    return pat.sub(lambda _m: body_text, wrapper, count=1)


def main() -> int:
    if len(sys.argv) != 4:
        print(
            f"usage: {sys.argv[0]} <wrapper.tex> <body_basename> <out.tex>",
            file=sys.stderr,
        )
        return 2
    wrapper_path = Path(sys.argv[1])
    body_basename = sys.argv[2]
    out_path = Path(sys.argv[3])
    body_text = (wrapper_path.parent / f"{body_basename}.tex").read_text()
    flat = flatten(wrapper_path, body_basename, body_text)
    out_path.write_text(flat)
    print(f"flattened {wrapper_path} + {body_basename}.tex -> {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
