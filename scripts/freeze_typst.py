"""Freeze a Typst manuscript for Pandoc->LaTeX conversion.

Pandoc's typst-hs reader chokes on:
  - `context { ... }` blocks (counter access)
  - `calc.pow(10.0, negative_int)` (negative exponents)
  - some `#show` rules with closures

So we evaluate the dynamic content in Python BEFORE handing the file to
Pandoc. We precompute every `#fmtk(...)`, `#sci(...)`, `#str(calc.round(...))`,
`#calc.round(...)`, and bare `#var.path` access against the loaded JSON,
and substitute literal text into the source.

Removes:
  - `#set page(...)` / `#set text(...)` / `#set par(...)` / `#set heading(...)`
  - `#set figure(...)` / `#set list(...)` / `#set math.equation(...)`
  - all `#show ...:` rules
  - `#let stats = json("...")` and the alias chain
  - the `#let fmtk(...)` and `#let sci(...)` helper definitions
  - `#bibliography(...)` directive (elsarticle template injects this)

Keeps and resolves:
  - all `#fmtk(EXPR)` -> "10,000"
  - all `#sci(EXPR)` -> "$3.14 times 10^(-8)$" (Typst math mode preserved)
  - `#str(calc.round(EXPR, digits: N))` -> "0.05"
  - `#str(int(calc.round(EXPR, digits: 0)))` -> "10000"
  - `#str(EXPR)` -> string
  - `#calc.round(EXPR, digits: N)` -> rounded number
  - bare `#var.path` -> literal value
  - `#stats._version`, `#pt.n_tables`, etc.

Keeps untouched:
  - all body content, math, tables, figures, citations, lists
  - inline `#let _exp_breaks = ...` are evaluated and substituted at the
    point of every use site
"""

from __future__ import annotations

import json
import math
import re
import sys
from pathlib import Path
from typing import Any


# ---------------------------------------------------------------------------
# 1. Top-level constructs to strip


SET_PREFIXES_TO_STRIP = (
    "#set page(",
    "#set text(",
    "#set par(",
    "#set heading(",
    "#set figure(",
    "#set list(",
    "#set math.equation(",
)


def _consume_balanced(src: str, start: int, open_c: str, close_c: str) -> int:
    """Return index just past the matched closing delimiter at src[start]."""
    assert src[start] == open_c
    depth = 0
    i = start
    in_str = False
    str_quote = ""
    while i < len(src):
        c = src[i]
        if in_str:
            if c == "\\" and i + 1 < len(src):
                i += 2
                continue
            if c == str_quote:
                in_str = False
            i += 1
            continue
        if c in ('"', "'"):
            in_str = True
            str_quote = c
        elif c == open_c:
            depth += 1
        elif c == close_c:
            depth -= 1
            if depth == 0:
                return i + 1
        i += 1
    raise ValueError(f"unbalanced {open_c}{close_c} starting at {start}")


def strip_set_directives(src: str) -> str:
    out = []
    i = 0
    while i < len(src):
        matched = False
        if src[i] == "#":
            for prefix in SET_PREFIXES_TO_STRIP:
                if src.startswith(prefix, i):
                    paren_idx = i + len(prefix) - 1
                    end = _consume_balanced(src, paren_idx, "(", ")")
                    if end < len(src) and src[end] == "\n":
                        end += 1
                    matched = True
                    i = end
                    break
        if not matched:
            out.append(src[i])
            i += 1
    return "".join(out)


def strip_show_rules(src: str) -> str:
    out = []
    i = 0
    n = len(src)
    while i < n:
        if src.startswith("#show ", i):
            colon = src.find(":", i)
            if colon == -1:
                out.append(src[i])
                i += 1
                continue
            j = colon + 1
            while j < n and src[j] not in ("{", "\n"):
                j += 1
            if j < n and src[j] == "{":
                end = _consume_balanced(src, j, "{", "}")
                if end < n and src[end] == "\n":
                    end += 1
                i = end
            else:
                end = j
                if end < n and src[end] == "\n":
                    end += 1
                i = end
        else:
            out.append(src[i])
            i += 1
    return "".join(out)


def strip_bibliography_directive(src: str) -> str:
    return re.sub(r"^#bibliography\([^\n]*\n", "", src, flags=re.MULTILINE)


def strip_set_par_first_line_indent(src: str) -> str:
    return re.sub(r"^#set par\([^)]*\)\n", "", src, flags=re.MULTILINE)


def strip_let_helpers_and_aliases(src: str) -> str:
    """Remove only the multi-line helper function definitions (fmtk, sci);
    the JSON load and alias chain (`#let mm = stats.metrics`, etc.) are
    left in place so typst-hs can resolve any residual content-mode
    conditionals like `#if mi.p < 0.001 [< 0.001] else [...]` that our
    Python evaluator doesn't substitute.
    """
    # Keep `#let fmtk(n) = ...` so typst-hs can evaluate residual calls
    # inside `#for` loops where the loop variable is unknown to our Python
    # evaluator. `fmtk` has no negative-exponent / calc.pow issues.
    # Replace `#let sci(n, sig: 2) = ...` with a typst-hs-friendly version
    # that avoids `calc.pow(10.0, negative_int)` (which typst-hs rejects).
    src = _replace_sci_helper(src)
    return src


_TYPST_HS_FRIENDLY_SCI = """\
#let sci(n, sig: 2) = {
  if n == 0 {
    [0]
  } else if calc.abs(n) >= 0.001 and calc.abs(n) < 1 {
    str(calc.round(n, digits: sig + 1))
  } else {
    let exp = int(calc.floor(calc.log(calc.abs(n), base: 10)))
    let abs_exp = calc.abs(exp)
    let scale = calc.pow(10.0, abs_exp)
    let mant = if exp >= 0 {
      calc.round(n / scale, digits: sig)
    } else {
      calc.round(n * scale, digits: sig)
    }
    $#mant times 10^(#exp)$
  }
}
"""


def _replace_sci_helper(src: str) -> str:
    """Strip the original `#let sci(n, sig: 2) = { ... }` block and reinsert
    a typst-hs-friendly equivalent at the same location."""
    prefix = "#let sci(n, sig: 2) = "
    idx = src.find(prefix)
    if idx == -1:
        return src
    brace = src.find("{", idx)
    if brace == -1:
        return src
    end = _consume_balanced(src, brace, "{", "}")
    if end < len(src) and src[end] == "\n":
        end += 1
    return src[:idx] + _TYPST_HS_FRIENDLY_SCI + src[end:]


def _strip_let_block(src: str, prefix: str) -> str:
    idx = src.find(prefix)
    if idx == -1:
        return src
    brace = src.find("{", idx)
    if brace == -1:
        return src
    end = _consume_balanced(src, brace, "{", "}")
    if end < len(src) and src[end] == "\n":
        end += 1
    return src[:idx] + src[end:]


# ---------------------------------------------------------------------------
# 2. Expression evaluator -- a tiny subset of Typst expressions


class AttrDict(dict):
    """dict with attribute access; unknown keys raise KeyError."""

    def __getattr__(self, name: str) -> Any:
        try:
            return self[name]
        except KeyError as e:
            raise AttributeError(name) from e


def _to_attr(o: Any) -> Any:
    if isinstance(o, dict):
        return AttrDict({k: _to_attr(v) for k, v in o.items()})
    if isinstance(o, list):
        return [_to_attr(v) for v in o]
    return o


# Token types -- order matters: longer comparison operators first so e.g.
# `<=` is one token, not `<` then `=`.
_TOKEN_RE = re.compile(
    r"""
    \s+                                |  # whitespace
    (?P<NUM>-?\d+(?:\.\d+)?(?:e-?\d+)?)|
    (?P<STR>"(?:[^"\\]|\\.)*")         |
    (?P<CMP><=|>=|==|!=|<|>)           |
    (?P<IDENT>[A-Za-z_][A-Za-z_0-9]*)  |
    (?P<DOT>\.)                        |
    (?P<COMMA>,)                       |
    (?P<COLON>:)                       |
    (?P<LPAREN>\()                     |
    (?P<RPAREN>\))                     |
    (?P<LBRACE>\{)                     |
    (?P<RBRACE>\})                     |
    (?P<LBRACK>\[)                     |
    (?P<RBRACK>\])                     |
    (?P<OP>[+\-*/])
    """,
    re.VERBOSE,
)


class Parser:
    """Recursive-descent parser for Typst expressions used in our manuscript."""

    def __init__(self, src: str, env: dict[str, Any]):
        self.src = src
        self.env = env
        self.tokens: list[tuple[str, str]] = []
        self.pos = 0
        self._tokenize(src)

    def _tokenize(self, src: str) -> None:
        i = 0
        while i < len(src):
            m = _TOKEN_RE.match(src, i)
            if not m:
                raise ValueError(f"tokenize failed at {i}: {src[i : i + 20]!r}")
            kind = m.lastgroup
            if kind is not None:
                self.tokens.append((kind, m.group(kind)))
            i = m.end()

    def _peek(self) -> tuple[str, str] | None:
        return self.tokens[self.pos] if self.pos < len(self.tokens) else None

    def _eat(self, kind: str) -> str:
        t = self._peek()
        if t is None or t[0] != kind:
            raise ValueError(f"expected {kind}, got {t}")
        self.pos += 1
        return t[1]

    def _try_eat(self, kind: str, value: str | None = None) -> bool:
        t = self._peek()
        if t is None or t[0] != kind:
            return False
        if value is not None and t[1] != value:
            return False
        self.pos += 1
        return True

    def parse(self) -> Any:
        result = self._expr()
        if self.pos != len(self.tokens):
            raise ValueError(
                f"trailing tokens after {result!r}: {self.tokens[self.pos :]}"
            )
        return result

    # expr := compare ((</<=/>/>=/==/!=) compare)?
    def _expr(self) -> Any:
        v = self._compare()
        t = self._peek()
        if t is None or t[0] != "CMP":
            return v
        op = t[1]
        self.pos += 1
        rhs = self._compare()
        if op == "<":
            return v < rhs
        if op == "<=":
            return v <= rhs
        if op == ">":
            return v > rhs
        if op == ">=":
            return v >= rhs
        if op == "==":
            return v == rhs
        if op == "!=":
            return v != rhs
        raise ValueError(f"unknown comparison {op}")

    # compare := term ((+|-) term)*
    def _compare(self) -> Any:
        v = self._term()
        while True:
            t = self._peek()
            if t is None or t[0] != "OP" or t[1] not in ("+", "-"):
                break
            op = t[1]
            self.pos += 1
            rhs = self._term()
            v = (v + rhs) if op == "+" else (v - rhs)
        return v

    # term := factor ((*|/) factor)*
    def _term(self) -> Any:
        v = self._factor()
        while True:
            t = self._peek()
            if t is None or t[0] != "OP" or t[1] not in ("*", "/"):
                break
            op = t[1]
            self.pos += 1
            rhs = self._factor()
            v = (v * rhs) if op == "*" else (v / rhs)
        return v

    # factor := atom (.IDENT(args)? | .IDENT)*
    def _factor(self) -> Any:
        v = self._atom()
        while True:
            t = self._peek()
            if t is None:
                break
            if t[0] == "DOT":
                self.pos += 1
                name = self._eat("IDENT")
                # method call?
                if self._try_eat("LPAREN"):
                    args, kwargs = self._args()
                    self._eat("RPAREN")
                    v = self._method(v, name, args, kwargs)
                else:
                    if isinstance(v, dict):
                        v = v[name]
                    elif hasattr(v, name):
                        v = getattr(v, name)
                    else:
                        raise KeyError(name)
            else:
                break
        return v

    def _atom(self) -> Any:
        t = self._peek()
        if t is None:
            raise ValueError("unexpected end of expression")
        kind, val = t
        if kind == "NUM":
            self.pos += 1
            return float(val) if "." in val or "e" in val else int(val)
        if kind == "STR":
            self.pos += 1
            return val[1:-1]  # strip quotes
        if kind == "OP" and val == "-":
            self.pos += 1
            inner = self._factor()
            return -inner
        if kind == "LPAREN":
            self.pos += 1
            # Typst empty array literal `()`
            if self._try_eat("RPAREN"):
                return []
            # Typst empty dict literal `(:)`
            if self._try_eat("COLON"):
                self._eat("RPAREN")
                return {}
            # Typst dict literal `(key: value, ...)` -- distinguish from
            # parenthesized expression by peeking past IDENT for COLON.
            tok = self._peek()
            if tok and tok[0] == "IDENT":
                next_tok = (
                    self.tokens[self.pos + 1]
                    if self.pos + 1 < len(self.tokens)
                    else None
                )
                if next_tok and next_tok[0] == "COLON":
                    return self._dict_literal()
            v = self._expr()
            self._eat("RPAREN")
            return v
        if kind == "IDENT":
            self.pos += 1
            # function call?
            if self._try_eat("LPAREN"):
                args, kwargs = self._args()
                self._eat("RPAREN")
                return self._call_func(val, args, kwargs)
            # bare identifier -- look up in env
            if val in self.env:
                return self.env[val]
            raise KeyError(f"unknown identifier {val!r}")
        raise ValueError(f"unexpected token {t!r}")

    def _args(self) -> tuple[list[Any], dict[str, Any]]:
        args: list[Any] = []
        kwargs: dict[str, Any] = {}
        if self._peek() and self._peek()[0] == "RPAREN":
            return args, kwargs
        while True:
            # kwarg?  IDENT COLON expr
            saved = self.pos
            tok = self._peek()
            if tok and tok[0] == "IDENT":
                # peek ahead for ":"
                next_tok = (
                    self.tokens[self.pos + 1]
                    if self.pos + 1 < len(self.tokens)
                    else None
                )
                if next_tok and next_tok[0] == "COLON":
                    name = tok[1]
                    self.pos += 2  # skip IDENT and COLON
                    value = self._kwarg_value()
                    kwargs[name] = value
                    if not self._try_eat("COMMA"):
                        break
                    continue
            self.pos = saved
            args.append(self._expr())
            if not self._try_eat("COMMA"):
                break
        return args, kwargs

    def _kwarg_value(self) -> Any:
        """Kwarg value is any expression (incl. `(:)` empty-dict literal,
        which `_atom()` already handles)."""
        return self._expr()

    def _dict_literal(self) -> dict[str, Any]:
        """Parse `key: value, key2: value2, ...)`. The leading `(` was
        already consumed by `_atom()`; we consume up through the closing `)`.
        """
        result: dict[str, Any] = {}
        while True:
            tok = self._peek()
            if tok is None:
                raise ValueError("unterminated dict literal")
            if tok[0] == "RPAREN":
                self.pos += 1
                return result
            key = self._eat("IDENT")
            self._eat("COLON")
            value = self._expr()
            result[key] = value
            if self._try_eat("COMMA"):
                continue
            self._eat("RPAREN")
            return result

    def _call_func(self, name: str, args: list[Any], kwargs: dict[str, Any]) -> Any:
        if name == "calc":
            raise ValueError("bare 'calc' identifier; use calc.round / calc.abs / etc.")
        if name == "int":
            return int(args[0])
        if name == "str":
            v = args[0]
            return _typst_str(v)
        if name == "float":
            return float(args[0])
        raise ValueError(f"unknown function {name}")

    def _method(
        self, recv: Any, name: str, args: list[Any], kwargs: dict[str, Any]
    ) -> Any:
        # `calc.X(...)` -- recv is the symbol `calc`
        if recv is _CALC:
            return _calc_dispatch(name, args, kwargs)
        # `<dict>.at("key", default: D)` or `<list>.at(idx, default: D)`
        if name == "at":
            key = args[0]
            default = kwargs.get("default")
            try:
                if isinstance(recv, list):
                    return recv[int(key)]
                if isinstance(recv, dict):
                    return recv[key]
                return getattr(recv, key)
            except (KeyError, IndexError, AttributeError):
                if "default" in kwargs:
                    return default
                raise
        if name == "len":
            return len(recv)
        raise ValueError(f"unknown method {name} on {type(recv).__name__}")


class _CalcSentinel:
    """Marker for the `calc` namespace identifier."""

    pass


_CALC = _CalcSentinel()


def _calc_dispatch(name: str, args: list[Any], kwargs: dict[str, Any]) -> Any:
    digits = kwargs.get("digits", 0)
    base = kwargs.get("base")
    if name == "round":
        return round(args[0], int(digits))
    if name == "abs":
        return abs(args[0])
    if name == "min":
        return min(args)
    if name == "max":
        return max(args)
    if name == "floor":
        return math.floor(args[0])
    if name == "ceil":
        return math.ceil(args[0])
    if name == "log":
        if base is not None:
            return math.log(args[0], float(base))
        return math.log(args[0])
    if name == "pow":
        return math.pow(args[0], args[1])
    if name == "rem":
        return args[0] % args[1]
    raise ValueError(f"unknown calc.{name}")


def _typst_str(v: Any) -> str:
    """Mimic Typst's str() conversion."""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        # Typst prints floats without trailing zeros if integer-valued,
        # otherwise with full precision. round() should have made these tidy.
        if v == int(v):
            return str(int(v))
        return repr(v)
    return str(v)


def evaluate(expr: str, env: dict[str, Any]) -> Any:
    return Parser(expr, env).parse()


# ---------------------------------------------------------------------------
# 3. Evaluator helpers for the manuscript's specific call patterns


def _fmtk(n: float) -> str:
    """Mirror of the Typst #fmtk helper: round to int, comma-separate thousands."""
    v = int(round(float(n)))
    return f"{v:,}"


def _sci(n: float, sig: int = 2) -> str:
    """Mirror of Typst #sci: emit a Typst math fragment.

    The original helper:
      - n == 0:                     "0"
      - 0.001 <= |n| < 1:           round(n, sig+1) as a decimal string
      - else:                       `$<mant> times 10^(<exp>)$` (math mode)
    """
    if n == 0:
        return "0"
    a = abs(n)
    if 0.001 <= a < 1:
        return _typst_str(round(n, sig + 1))
    exp = int(math.floor(math.log10(a)))
    mant = round(n / (10**exp), sig)
    mant_s = _typst_str(mant)
    return f"$#({mant_s}) times 10^({exp})$"


# ---------------------------------------------------------------------------
# 4. Substitution pass over the source


# Patterns to target -- order matters because some are specializations of others.
# Each is a (compiled_regex, handler_name) pair where the regex captures the
# beginning of the call and we use balanced-paren consumption to find the end.

CALL_PREFIXES = (
    "#fmtk(",
    "#sci(",
    "#str(",
    "#calc.round(",
    "#calc.min(",
    "#calc.max(",
    "#calc.abs(",
    "#calc.floor(",
)


def _find_call(src: str, start: int) -> tuple[str, int, int] | None:
    """Find the next call site at or after `start`. Returns (prefix, start_idx, paren_idx)."""
    n = len(src)
    i = start
    while i < n:
        if src[i] != "#":
            i += 1
            continue
        for prefix in CALL_PREFIXES:
            if src.startswith(prefix, i):
                paren_idx = i + len(prefix) - 1
                return prefix, i, paren_idx
        i += 1
    return None


def _evaluate_call(prefix: str, inner: str, env: dict[str, Any]) -> str:
    """Evaluate a call given its prefix and the text inside the outer parens."""
    if prefix == "#fmtk(":
        v = evaluate(inner, env)
        return _fmtk(v)
    if prefix == "#sci(":
        # Possible signatures: sci(n) or sci(n, sig: K)
        # Parse arguments using the parser's _args via a small wrapper.
        p = Parser("(" + inner + ")", env)
        # Re-tokenize via a pseudo-call
        # Easier: split by top-level comma manually.
        parts = _split_top_level_comma(inner)
        n = evaluate(parts[0], env)
        sig = 2
        for p in parts[1:]:
            if ":" in p:
                key, val = p.split(":", 1)
                if key.strip() == "sig":
                    sig = int(evaluate(val.strip(), env))
        return _sci(n, sig=sig)
    if prefix == "#str(":
        v = evaluate(inner, env)
        return _typst_str(v)
    if prefix == "#calc.round(":
        # calc.round(EXPR, digits: N)
        parts = _split_top_level_comma(inner)
        n = evaluate(parts[0], env)
        digits = 0
        for p in parts[1:]:
            if ":" in p:
                key, val = p.split(":", 1)
                if key.strip() == "digits":
                    digits = int(evaluate(val.strip(), env))
        return _typst_str(round(n, digits))
    if prefix == "#calc.min(":
        parts = _split_top_level_comma(inner)
        return _typst_str(min(evaluate(p, env) for p in parts))
    if prefix == "#calc.max(":
        parts = _split_top_level_comma(inner)
        return _typst_str(max(evaluate(p, env) for p in parts))
    if prefix == "#calc.abs(":
        return _typst_str(abs(evaluate(inner, env)))
    if prefix == "#calc.floor(":
        return _typst_str(math.floor(evaluate(inner, env)))
    raise ValueError(f"unknown prefix {prefix}")


def _split_top_level_comma(src: str) -> list[str]:
    """Split a string on top-level commas (ignoring commas inside parens/brackets/strings)."""
    out: list[str] = []
    depth = 0
    current = []
    in_str = False
    str_quote = ""
    for c in src:
        if in_str:
            if c == str_quote:
                in_str = False
            current.append(c)
            continue
        if c in ('"', "'"):
            in_str = True
            str_quote = c
            current.append(c)
        elif c in "([{":
            depth += 1
            current.append(c)
        elif c in ")]}":
            depth -= 1
            current.append(c)
        elif c == "," and depth == 0:
            out.append("".join(current))
            current = []
        else:
            current.append(c)
    if current:
        out.append("".join(current))
    return out


def substitute_calls(src: str, env: dict[str, Any]) -> str:
    """Walk the source, substituting every supported call expression."""
    out: list[str] = []
    i = 0
    while True:
        match = _find_call(src, i)
        if match is None:
            out.append(src[i:])
            break
        prefix, start, paren_idx = match
        out.append(src[i:start])
        end = _consume_balanced(src, paren_idx, "(", ")")
        inner = src[paren_idx + 1 : end - 1]
        try:
            replacement = _evaluate_call(prefix, inner, env)
        except Exception as e:
            print(
                f"WARN: failed to evaluate {prefix}{inner}): {e}; leaving as-is",
                file=sys.stderr,
            )
            out.append(src[start:end])
        else:
            out.append(replacement)
        i = end
    return "".join(out)


# ---------------------------------------------------------------------------
# 5. Inline `#let` resolution and bare `#var.path` substitution


_INLINE_LET_RE = re.compile(
    r"^#let\s+([A-Za-z_][A-Za-z_0-9]*)\s*=\s*(.+)$", re.MULTILINE
)


def resolve_inline_lets(src: str, env: dict[str, Any]) -> str:
    """For every top-level `#let NAME = EXPR` whose EXPR can be evaluated in
    the current env, evaluate it and bind NAME in env. The line is left in
    the source so that typst-hs (Pandoc's evaluator) can also resolve it
    for any residual expressions our Python evaluator does not handle
    (e.g. `#if mi.p < 0.001 [...]` content-mode conditionals).

    Multiple passes handle let-chains where one binding depends on another.
    """
    while True:
        progress = False
        for line in src.split("\n"):
            m = _INLINE_LET_RE.match(line)
            if not m:
                continue
            rhs = m.group(2).strip()
            if rhs.endswith("{") or rhs == "{":
                continue
            name = m.group(1)
            if name in env:
                continue
            try:
                value = evaluate(rhs, env)
            except Exception:
                continue
            env[name] = value
            progress = True
        if not progress:
            break
    return src


# Bare references like `#tk43.observed_breaks`, `#stats._version`, `#pos2_pct`.
# We match `#IDENT(.IDENT)*` not followed by `(` (function call) or `[` (content).
_BARE_REF_RE = re.compile(
    r"#([A-Za-z_][A-Za-z_0-9]*(?:\.[A-Za-z_][A-Za-z_0-9]*)*)(?![\(\[A-Za-z_0-9])"
)


# Names we deliberately leave alone -- they're Typst content/markup, not data.
_BARE_REF_SKIP = {
    "v",
    "h",
    "let",
    "set",
    "show",
    "context",
    "text",
    "align",
    "block",
    "table",
    "figure",
    "cite",
    "ref",
    "image",
    "pagebreak",
    "heading",
    "super",
    "sub",
    "sym",
    "calc",
    "metadata",
    "smallcaps",
    "rect",
    "circle",
    "place",
    "stack",
    "grid",
    "linebreak",
    "if",
    "for",
    "while",
    "import",
    "include",
    "underline",
    "emph",
    "strong",
    "raw",
    "lorem",
    "footnote",
    "highlight",
    "list",
    "enum",
    "outline",
    "bibliography",
    "page",
    "par",
    "math",
}


def substitute_bare_refs(src: str, env: dict[str, Any]) -> str:
    """Replace `#var.path` with its literal value when var is in env."""

    def repl(m: re.Match) -> str:
        path = m.group(1)
        head = path.split(".")[0]
        if head in _BARE_REF_SKIP or head not in env:
            return m.group(0)
        try:
            value = evaluate(path, env)
        except Exception:
            return m.group(0)
        return _typst_str(value)

    return _BARE_REF_RE.sub(repl, src)


# ---------------------------------------------------------------------------
# 5b. Expand `..LIST.map(VAR => ( cell, cell, ... )).flatten()?,?` spreads.


_MAP_SPREAD_RE = re.compile(
    r"\.\.([A-Za-z_][A-Za-z_0-9]*(?:\.[A-Za-z_][A-Za-z_0-9]*)*)\.map\(\s*"
    r"([A-Za-z_][A-Za-z_0-9]*)\s*=>\s*"
    r"\("  # opening paren of the tuple body
)


def expand_map_spreads(src: str, env: dict[str, Any]) -> str:
    """Expand `..LIST.map(VAR => (CELL, CELL, ...)).flatten()?,?` table-row
    spreads. The tuple body is repeated once per item in LIST with VAR
    bound to that item; substitute_calls + substitute_bare_refs are then
    run on each instance to resolve `#VAR.path`, `#fmtk(VAR.X)`, etc.

    Pandoc's typst-hs cannot evaluate method chains like `r.x.map(str).join(...)`,
    so we materialise the rows in Python before handing the file to Pandoc.
    """
    out: list[str] = []
    i = 0
    n = len(src)
    while i < n:
        m = _MAP_SPREAD_RE.search(src, i)
        if not m:
            out.append(src[i:])
            break
        out.append(src[i : m.start()])
        list_path = m.group(1)
        var = m.group(2)
        # Find the tuple body's closing paren. The structure is:
        #   .map(VAR => (cell, cell))[.flatten()]
        #               ^body_start ^body_end
        #                            ^^need to also consume the closing
        #                              `)` of the map( call.
        body_start = m.end() - 1  # index of `(`  (tuple-body open)
        body_end = _consume_balanced(src, body_start, "(", ")")
        body = src[body_start + 1 : body_end - 1]

        j = body_end
        # Consume the closing `)` of `.map(`.
        if j < n and src[j] == ")":
            j += 1
        # Optionally consume `.flatten()` and trailing comma.
        if src.startswith(".flatten()", j):
            j += len(".flatten()")

        # Look up the list in env.
        try:
            items = evaluate(list_path, env)
        except Exception as e:
            print(
                f"WARN: cannot expand map spread on {list_path}: {e}; leaving as-is",
                file=sys.stderr,
            )
            out.append(src[m.start() : j])
            i = j
            continue
        if not isinstance(items, list):
            print(
                f"WARN: {list_path} is not a list; leaving as-is",
                file=sys.stderr,
            )
            out.append(src[m.start() : j])
            i = j
            continue

        # Split body by top-level commas to get individual cells.
        cells = _split_top_level_comma(body)
        # Strip whitespace; drop trailing empty (from trailing comma).
        cells = [c.strip() for c in cells if c.strip()]

        rendered_rows: list[str] = []
        for item in items:
            local_env = dict(env)
            local_env[var] = item
            rendered_cells = []
            for cell in cells:
                expanded = _resolve_cell(cell, local_env)
                rendered_cells.append(expanded)
            rendered_rows.append(", ".join(rendered_cells))

        replacement = ", ".join(rendered_rows)
        # Empty list -> empty expansion. Left as-is, this leaves a lone `,`
        # in the enclosing argument list (`table(header(...),  ,)`) which
        # pandoc's typst-hs rejects. Swallow the trailing `,` after
        # `.flatten()` (and any preceding whitespace we already emitted)
        # so the enclosing list stays well-formed. Non-empty expansions
        # keep the existing behavior -- the spread becomes real cells and
        # the trailing comma is desirable.
        if not rendered_rows:
            if j < n and src[j] == ",":
                j += 1
            # Trim trailing whitespace on the last emitted chunk so we
            # don't leave a dangling indented blank line.
            if out and out[-1] and out[-1][-1] in " \t":
                out[-1] = out[-1].rstrip(" \t")
        out.append(replacement)
        i = j
    return "".join(out)


def _resolve_cell(cell: str, env: dict[str, Any]) -> str:
    """Resolve all dynamic references inside a cell of an expanded map spread.

    Order matters:
      1. user helpers like `_format_clade_label(r.X)` first -- they take
         `r.X` as input and return a literal string
      2. `#r.x.map(str).join(", ")` next -- consumes a path-with-method-chain
         that bare-ref substitution would otherwise prematurely replace
      3. `#fmtk(...)` / `#sci(...)` / `#str(...)` calls (including any
         nested ones inside the branches of a content-mode `#if`)
      4. bare `#r.x` references
      5. content-mode `#if EXPR [true] else [false]` conditionals using
         resolved arguments
      6. `#{EXPR}` block splices
    """
    cell = _resolve_user_helpers(cell, env)
    cell = _resolve_map_str_join(cell, env)
    cell = substitute_calls(cell, env)
    cell = substitute_bare_refs(cell, env)
    cell = _resolve_if_else(cell, env)
    cell = _resolve_block_splice(cell, env)
    return cell


def _resolve_if_else(text: str, env: dict[str, Any]) -> str:
    """Replace `#if EXPR [TRUE_BRANCH] else [FALSE_BRANCH]` with the
    appropriate branch when EXPR can be evaluated. EXPR may use loop
    variables bound in env; the branches stay verbatim (already had their
    inner calls/refs resolved by earlier steps in _resolve_cell).
    """
    out: list[str] = []
    i = 0
    n = len(text)
    while i < n:
        if not text.startswith("#if ", i):
            out.append(text[i])
            i += 1
            continue
        # find the first `[` that isn't preceded by an operator continuation
        bracket = text.find("[", i + 4)
        if bracket == -1:
            out.append(text[i])
            i += 1
            continue
        cond_src = text[i + 4 : bracket].strip()
        try:
            cond_value = evaluate(cond_src, env)
        except Exception:
            out.append(text[i])
            i += 1
            continue
        true_end = _consume_balanced(text, bracket, "[", "]")
        true_branch = text[bracket + 1 : true_end - 1]
        # expect ` else [...]`
        rest = text[true_end:]
        m = re.match(r"\s*else\s*\[", rest)
        if not m:
            # no else branch; emit true if cond_value, otherwise empty
            chosen = true_branch if cond_value else ""
            out.append(chosen)
            i = true_end
            continue
        false_open = true_end + m.end() - 1
        false_end = _consume_balanced(text, false_open, "[", "]")
        false_branch = text[false_open + 1 : false_end - 1]
        chosen = true_branch if cond_value else false_branch
        out.append(chosen)
        i = false_end
    return "".join(out)


# Hand-coded emulators for user-defined helpers in the manuscript/supplement.
# Each takes the resolved Python argument(s) and returns the literal text
# that Typst would have produced.


def _format_clade_label(s: str) -> str:
    """Mirror of `#let _format_clade_label(s) = parts.slice(1).join(" ")`
    (drops the first underscore-separated token)."""
    parts = s.replace("_", " ").split(" ")
    return " ".join(parts[1:])


def _phylo_label(s: str) -> str:
    """Same shape as `_format_clade_label` -- supplement's phylo regimes."""
    parts = s.replace("_", " ").split(" ")
    return " ".join(parts[1:])


def _abbrev_species(name: str) -> str:
    """Python mirror of the Typst `abbrev_species()` helper: "Saccharomyces
    cerevisiae" -> italicised "S. cerevisiae"; strips trailing "_mito" and
    "_nuclear" suffixes; leaves single-word names in italics as-is."""
    clean = name.replace("_mito", "").replace("_nuclear", "").strip()
    parts = clean.split(" ")
    if len(parts) >= 2 and parts[0]:
        body = f"{parts[0][0]}. {' '.join(parts[1:])}"
    else:
        body = clean
    # Typst emph in the pandoc conversion path -> LaTeX \emph{...}. We emit
    # the Typst source form so the pandoc reader sees emph consistently.
    return f"#emph[{body}]"


_USER_HELPERS = {
    "_format_clade_label": _format_clade_label,
    "_phylo_label": _phylo_label,
    "abbrev_species": _abbrev_species,
}


_USER_HELPER_RE = re.compile(
    r"#(" + "|".join(re.escape(name) for name in _USER_HELPERS) + r")\(([^()]*)\)"
)


def _resolve_user_helpers(text: str, env: dict[str, Any]) -> str:
    def repl(m: re.Match) -> str:
        name = m.group(1)
        arg_src = m.group(2).strip()
        try:
            arg = evaluate(arg_src, env)
        except Exception:
            return m.group(0)
        try:
            result = _USER_HELPERS[name](arg)
        except Exception:
            return m.group(0)
        return _typst_str(result)

    return _USER_HELPER_RE.sub(repl, text)


_MAP_STR_JOIN_RE = re.compile(
    r"#([A-Za-z_][A-Za-z_0-9]*(?:\.[A-Za-z_][A-Za-z_0-9]*)*)\.map\(str\)\.join\(\"([^\"]*)\"\)"
)


def _resolve_map_str_join(text: str, env: dict[str, Any]) -> str:
    """Replace `#path.map(str).join("SEP")` with the string-joined value."""

    def repl(m: re.Match) -> str:
        path = m.group(1)
        sep = m.group(2)
        try:
            value = evaluate(path, env)
        except Exception:
            return m.group(0)
        if not isinstance(value, list):
            return m.group(0)
        return sep.join(_typst_str(v) for v in value)

    return _MAP_STR_JOIN_RE.sub(repl, text)


_BLOCK_SPLICE_RE = re.compile(r"#\{([^{}]+)\}")


def _resolve_block_splice(text: str, env: dict[str, Any]) -> str:
    """Replace `#{EXPR}` with the evaluated expression (no nested braces)."""

    def repl(m: re.Match) -> str:
        expr = m.group(1).strip()
        try:
            value = evaluate(expr, env)
        except Exception:
            return m.group(0)
        return _typst_str(value)

    return _BLOCK_SPLICE_RE.sub(repl, text)


# ---------------------------------------------------------------------------
# 6. End-to-end driver


def freeze(src: str, stats_path: Path) -> str:
    src = strip_set_directives(src)
    src = strip_show_rules(src)
    src = strip_set_par_first_line_indent(src)
    src = strip_bibliography_directive(src)
    src = strip_let_helpers_and_aliases(src)

    # Build evaluation env from JSON.
    raw = json.loads(stats_path.read_text())
    stats = _to_attr(raw)
    env: dict[str, Any] = {
        "calc": _CALC,
        "stats": stats,
        "mm": stats.metrics,
        "tq6": stats.topology_avoidance_q6,
        "tk43": stats.topology_avoidance_k43,
        "pt": stats.per_table,
        "rho_data": stats.rho_sweep,
        "dec": stats.decomposition,
        "cl": stats.condlogit,
        "trna": stats.trna,
    }

    # Pre-resolve inline #let lines so calls below see them.
    src = resolve_inline_lets(src, env)
    # Expand `..LIST.map(VAR => (...))` table-row spreads in Python (typst-hs
    # cannot handle method chains like `.map(str).join(...)`).
    src = expand_map_spreads(src, env)
    # Substitute calls (#fmtk, #sci, #str, #calc.*).
    src = substitute_calls(src, env)
    # And bare references (#tk43.observed_breaks, #pos2_pct, ...).
    src = substitute_bare_refs(src, env)
    return src


def main() -> int:
    if len(sys.argv) != 4:
        print(
            f"usage: {sys.argv[0]} <input.typ> <manuscript_stats.json> <output.typ>",
            file=sys.stderr,
        )
        return 2
    in_path = Path(sys.argv[1])
    stats_path = Path(sys.argv[2])
    out_path = Path(sys.argv[3])
    src = in_path.read_text()
    frozen = freeze(src, stats_path)
    out_path.write_text(frozen)
    print(f"frozen {in_path} -> {out_path} ({len(src)} -> {len(frozen)} bytes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
