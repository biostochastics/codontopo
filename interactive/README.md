# Logic Engines — Clayworth Proof Bench (React)

`ClayworthProofBench.jsx` is a single-file, self-contained React artifact bundling
**eight verified reasoning engines**, each a faithful JavaScript port of its Python
original (no Pyodide, no network, no build-time dependencies beyond React):

| Tab | Engine | What it checks |
|-----|--------|----------------|
| Core logics | `prover_core` | Propositional ND (truth-table backstop), modal **K/T/S4/S5** (Kripke countermodel search), Fano/octonion units (real Cayley–Dickson associator) |
| Equational algebra | `equational` | Birkhoff equational logic over 5 theories: group, octonion (7 Fano triples as defining equations), modular ℤ/5ℤ, Boolean algebra, prime field GF(7) |
| Linear (MLL) | `linear` | Multiplicative linear logic with backtracking resource search (no weakening / no contraction) |
| Discharge / subproofs | `discharge` | Natural deduction with assumption discharge and scope enforcement; flags genuine theorems |
| Incidence geometry | `incidence` | Finite incidence structures (Fano plane PG(2,2), triangle); projective closure generalizes the Fano rule |
| Finite field GF(8) | `finite_field` | 𝔽₂[x]/(x³+x+1): add (XOR), multiply, inverse, powers of the primitive element, multiplicative order |

A collapsible **self-test panel** runs every engine's checks live in-artifact, so a
broken bundle is immediately visible.

## Usage

This is a standalone component with a default export. Drop it into any React
toolchain (Vite, Next.js, CRA, etc.):

```jsx
import ClayworthProofBench from "./ClayworthProofBench.jsx";

export default function App() {
  return <ClayworthProofBench />;
}
```

It styles itself (dark theme, inline styles) and needs only `react` at runtime.

## Verification

The pure-JS engine layer (everything above the React section) is verified
against the Python originals case-for-case. The bundled `runAllSelfTests()`
harness — 40+ assertions spanning all eight engines — is exercised both in the
UI panel and out-of-band during development; all pass.
