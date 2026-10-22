# sympy_matrix_tools

<!-- Time-stamp: "2026-10-22T18:27:21Z" -->

Tools for symbolic matrix calculus and checked derivations with
[SymPy](https://www.sympy.org/).

  * **Undefined matrix-valued functions** (`MatrixFunction`) with the
    chain rule, e.g. `F(X*Y).diff(X)`, `F(t*X).diff(t)`.
  * **Algebra of non-commutative matrix expressions**: coefficients,
    collection and division (`mat_coeff`, `mat_collect`,
    `mat_divide`).
  * **Sums and products of matrices** (`MatSum`, `MatProduct`) and
    manipulation of symbolic series (`Sum_step_forward`, ...).
  * **A derivation ledger** (`Derivation`) that records why each step of
    a derivation is believed (a computation by SymPy, a rewrite by a
    named fact, or a bare assertion) and checks every step exactly with
    random rational matrices.
  * **Logic tools**: higher-order pattern unification, rewriting, and
    an Isabelle-style prover with first-order logic built on a minimal
    core.

As of SymPy 1.14, SymPy has nothing like `MatrixFunction`,
`MatrixWild` or `MatSum`.

**Status (2026-10):** maintained for SymPy 1.14. See
[CHANGELOG.md](CHANGELOG.md).


## Installation

You need SymPy >= 1.14.

```sh
pip install "sympy>=1.14"
pip install git+https://github.com/JRF-2018/sympy_matrix_tools
```

In Google Colab, put `!` before `pip`. To fix the version, specify a
tag such as `...sympy_matrix_tools@v0.8.0`.


## Quick tour

The chain rule for an undefined matrix-valued function:

```python
>>> from sympy import *
>>> from sympy_matrix_tools import *
>>> x1 = Symbol("x1")
>>> N = Symbol("N", integer=True)
>>> M1 = MatrixSymbol("M1", N, N)
>>> M2 = MatrixSymbol("M2", N, N)
>>> Mf = MatrixFunction("Mf", N, N)
>>> (Mf(2 * x1 + 1) + M1).diff(x1)
2*Subs(Derivative(Mf(_xi_1), _xi_1), _xi_1, 2*x1 + 1) + 0
>>> (Mf(2 * x1 + 1) + M1).diff(x1).subs(Mf, Lambda(x1, x1 * M1 + M2)).doit()
2*M1

```

Coefficients and collection in non-commutative expressions:

```python
>>> z = (M1 + M2) * M1 + 2 * x1 * M1 * M2 + M2 * (M1 ** 2) + 3 * M2 * M1
>>> mat_coeff(z, M1)
M1 + M2
>>> mat_coeff(z, M1, right=True)
2*x1*M2
>>> mat_collect(M1 * M2 + M2 * M2 + 3 * M2, M2, expand_pow=True)
(3*I + M1 + M2)*M2

```

A derivation ledger: each step has a justification, and `check()`
tests every step with random rational matrices.

```python
>>> from sympy_matrix_tools.derivation import Derivation, cas, rule
>>> n = Symbol("n", integer=True, positive=True)
>>> A = MatrixSymbol("A", n, n)
>>> X = MatrixSymbol("X", n, n)
>>> d = Derivation("small example")
>>> d.hypothesis(Eq(A.T, A), "A_sym")
>>> d.define(X, A * A.T, "X_def")
>>> d.calc("X_sq", X,
...        (A * A.T, rule("X_def")),
...        (A * A, rule("A_sym")),
...        (A ** 2, cas("doit")),
...        2 * A)                     # a wrong step without justification
>>> d.check(trials=2)
False
>>> print(d.report())
Derivation: small example
[hypothesis] A_sym: A.T = A
[definition] X_def: X = A*A.T
[derived] X_sq: X = 2*A
  X
  = A*A.T    by rule X_def    symbolic: ok    numeric: ok 2/2
  = A**2    by rule A_sym    symbolic: ok    numeric: ok 2/2
  = A**2    by cas doit    symbolic: ok    numeric: ok 2/2
  = 2*A    by asserted    symbolic: none    numeric: FAILED 0/2
Summary: 4 steps, 1 asserted, 0 symbolic failures, 1 numeric failure, 0 open obligations.

```


## Documentation

| Topic | Document | Main names |
|---|---|---|
| Matrix expressions | [docs/matrix_tools.md](docs/matrix_tools.md) | `mat_coeff`, `mat_collect`, `mat_divide`, `mat_trivial_divide`, `mat_expand`, `partial_apply` |
| Matrix functions and the chain rule | [docs/matrix_function.md](docs/matrix_function.md) | `MatrixFunction`, `mat_diff`, `freeze_matrix_function`, `melt_matrix_function` |
| Sequences and sums of matrices | [docs/sequences.md](docs/sequences.md) | `Sum_expand`, `Sum_collect`, `Sum_step_forward`, `MatSum`, `MatProduct` |
| Wild symbols and application | [docs/wild_apply.md](docs/wild_apply.md) | `MatrixWild`, `MatrixDummy`, `WildMatrixFunction`, `Apply`, `MatApply`, `PredApply` |
| Derivation ledger | [docs/derivation.md](docs/derivation.md) | `Derivation`, `Grad`, `cas`, `rule` |
| Unification, rewriting and proofstates | [docs/logic.md](docs/logic.md) | `ho_unify`, `rewriterules_from_eqs`, `resolve_mp`, `make_proofstate`, `resolve_implications`, `auto` |
| First-order logic | [docs/fol.md](docs/fol.md) | `Theorem`, `parse_formula`, `Exists`, `Bot`, `neg`, `exI`, `exE`, `by_auto` |
| Fixes of SymPy and known SymPy bugs | [docs/fixes.md](docs/fixes.md) | `fix_ExpectationMatrix_expand`, `fix_Basic_MatrixSymbol_about_dummy`, `fix_Transpose_eval_derivative` |
| For AI agents | [docs/ai_guide.md](docs/ai_guide.md), [llms.txt](llms.txt) | when to use, recipes and pitfalls |
| Development | [docs/development.md](docs/development.md) | tests, conventions |

The examples in the documents are tested by doctest.


## Examples (notebooks)

  * [examples/lq_control.ipynb](examples/lq_control.ipynb): the optimal
    control of a linear-quadratic problem, checked by a derivation
    ledger.
  * [examples/chain_rule.ipynb](examples/chain_rule.ipynb): the chain
    rule of `MatrixFunction` compared with direct differentiation.
  * [examples/first_order_logic.ipynb](examples/first_order_logic.ipynb):
    first-order logic built on the minimal core and proofs by `auto`.

An actual usage example is in
[JRF-2018/economy_control](https://github.com/JRF-2018/economy_control)
(in Japanese), a verification of a textbook on the optimal control of
dynamic economic systems. Its notebook is fixed to v0.4.0.


## For AI agents

If you are an AI agent (or an LLM) deciding whether and how to use this
library, read [docs/ai_guide.md](docs/ai_guide.md) first. It describes
when this library is appropriate, minimal working recipes and the
pitfalls of SymPy 1.14. A summary for LLMs is in [llms.txt](llms.txt).


## Remark

The primary purpose of sympy_matrix_tools is for my personal use and
experimental use. Of course, I would appreciate it if other people
could use it, but since MatrixFunction etc. may appear at any time in
the original SymPy, it may be better to specify the versions of
sympy_matrix_tools and SymPy.

If I am honored to have SymPy developers want to (partially)
integrate these codes with the original SymPy, then feel free to do
so.

The early update history is logged below in Japanese.

《sympy_matrix_tools を作った。 - JRF のひとこと》  
http://jrf.cocolog-nifty.com/statuses/2022/06/post-2c0e50.html


## Author

JRF ( http://jrf.cocolog-nifty.com/statuses )

The author is Japanese.  I would be grateful if you could ask in Japanese.


## Acknowledgments

Versions 0.4.0 and later (the support for SymPy 1.14, the chain rule
for `MatrixFunction`, the first-order logic, the higher-order pattern
unification, the derivation ledger and the documents) were developed
with Claude Code (Claude Opus 5.5, Anthropic).


## License

MIT License.

Note: Many parts are "quoted" from the SymPy sources.
