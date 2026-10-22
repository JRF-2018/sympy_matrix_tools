# Guide for AI agents

[Back to README](../README.md)

This page is written for AI agents (and LLMs) that consider using
sympy_matrix_tools.  Every code block is a tested doctest, so it can
be copied as it is.


## When to use this library

  * You derive formulas with **symbolic matrices** in SymPy and need an
    **undefined matrix-valued function** `F(X)` or `x(t)` and its
    derivatives (the chain rule).  SymPy 1.14 has no such thing.
  * You want to **check a derivation step by step**, especially one
    written by you (an AI) or taken from a textbook: `Derivation`
    checks every step exactly with random rational matrices and tells
    which step is wrong.
  * You need coefficients or collection in **non-commutative** matrix
    expressions (`mat_coeff`, `mat_collect`, `mat_divide`), or sums of
    matrices (`MatSum`).

## When not to use it

  * Numerical linear algebra: use NumPy or SciPy.
  * Rigorous proofs: use Lean or Isabelle.  The logic tools here are
    experimental, and SymPy's automatic evaluation is part of the
    trusted base.  A `Derivation` gives evidence, not proofs.
  * Large expressions: the tools are written for derivations of the
    size of a textbook page.


## Setup

```python
>>> from sympy import *
>>> from sympy_matrix_tools import *
>>> from sympy_matrix_tools.derivation import Derivation, Grad, cas, rule

```

`sympy_matrix_tools.fol` (first-order logic) and
`sympy_matrix_tools.derivation` are imported separately.  Importing
`sympy_matrix_tools` fixes a bug of SymPy 1.14 (see
[fixes.md](fixes.md)).


## Recipe 1: differentiate an undefined matrix function

```python
>>> t = Symbol("t")
>>> X = MatrixSymbol("X", 2, 2)
>>> A = MatrixSymbol("A", 2, 2)
>>> D = MatrixSymbol("D", 2, 2)
>>> F = MatrixFunction("F", 2, 2)

>>> F(t * X).diff(t)          # by a scalar: a directional derivative
Subs(Derivative(F(t*X + _s_1*X), _s_1), _s_1, 0)
>>> F(t * X).diff(t).subs(F, Lambda(D, D * A * D)).doit()
(2*t)*X*A*X

>>> r = mat_diff(F(X), X)     # by a matrix: an array expression
>>> r
Subs(Derivative(F(_xi_1), _xi_1), _xi_1, X)
>>> r.shape                   # the indices of X come first
(2, 2, 2, 2)

```

Use `mat_diff(expr, X)` instead of `expr.diff(X)` for a matrix `X`:
`.diff` uses SymPy's `convert_array_to_matrix`, which returns wrong
results for some arrays of rank 4.  To check a derivative, substitute
a concrete `Lambda` for the function and call `doit()`.


## Recipe 2: check a derivation

Write each intermediate expression and its justification.  Use
`cas(...)` for a computation by SymPy, `rule(name)` for a rewrite by a
named fact, and nothing for a step you only assert.

```python
>>> n = Symbol("n", integer=True, positive=True)
>>> W = MatrixSymbol("W", n, n)
>>> x = MatrixSymbol("x", n, 1)
>>> b = MatrixSymbol("b", n, 1)
>>> d = Derivation("gradient of a quadratic form")
>>> d.hypothesis(Eq(W.T, W), "W_sym")
>>> d.calc("grad", Grad(x.T * W * x / 2 + b.T * x, x),
...        (W * x / 2 + W.T * x / 2 + b, cas("doit")),
...        (W * x + b, rule("W_sym")))
>>> d.check(trials=3)
True
>>> d.summary()
{'steps': 2, 'asserted': 0, 'symbolic_failures': 0, 'numeric_failures': 0, 'open_obligations': 0}

```

How to read the result:

  * `numeric: FAILED` means the step is **wrong** (a counterexample is
    in `step.counterexample`).
  * `symbolic: FAILED` with `numeric: ok` usually means that the
    symbolic check is too weak; the step is probably right.
  * `numeric: n/a` means that the step uses a hypothesis that random
    values do not satisfy (e.g. a first-order condition).
  * `asserted` steps have no justification; prefer to give one.
  * An `OPEN` obligation (e.g. `Q.invertible(K)`) must be assumed or
    proved.

See [derivation.md](derivation.md) for `solve_linear`, obligations
and the LaTeX output.


## Recipe 3: coefficients in non-commutative expressions

```python
>>> M1 = MatrixSymbol("M1", n, n)
>>> M2 = MatrixSymbol("M2", n, n)
>>> z = M1 * M2 + 3 * M2 * M1 + M1 * M1
>>> z
M1**2 + M1*M2 + 3*M2*M1
>>> mat_coeff(z, M1)                # the first term of the form C*M1
3*M2
>>> mat_coeff(z, M1, nth=1, expand_pow=True)   # M1**2 as M1*M1
M1
>>> mat_coeff(z, M1, right=True)    # the first term of the form M1*C
M2
>>> mat_collect(z, M1, expand_pow=True)        # collect all C*M1
(M1 + 3*M2)*M1 + M1*M2
>>> mat_collect(z, M1, right=True, expand_pow=True)
M1*(M1 + M2) + 3*M2*M1

```

`mat_coeff` returns the coefficient of **one** matching term (`nth`
selects it), not the sum; use `mat_collect` to collect all terms.
SymPy turns `M1*M1` into `M1**2`, so pass `expand_pow=True` to treat
powers as products.


## Recipe 4: prove a small first-order statement

```python
>>> from sympy_matrix_tools.fol import Theorem, Exists, by_auto, exI, SAFE
>>> x, y = symbols("x y")
>>> P = Symbol("P")
>>> g = Implies(Exists(Lambda(x, ForAll(Lambda(y, PredApply(P, x, y))))),
...             ForAll(Lambda(y, Exists(Lambda(x, PredApply(P, x, y))))))
>>> Theorem("ex_all_swap", g, by_auto(exI, safe=SAFE)).verify()
True

```

`verify()` raises `ValueError` if the proof search fails.  See
[fol.md](fol.md).


## Pitfalls (SymPy 1.14)

  * `Eq(Grad(...), ZeroMatrix(...))` evaluates the gradient; pass
    `evaluate=False`.
  * `Grad(f, x)` needs `f` written explicitly in `x`; a symbol defined
    by `Derivation.define` is independent of `x` inside `Grad`.
  * `Q` is SymPy's assumption namespace (`Q.invertible(K)`); do not
    shadow it by a matrix named `Q`.
  * `Predicate` takes one argument; use `PredApply(P, x, y)`.
  * Write `.diff(t).diff(t)` instead of `.diff(t, 2)`, and `X*X`
    instead of `X**2`, for matrix expressions.
  * `MatMul.expand()` leaves scalar coefficients; use `mat_expand`.
  * Write the derivative of a `MatrixFunction` as `x(t).diff(t)`, not
    `Derivative(x(t), t)` (which is a scalar `Derivative`).
  * `derivative_str(expr)` shows derivatives as `x'(t)`, `F'(a)` and
    `DF(a)[d]`, which is easier to read than the default `Subs(...)`.
  * Dummy names such as `_xi_1` and `_s_1` in derivatives are bound
    variables of `Subs`; they are not free symbols.
  * `fix_ExpectationMatrix_expand()` and
    `fix_Basic_MatrixSymbol_about_dummy()` change SymPy globally.

The full list of known SymPy bugs is in [fixes.md](fixes.md).
