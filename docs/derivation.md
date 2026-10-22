# Derivation ledger

[Back to README](../README.md)

```python
>>> from sympy import *
>>> from sympy_matrix_tools import *

```


**Experimental.**

When a formula is derived by hand, by an AI or with a CAS, a wrong
step is easy to make and hard to find.  `sympy_matrix_tools.derivation`
provides `Derivation`, a ledger of a derivation.  It is not a theorem
prover: it records why each step is believed, and checks every step
with random values, so that a wrong step is found at once.

  * Each step is justified by `cas(...)` (a computation by SymPy),
    `rule(name)` (a rewrite by a named hypothesis, definition or
    derived fact), or nothing (only asserted).
  * As in the calc of Lean, the intermediate expressions are written
    in the script, so the script can be read without running it.
  * `check()` replaces the symbols by random rational numbers and
    matrices that satisfy the hypotheses (symmetric, positive
    definite, invertible and the definitions), and compares both sides
    of each step exactly. It finds wrong steps whoever wrote them: a
    human, an AI or a buggy CAS.
  * `solve_linear` records the invertibility it needs as an
    obligation, which is discharged by a hypothesis or left open.

| Name | Purpose |
|---|---|
| `Derivation(title)` | a new ledger |
| `d.hypothesis(prop, name)` | add a hypothesis: an `Eq`, or a predicate such as `Q.invertible(K)` |
| `d.define(symbol, expr, name)` | define `symbol` as an abbreviation of `expr` |
| `d.calc(name, start, *steps)` | derive `start = last` by a chain of equalities; each step is `(expr, justification, ...)` or a bare `expr` |
| `cas(func, *args)` | a justification: `func(previous, *args)`, or a method name such as `"doit"` or `"expand"` |
| `rule(name)` | a justification: rewrite by the named fact, in either direction |
| `d.solve_linear(name, fact_name, x)` | derive `x = -K**-1*c` from the fact `K*x + c = 0`, with the obligation that `K` is invertible |
| `d.alias(expr, name)` | show `expr` as `name` in the reports (only for display) |
| `d.check(trials=3, seed=0)` | check every step numerically; return `True` if no failure is found |
| `d.report()`, `d.summary()` | a text report, and the counts of steps, assertions, failures and open obligations |
| `Grad(f, x)` | the unevaluated gradient of a scalar `f` by a matrix `x` (`doit()` computes it) |

Hypotheses of symmetry (`Eq(X.T, X)` or `Q.symmetric(X)`), positive
definiteness (`Q.positive_definite(X)`) and invertibility
(`Q.invertible(expr)`) are respected when random values are generated,
and definitions are computed from the other values.  Other hypotheses
(such as a first-order condition) are constraints that random values
do not satisfy, so the steps that use them are not checked
numerically.

The following derives the optimal control of a one-period
linear-quadratic problem (see
[examples/lq_control.ipynb](../examples/lq_control.ipynb)).

```python
>>> from sympy import Q    # The assumption predicates of SymPy.
>>> from sympy_matrix_tools.derivation import Derivation, Grad, cas, rule
>>> n, m = symbols("n m", integer=True, positive=True)
>>> A = MatrixSymbol("A", n, n)
>>> B = MatrixSymbol("B", n, m)
>>> W = MatrixSymbol("W", n, n)
>>> y0 = MatrixSymbol("y0", n, 1)
>>> z = MatrixSymbol("z", n, 1)
>>> p = MatrixSymbol("p", n, 1)
>>> x = MatrixSymbol("x", m, 1)
>>> K = MatrixSymbol("K", m, m)

>>> d = Derivation("One-period LQ control")
>>> d.hypothesis(Eq(W.T, W), "W_sym")
>>> d.hypothesis(Q.invertible(K), "K_inv")
>>> d.define(K, B.T*W*B, "K_def")
>>> y = A*y0 + B*x + z
>>> J = y.T*W*y/2 + p.T*y
>>> d.alias(J, "J")    # only for display
>>> d.calc("grad", Grad(J, x),
...        (B.T*W*y/2 + B.T*W.T*y/2 + B.T*p, cas("doit")),
...        (B.T*W*y + B.T*p, rule("W_sym")),
...        (K*x + B.T*(W*(A*y0 + z) + p), rule("K_def")))
>>> d.hypothesis(Eq(Grad(J, x), ZeroMatrix(m, 1), evaluate=False), "foc")
>>> d.calc("foc_linear", K*x + B.T*(W*(A*y0 + z) + p),
...        (Grad(J, x), rule("grad")),
...        (ZeroMatrix(m, 1), rule("foc")))
>>> d.solve_linear("x_opt", "foc_linear", x)
>>> d.check(trials=3)
True
>>> print(d.report())
Derivation: One-period LQ control
[hypothesis] W_sym: W.T = W
[hypothesis] K_inv: Q.invertible(K)
[definition] K_def: K = B.T*W*B
[derived] grad: Grad(J, x) = K*x + B.T*(W*(A*y0 + z) + p)
  Grad(J, x)
  = (1/2)*B.T*W*(A*y0 + B*x + z) + B.T*p + (1/2)*B.T*W.T*(A*y0 + B*x + z)    by cas doit    symbolic: ok    numeric: ok 3/3
  = B.T*W*(A*y0 + B*x + z) + B.T*p    by rule W_sym    symbolic: ok    numeric: ok 3/3
  = K*x + B.T*(W*(A*y0 + z) + p)    by rule K_def    symbolic: ok    numeric: ok 3/3
[hypothesis] foc: Grad(J, x) = 0
[derived] foc_linear: K*x + B.T*(W*(A*y0 + z) + p) = 0
  K*x + B.T*(W*(A*y0 + z) + p)
  = Grad(J, x)    by rule grad    symbolic: ok    numeric: ok 3/3
  = 0    by rule foc    symbolic: ok    numeric: n/a
[derived] x_opt: x = -K**(-1)*B.T*(W*(A*y0 + z) + p)
  by solve_linear(foc_linear, x)    numeric: ok 3/3
  obligation Q.invertible(K): discharged by K_inv
Summary: 5 steps, 0 asserted, 0 symbolic failures, 0 numeric failures, 0 open obligations.

```

The step by the first-order condition `foc` is not checked
numerically ("n/a"), because random values do not satisfy it. The
solution of `solve_linear` is checked by substituting it into the
equation.

A wrong step is found even if it claims to be a computation by SymPy:

```python
>>> d2 = Derivation("A wrong step")
>>> d2.calc("grad", Grad(J, x),
...         (B.T*W.T*y + B.T*p/2, cas("doit")))
>>> d2.check(trials=3)
False
>>> print(d2.report().splitlines()[3])
  = (1/2)*B.T*p + B.T*W.T*(A*y0 + B*x + z)    by cas doit    symbolic: FAILED    numeric: FAILED 0/3

```

In a notebook, a `Derivation` is shown as a LaTeX table by
`_repr_latex_`.


## Undefined functions in a ledger

Expressions may contain undefined functions (`MatrixFunction` and
SymPy's `Function`).  For the numerical checks, each function is
replaced by a random polynomial of its arguments (the same one in all
steps of a trial), and the derivatives are computed after the
replacement.  In reports, derivatives are shown by `derivative_str`:
`x'(t)`, `F'(a)`, `D_2F(a, b)` (by the second argument) and `DF(a)[d]`
(the directional derivative at `a` in the direction `d`).

```python
>>> t = Symbol("t")
>>> V = MatrixSymbol("V", n, n)
>>> xf = MatrixFunction("x", n, 1)
>>> d3 = Derivation("time derivative")
>>> d3.hypothesis(Eq(V.T, V), "V_sym")
>>> d3.calc("dq", (xf(t).T*V*xf(t)).diff(t),
...         (2*xf(t).diff(t).T*V*xf(t), rule("V_sym")))
>>> d3.check(trials=3)
True
>>> print(d3.report())
Derivation: time derivative
[hypothesis] V_sym: V.T = V
[derived] dq: x(t).T*V*x'(t) + x'(t).T*V*x(t) = 2*x'(t).T*V*x(t)
  x(t).T*V*x'(t) + x'(t).T*V*x(t)
  = 2*x'(t).T*V*x(t)    by rule V_sym    symbolic: ok    numeric: ok 3/3
Summary: 1 step, 0 asserted, 0 symbolic failures, 0 numeric failures, 0 open obligations.

```

Write the derivative of `x(t)` as `x(t).diff(t)`.  `Derivative(x(t),
t)` written directly is SymPy's scalar `Derivative`, not a matrix.


## The symbolic check

The symbolic check applies the justifications to the previous
expression and compares the result with the written expression after
`mat_expand`.  A rule rewrites in both directions, also the transposed
form of the fact, a part of a sum and a contiguous part of a product.
For (1, 1) matrices, a term and its transpose are identified, using the
symmetry `X.T = X` of the facts cited in the step.

Remarks:

  * `Grad(f, x)` needs `f` written explicitly in `x`. A symbol defined
    by `define` is treated as independent of `x` in `Grad`.
  * Hypotheses other than the symmetry, positive definiteness,
    invertibility and definitions are constraints that random values
    do not satisfy, and steps using them are not checked numerically.
  * The symbolic check is incomplete; a step that is correct may show
    "symbolic: FAILED" while "numeric: ok".
  * A `MatrixFunction` must be applied to arguments of the same kinds
    (scalar or matrix) everywhere in a ledger.
