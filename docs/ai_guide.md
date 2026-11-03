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
  * You are asked to **prove or check an argument** that mixes logic
    (∀, ∃, →) and matrix algebra, e.g. properties of an order of
    matrices defined by the user: see "Proving or checking an
    argument" below.

## When not to use it

  * Numerical linear algebra: use NumPy or SciPy.
  * Rigorous proofs: use Lean or Isabelle.  The logic tools here are
    experimental, and SymPy's automatic evaluation is part of the
    trusted base.  A `Derivation` gives evidence, not proofs.  The
    logic tools are good for checking the structure of an argument of
    the size of a textbook page, with the algebra checked by a
    `Derivation` (see the workflow below).
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

Write the statement in the logical notation (or build it with
`ForAll`, `Exists`, `Implies`, ... and `PredApply`):

```python
>>> from sympy_matrix_tools.fol import Theorem, by_auto, exI, SAFE
>>> Theorem("ex_all_swap", "(∃x. ∀y. P(x, y)) → (∀y. ∃x. P(x, y))",
...         by_auto(exI, safe=SAFE)).verify()
True

```

`verify()` raises `ValueError` if the proof search fails.  See
[fol.md](fol.md).


## Recipe 5: prove a statement about matrices

Definitions are axioms of your own, written in the logical notation.
The algebra that the logic needs is derived in a ledger and used as a
theorem.  Give `symbols` for every name, including the dimension.

```python
>>> from sympy_matrix_tools.fol import Proof, by_rule, by_rewrite, by_steps
>>> S = {"n": n}
>>> S.update({k: MatrixSymbol(k, n, n) for k in ["A", "B", "M", "N"]})
>>> S.update({"x": MatrixSymbol("x", n, 1),
...           "a": MatrixSymbol("a", 1, 1), "b": MatrixSymbol("b", 1, 1)})
>>> M, N, x = S["M"], S["N"], S["x"]

>>> nonneg_add = Theorem("nonneg_add",
...     "∀a b. nonneg(a) ∧ nonneg(b) → nonneg(a + b)", axiom=True, symbols=S)
>>> psd_def = Theorem("psd_def", "∀M. psd(M) ↔ (∀x. nonneg(x.T*M*x))",
...                   axiom=True, symbols=S)

>>> d = Derivation("algebra")
>>> d.calc("quad_add", x.T*(M + N)*x, (x.T*M*x + x.T*N*x, cas(mat_expand)))
>>> quad_add = d.theorem("quad_add")

>>> psd_unfold = Theorem("psd_unfold", "∀M x. psd(M) → nonneg(x.T*M*x)",
...                      by_auto(psd_def), symbols=S)
>>> psd_add = Theorem("psd_add", "∀M N. psd(M) ∧ psd(N) → psd(M + N)",
...                   by_steps(by_rule(psd_def), by_rewrite(quad_add),
...                            by_auto(nonneg_add, psd_unfold)),
...                   symbols=S)
>>> psd_add.verify()
True
>>> print(psd_add.explain())
theorem psd_add: ∀M N. (psd(M) ∧ psd(N)) → psd(M + N)
  proof: rule(psd_def); rewrite(quad_add); auto(nonneg_add, psd_unfold)
  axioms: psd_def, nonneg_add
  ledger facts: quad_add
  theorems: psd_unfold
  verified: yes

```

A longer argument is easier as a forward `Proof`, where each `have`
is checked where it is written:

```python
>>> pf = Proof("psd_add3", "∀A B M. psd(A) ∧ psd(B) ∧ psd(M) → psd(A + B + M)",
...            symbols=S)
>>> pf.have("ab", "psd(A + B)", by_auto(psd_add), using=["h1", "h2"])
>>> th = pf.show(by_auto(psd_add), using=["ab", "h3"])
>>> print(pf)
proof psd_add3: ∀A B M. (psd(A) ∧ psd(B) ∧ psd(M)) → psd(A + B + M)
  fix A B M
  assume h1: psd(A)
  assume h2: psd(B)
  assume h3: psd(M)
  have ab: psd(A + B)
      using h1, h2  by auto(psd_add)
  show psd(A + B + M)
      using ab, h3  by auto(psd_add)
qed

```

See [fol.md](fol.md) and the long example
[examples/loewner_order.ipynb](../examples/loewner_order.ipynb) (the
Loewner order, with two ledgers and backward and forward proofs).


## Axioms are fine: what this framework is for

The logic of this library does not aim at proofs without gaps from
first principles, as Lean or Isabelle/HOL do.  It aims at **plausible
arguments whose every step is checked**: the structure of the argument
by the logic, and the algebra by SymPy and the derivation ledger.

So introducing axioms is normal, not a "sin": definitions, facts
about numbers (`nonneg_add`: a ≥ 0 ∧ b ≥ 0 → a + b ≥ 0), induction
(`nat_induct` below), or a known theorem of a textbook that you do not
want to prove again.  What matters is that

  * each axiom is plausible (a definition, a standard fact or a cited
    theorem) and is shown to the user; `explain()` lists every axiom
    that a theorem rests on;
  * the claim itself, or a step of it, is **never** an axiom.


## Recipe 6: induction

Induction is an axiom with a predicate variable `Phi`.  `by_induct(rule,
var)` instantiates `Phi` with the goal, as `induct` of Isabelle does,
and leaves the base case and the step as subgoals.  A sequence is a
function of a natural number (a `Function`, or a `MatrixFunction` for a
sequence of matrices), defined by axioms.  `by_rewrite_asm()` rewrites
by an equation among the premises of the subgoal, such as the
induction hypothesis.

```python
>>> from sympy_matrix_tools.fol import by_induct, by_rewrite_asm, by_cas
>>> k = Symbol("k", integer=True, nonnegative=True)
>>> s = Function("s")
>>> T = {"k": k, "s": s}
>>> nat_induct = Theorem("nat_induct",
...     "∀Phi k. Phi(0) ∧ (∀k. Phi(k) → Phi(k + 1)) → Phi(k)",
...     axiom=True, symbols=T)
>>> s_0 = Theorem("s_0", "s(0) = 0", axiom=True, symbols=T)
>>> s_step = Theorem("s_step", "∀k. s(k + 1) = s(k) + 2*k + 1",
...                  axiom=True, symbols=T)
>>> sq = Theorem("sq", "∀k. s(k) = k**2",
...     by_steps(by_induct(nat_induct, "k"),
...              by_rewrite(s_0),        # the base case: s(0) = 0
...              by_rewrite(s_step),     # the step: s(k + 1) = ...
...              by_rewrite_asm(),       # s(k) = k**2
...              by_cas("expand")),      # k**2 + 2*k + 1 = (k + 1)**2
...     symbols=T)
>>> sq.verify()
True
>>> print(sq.explain())
theorem sq: ∀k. s(k) = k**2
  proof: induct(nat_induct, k); rewrite(s_0); rewrite(s_step); rewrite_asm(0); cas(expand)
  axioms: nat_induct, s_0, s_step
  verified: yes

```

The Loewner example proves `V_mono` (Vk ≼ Vk+1 for all k) and
`V_nonneg` in the same way, for a sequence of matrices.


## Proving or checking an argument: a workflow

When you are asked to "prove (or check) this argument with
sympy_matrix_tools", work as follows.

  1. **Install and import.**
     `pip install "sympy>=1.14" git+https://github.com/JRF-2018/sympy_matrix_tools`,
     then `from sympy import *`, `from sympy_matrix_tools import *`,
     `from sympy_matrix_tools.fol import *` and
     `from sympy_matrix_tools.derivation import Derivation, cas, rule`.
  2. **Split the argument** into three kinds of parts:
       * definitions and assumptions: axioms, `Theorem(name, "...",
         axiom=True, symbols=S)`.  Axioms are fine here (see "Axioms
         are fine" above): definitions, facts about numbers,
         induction, cited theorems.  Keep them plausible and show
         them to the user.  **Never** state the claim (or a step of it)
         as an axiom.
       * algebra (identities of matrices and scalars): a `Derivation`,
         each step justified by `cas(...)` or `rule(...)`, then
         `check()` and `d.theorem(name)`.
       * logical steps: `Theorem` with tactics (`by_auto`, `by_rule`,
         `by_rewrite`, `by_cas`, `by_steps`) or a forward `Proof`.
  3. **Build up from small lemmas.**  For a long argument, follow the
     text of the argument with a `Proof`: one `have` per sentence.  A
     wrong step raises `ValueError` where it is written.
  4. **When a step fails**, try in this order:
       * check that `symbols` has every name and dimension (a missing
         name becomes a new scalar `Symbol`, and nothing matches);
       * give an instance by `th.instance(x=...)` when the unifier has
         to find a term inside a product (e.g. `P*x` in `x.T*P.T*M*P*x`);
       * rewrite by a ledger equation (`by_rewrite(eq)`, or
         `reverse=True`, or only one `occurrence`) when SymPy does not
         apply the algebra by itself, by an equation among the
         premises (`by_rewrite_asm()`), or compute by
         `by_cas("expand")`;
       * for a statement about all natural numbers, use induction
         (Recipe 6);
       * split the step into smaller `have`s, pass fewer facts in
         `using`, or raise `depth` and `max_nodes` of `by_auto`.
     If the step is still not proved, it may be false: check it
     numerically (e.g. with a `Derivation`) before trying harder.
  5. **Report honestly.**  Print each `Proof` and the `explain()` of the
     main theorem, and say what the result rests on: the axioms (which
     are assumed, not proved), the ledger facts (trusted to SymPy's
     computation and checked with random values) and SymPy's automatic
     evaluation.  If `verify()` raised, say that the claim was not
     proved; do not present a partial proof as a proof.

Limits to tell the user about:

  * Induction, arithmetic and inequalities of numbers are not built
    in; state what you need (such as `nat_induct` and `nonneg_add`)
    as axioms, and say so.
  * `auto` is incomplete and takes about 0.1 seconds per proofstate.
  * For a proof that must be rigorous, translate the result to Lean or
    Isabelle.


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
    Importing `sympy_matrix_tools.fol` applies the latter.
  * In formulas, SymPy reorders and simplifies terms: `B - A` is
    printed as `-A + B`, and `(P*x).T` is `x.T*P.T`.  Write formulas in
    any form; compare them after parsing, not as text.
  * A quadratic form `x.T*M*x` is a 1×1 matrix, not a scalar; declare
    the "scalars" of axioms about it as 1×1 `MatrixSymbol`s.

The full list of known SymPy bugs is in [fixes.md](fixes.md).
