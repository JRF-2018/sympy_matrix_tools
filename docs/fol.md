# First-order logic on the minimal core

[Back to README](../README.md)

**Experimental.**

The proofstates of [logic.md](logic.md) have only `ForAll`, `Implies`
and `And` as connectives.  `sympy_matrix_tools.fol` builds the
classical first-order logic on this minimal core, as in higher-order
logic: the other connectives are **defined** from `ForAll` and
`Implies`, and their rules are **proved** as theorems from the
definitions, not assumed.  This shows that the core is enough for
first-order logic, and tests the soundness of the core and of the
unification.

It is not meant to replace Lean or Isabelle; the trusted base includes
SymPy's automatic evaluation.

```python
>>> from sympy import *
>>> from sympy_matrix_tools import *
>>> from sympy_matrix_tools.fol import *

```


## Writing formulas

| Formula | Written as |
|---|---|
| ∀x. φ | `ForAll(Lambda(x, φ))` |
| ∃x. φ | `Exists(Lambda(x, φ))` |
| φ → ψ, φ ∧ ψ, φ ∨ ψ | `Implies(φ, ψ)`, `And(φ, ψ)`, `Or(φ, ψ)` |
| ⊥, ¬φ | `Bot`, `neg(φ)` (which is `Implies(φ, Bot)`) |
| P(x), R(x, y) | `PredApply(P, x)`, `PredApply(R, x, y)` with symbols `P`, `R` |

A predicate is a plain symbol applied by `PredApply`, because SymPy's
`Predicate` takes only one argument and cannot be a variable.  A
symbol that stands for a formula needs `boolean=True`.
`formula_str(φ)` shows a formula in the logical notation.

SymPy's `false` and `Not` are not used for ⊥ and ¬, because SymPy
evaluates `Implies(A, false)` to `~A` at once.  SymPy's `Or` is used;
its automatic simplification is classically sound.


## Writing formulas as text

`parse_formula(text)` reads a formula written in the logical notation,
the same notation that `formula_str` and `print_proofstate(...,
pretty=True)` print, so that their output can be read back.  ASCII
alternatives are accepted: `ALL`, `EX`, `->`, `<->`, `|`, `&`, `~` and
`Bot`.

  * A binder `∀x y.` or `∃x.` extends as far to the right as possible.
    The precedence is `↔` < `→` (right associative) < `∨` < `∧` < `¬`.
  * `P(t, ...)` is `PredApply(P, t, ...)`; the arguments are SymPy
    expressions such as `f(x)` or `x + 1`.  A bare name `A` is a
    proposition `Symbol("A", boolean=True)`.
  * `?x` is a metavariable (`Wild`), `?P(x)` is
    `PredApply(Wild("P"), x)`, and `?f(x)` in a term is
    `Apply(Wild("f"), x)`.
  * `[expr]` is any SymPy expression, e.g. `[x < 3]`, and `*φ` is
    `FrozenGoal(φ)`, the goal of a proofstate.
  * `symbols={"t": t}` gives the objects for free names, and
    `predicates={"b": Predicate("b")}` applies SymPy predicates instead
    of `PredApply`.

`Theorem` accepts such a string as its statement.

```python
>>> f = parse_formula("ALL x. P(x) -> (EX y. R(f(x), y))")
>>> f
Q.ForAll(Lambda(x, Implies(PredApply(P, x), Q.Exists(Lambda(y, PredApply(R, f(x), y))))))
>>> print(formula_str(f))
∀x. P(x) → (∃y. R(f(x), y))
>>> parse_formula(formula_str(f)) == f
True
>>> A, B = symbols("A B", boolean=True)
>>> parse_formula("¬(A ∨ B) → ¬A ∧ ¬B") == Implies(neg(Or(A, B)), And(neg(A), neg(B)))
True
>>> Theorem("syllogism",
...         "(∀x. P(x) → Q(x)) ∧ (∀x. Q(x) → R(x)) → (∀x. P(x) → R(x))",
...         by_auto()).verify()
True

```


## The theory

The definitions are given as pairs of axioms (fold and unfold):

  * `Exists(Lambda(x, P(x)))` == ∀R. (∀x. P(x) → R) → R
    (`ex_fold`, `ex_unfold`)
  * `A | B` == ∀R. (A → R) ∧ (B → R) → R (`or_fold`, `or_unfold`)
  * `Bot` == ∀R. R (`bot_unfold`)

The classical logic is given by the axiom `ccontr`: ¬¬A → A.  These are
listed in `AXIOMS`.  The usual rules are theorems proved from them
(`THEOREMS`): `exI`, `exE`, `disjI1`, `disjI2`, `disjE`, `botE` and
`excluded_middle`.  `verify()` runs the proof of a theorem.

```python
>>> for th in [ex_fold, ex_unfold, ccontr, exI, exE, excluded_middle]:
...     print(th)
axiom ex_fold: ∀P. (∀R. (∀x. P(x) → R) → R) → (∃y. P(y))
axiom ex_unfold: ∀P R. ((∃y. P(y)) ∧ (∀x. P(x) → R)) → R
axiom ccontr: ∀A. ¬¬A → A
theorem exI: ∀P x. P(x) → (∃y. P(y))
theorem exE: ∀C P. ((∃y. P(y)) ∧ (∀x. P(x) → C)) → C
theorem excluded_middle: ∀A. A ∨ ¬A
>>> all(th.verify() for th in THEOREMS)
True

```


## Stating and proving a theorem

`Theorem(name, statement, proof)` is a named statement with its proof,
where `proof` is a function from the initial proofstate to a
proofstate without subgoals.  `verify()` runs the proof, checks that
it proves the statement, and raises `ValueError` otherwise.  A proof
must use only the theorems proved before it.

`by_auto(*rules, depth=6, max_nodes=2000, safe=())` makes a proof by
the tactic `auto` (see [logic.md](logic.md)):

  * `rules` are the theorems that `auto` may apply backward, e.g.
    `exI` to prove an `Exists`.
  * `safe` is a list of `(rule, kind)` of elimination rules that are
    applied to assumptions without backtracking.  `SAFE` contains the
    eliminations of `Exists` and `Or` (`exE` and `disjE`); pass
    `safe=SAFE` in most proofs.
  * `depth` limits the length of a proof, and `max_nodes` the number
    of proofstates explored.  A proof search that fails may take time
    until it reaches the limits.

```python
>>> x, y = symbols("x y")
>>> P = Symbol("P")
>>> g = Implies(Exists(Lambda(x, ForAll(Lambda(y, PredApply(P, x, y))))),
...             ForAll(Lambda(y, Exists(Lambda(x, PredApply(P, x, y))))))
>>> print(formula_str(g))
(∃x. ∀y. P(x, y)) → (∀y. ∃x. P(x, y))
>>> th = Theorem("ex_all_swap", g, by_auto(exI, safe=SAFE))
>>> th.verify()
True

```

The proof found by `auto` can be shown step by step.  Each `prem` is a
subgoal, `gl*` is the goal, and `?x` is a metavariable, which depends
on the bound variables in its arguments.

```python
>>> path = auto(make_proofstate(th.statement), [exI.statement],
...             safe=[(r.statement, i) for r, i in SAFE], trace=True)
>>> for z in path:
...     print("---")
...     print_proofstate(z, pretty=True)
---
prem 0: ∀y. (∃x. ∀y2. P(x, y2)) → (∃x1. P(x1, y))
gl*: ∀y1. (∃x2. ∀y3. P(x2, y3)) → (∃x3. P(x3, y1))
---
prem 0: ∀x y. (∀y2. P(x, y2)) → (∃x1. P(x1, y))
gl*: ∀y1. (∃x2. ∀y3. P(x2, y3)) → (∃x3. P(x3, y1))
---
prem 0: ∀x1 y. (∀y2. P(x1, y2)) → P(?x(x1, y), y)
gl*: ∀y1. (∃x2. ∀y3. P(x2, y3)) → (∃x3. P(x3, y1))
---
gl*: ∀y. (∃x. ∀y1. P(x, y1)) → (∃x1. P(x1, y))

```


## A classical theorem

The drinker paradox ∃x. (D(x) → ∀y. D(y)) needs the classical logic.
The witness of the outer ∃ remains a metavariable, which means that
the domain is not empty, as usual in first-order logic.

```python
>>> D = Symbol("D")
>>> drinker = Theorem("drinker",
...     Exists(Lambda(x, Implies(PredApply(D, x),
...                              ForAll(Lambda(y, PredApply(D, y)))))),
...     by_auto(excluded_middle, exI, botE, ccontr,
...             depth=12, max_nodes=4000, safe=SAFE))
>>> drinker.verify()
True

```

`tests/test_fol.py` proves more theorems (the syllogism, distributions
of ∃ over ∨ and of ∀ over ∧, ¬∃ ↔ ∀¬, De Morgan, Peirce's law, ¬∀ →
∃¬) and checks that non-theorems such as ∀y∃x R(x, y) → ∃x∀y R(x, y)
are rejected.  See also
[examples/first_order_logic.ipynb](../examples/first_order_logic.ipynb).


## Matrices, equations and algebra

The same logic works with matrices.  A theorem about matrices usually
needs some algebra, which the logic does not know; the algebra comes
from SymPy (its automatic evaluation and `by_cas`) and from a
derivation ledger (`Derivation.theorem`, see
[derivation.md](derivation.md)).  A long example, the Loewner order
of matrices, is in
[examples/loewner_order.ipynb](../examples/loewner_order.ipynb).

| Name | Purpose |
|---|---|
| `by_auto(*rules, ...)` | search a proof backward with the rules (see above) |
| `by_rule(th)` | one backward step with `th` (`rule` of Isabelle); it may leave subgoals |
| `by_rewrite(eq, reverse=False)` | rewrite the conclusion of the first subgoal by an equation theorem |
| `by_cas(func="doit")` | apply a computation of SymPy to the terms of the first subgoal |
| `by_steps(*tactics)` | apply tactics in order |
| `th.instance(x=t, ...)` | the instance of `th` at `x := t` |
| `th.cas(func="doit")` | `th` with its terms computed by SymPy |
| `d.theorem(name)` | the derived fact `name` of a `Derivation` `d` as a theorem |
| `Proof(name, statement)` | a structured proof in forward style: `have(...)`, then `show(...)` |
| `th.dependencies()`, `th.explain()` | what `th` rests on |

### Formulas about matrices

`parse_formula` (and `Theorem`) take `symbols`, a dict from names to
objects.  A bound variable whose name is in `symbols` is that object,
so `∀M. ...` binds a `MatrixSymbol` `M`.  Put the dimensions (here
`n`) in `symbols` too, for terms such as `ZeroMatrix(n, n)`.  `s = t`
is the equation `Eq(s, t)`.

Importing `sympy_matrix_tools.fol` applies
`fix_Basic_MatrixSymbol_about_dummy()` (see [fixes.md](fixes.md)),
because SymPy's canonical forms break formulas that bind matrices.

```python
>>> n = Symbol("n", integer=True, positive=True)
>>> S = {"n": n}
>>> S.update({k: MatrixSymbol(k, n, n) for k in ["A", "B", "M", "N", "P"]})
>>> S.update({k: MatrixSymbol(k, n, 1) for k in ["x"]})
>>> S.update({k: MatrixSymbol(k, 1, 1) for k in ["a", "b"]})
>>> A, B, M, N, P, x = [S[k] for k in ["A", "B", "M", "N", "P", "x"]]

>>> nonneg_add = Theorem("nonneg_add",
...                      "∀a b. nonneg(a) ∧ nonneg(b) → nonneg(a + b)",
...                      axiom=True, symbols=S)
>>> psd_def = Theorem("psd_def", "∀M. psd(M) ↔ (∀x. nonneg(x.T*M*x))",
...                   axiom=True, symbols=S)
>>> psd_def.statement.arguments[0].signature == (M,)
True

```

Here `nonneg(a)` means a ≥ 0 for a 1×1 matrix `a` (a quadratic form
`x.T*M*x` is a 1×1 matrix in SymPy), and `psd(M)` means that `M` is
positive semidefinite.  These are axioms of our own, not part of the
library.  A theorem whose conclusion is `↔` (or `∧`) is used by the
tactics as its parts, i.e. in both directions.

After each step, `MatMul`, `MatAdd` and `Transpose` are evaluated as
their operators do: `(P*x).T*M*(P*x)` is `x.T*P.T*M*P*x`, `x.T*0*x` is
`0` and `(B - A) + A` is `B`.  This automatic evaluation of SymPy is
part of the trusted base.

### Tactics

A tactic is a function from a proofstate to a proofstate (or `None` if
it fails); the proof of a `Theorem` is a tactic.  A tactic records the
theorems that it uses.

  * `by_rule(th, index=0)` resolves the subgoal of the index with `th`
    backward, like `rule` of Isabelle.
  * `by_rewrite(eq, index=0, reverse=False, num=0)` rewrites the
    conclusion of the subgoal by `eq`, a theorem `∀ws. C → l = r`: the
    first subterm that is an instance of `l` is replaced by the
    instance of `r` (with `reverse=True`, `r` by `l`), and the
    conditions `C` become subgoals.  The step is an instance of the
    substitution of equals, `s = t ∧ P(t) → P(s)`.
    (`logic_tools.rewrite_goal`.)
  * `by_cas(func="doit", index=0)` applies a computation of SymPy (a
    function or a method name such as `"expand"`) to the terms in the
    conclusion of the subgoal.  (`logic_tools.cas_goal`.)
  * `by_steps(*tactics)` applies the tactics in order; the last one
    must solve all the subgoals.

The steps can be run by hand to see the proofstates.  Here the
algebra x.T*(M + N)*x = x.T*M*x + x.T*N*x is assumed as an axiom; the
next section derives it in a ledger instead.

```python
>>> quad_add_ax = Theorem("quad_add_ax",
...     "∀M N x. x.T*(M + N)*x = x.T*M*x + x.T*N*x", axiom=True, symbols=S)
>>> psd_unfold = Theorem("psd_unfold", "∀M x. psd(M) → nonneg(x.T*M*x)",
...                      by_auto(psd_def), symbols=S)
>>> z = make_proofstate(parse_formula("∀M N. psd(M) ∧ psd(N) → psd(M + N)",
...                                   symbols=S))
>>> z = by_rule(psd_def)(z)
>>> print_proofstate(z, pretty=True)
prem 0: ∀M N x. (psd(M) ∧ psd(N)) → nonneg(x.T*(M + N)*x)
gl*: ∀M1 N1. (psd(M1) ∧ psd(N1)) → psd(M1 + N1)
>>> z = by_rewrite(quad_add_ax)(z)
>>> print_proofstate(z, pretty=True)
prem 0: ∀M N x. (psd(M) ∧ psd(N)) → nonneg(x.T*M*x + x.T*N*x)
gl*: ∀M1 N1. (psd(M1) ∧ psd(N1)) → psd(M1 + N1)
>>> z = by_auto(nonneg_add, psd_unfold)(z)
>>> print_proofstate(z, pretty=True)
gl*: ∀M N. (psd(M) ∧ psd(N)) → psd(M + N)

```

### Instances

`th.instance(x=t, ...)` replaces the outer bound variables named by
the keys (forall elimination) and binds the new free symbols of the
terms again (forall introduction), unless `over` gives the symbols to
bind.  Use it when unification cannot find a term, e.g. the vector
`P*x` in `x.T*P.T*M*P*x`.  `th.cas(func)` computes the terms of a
theorem by SymPy, like `simplified` of Isabelle.

```python
>>> psd_unfold.instance(x=P*x)
instance psd_unfold[x := P*x]: ∀M P x. psd(M) → nonneg(x.T*P.T*M*P*x)
>>> psd_congr = Theorem("psd_congr", "∀M P. psd(M) → psd(P.T*M*P)",
...                     by_auto(psd_def, psd_unfold.instance(x=P*x)),
...                     symbols=S)
>>> psd_congr.verify()
True

```

### Theorems from a derivation ledger

`d.theorem(name)` turns the derived fact `name` of a `Derivation` `d`
into a theorem of kind `"ledger"`, universally quantified over its
symbols; the hypotheses and definitions that the fact rests on become
premises.  Its `verify()` requires that every step is justified and
checked symbolically, and that no numerical check fails.  It is
trusted to SymPy's computation, with the numerical checks as further
evidence; it is not a proof in the logic.

```python
>>> from sympy_matrix_tools.derivation import Derivation, cas
>>> d = Derivation("Quadratic forms")
>>> d.calc("quad_add", x.T*(M + N)*x, (x.T*M*x + x.T*N*x, cas(mat_expand)))
>>> quad_add = d.theorem("quad_add")
>>> quad_add.verify()
True
>>> psd_add = Theorem("psd_add", "∀M N. psd(M) ∧ psd(N) → psd(M + N)",
...                   by_steps(by_rule(psd_def), by_rewrite(quad_add),
...                            by_auto(nonneg_add, psd_unfold)),
...                   symbols=S)
>>> psd_add.verify()
True

```

### Structured proofs

`Proof(name, statement, symbols=None, assume=None)` proves
`∀xs. H1 ∧ ... → C` forward, like Isar of Isabelle: `xs` are fixed,
the assumptions are named `h1`, `h2`, ... (or by `assume`), and

  * `have(name, formula, by, using=None)` derives `formula` from the
    facts named in `using` (all the facts if `None`) by the tactic
    `by`.  The step is a lemma `∀xs. (used facts) → formula`, proved at
    once, so a wrong step raises `ValueError` where it is written.
  * `show(by, using=None)` derives `C` and returns the proved
    `Theorem`, whose proof combines the lemmas.
  * `print(proof)` shows the proof as text.

In formulas of a `Proof`, the names of `xs` are the fixed variables.

```python
>>> le_def = Theorem("le_def", "∀A B. mat_le(A, B) ↔ psd(B - A)",
...                  axiom=True, symbols=S)
>>> pf = Proof("le_congr", "∀A B P. mat_le(A, B) → mat_le(P.T*A*P, P.T*B*P)",
...            symbols=S)
>>> pf.have("d", "psd(B - A)", by_auto(le_def), using=["h1"])
>>> pf.have("e", "psd(P.T*(B - A)*P)", by_auto(psd_congr), using=["d"])
>>> d.calc("congr_diff", P.T*(B - A)*P, (P.T*B*P - P.T*A*P, cas(mat_expand)))
>>> pf.have("f", "psd(P.T*B*P - P.T*A*P)",
...         by_steps(by_rewrite(d.theorem("congr_diff"), reverse=True),
...                  by_auto()), using=["e"])
>>> le_congr = pf.show(by_auto(le_def), using=["f"])
>>> print(pf)
proof le_congr: ∀A B P. mat_le(A, B) → mat_le(P.T*A*P, P.T*B*P)
  fix A B P
  assume h1: mat_le(A, B)
  have d: psd(-A + B)
      using h1  by auto(le_def)
  have e: psd(P.T*(-A + B)*P)
      using d  by auto(psd_congr)
  have f: psd(-P.T*A*P + P.T*B*P)
      using e  by rewrite_rev(congr_diff); auto()
  show mat_le(P.T*A*P, P.T*B*P)
      using f  by auto(le_def)
qed

```

### What a theorem rests on

`th.dependencies()` lists the theorems that `th` rests on, and
`th.explain()` shows the statement, the proof and the axioms, ledger
facts and theorems used.  When you report a result, report these
too: a theorem is only as good as its axioms.

```python
>>> print(le_congr.explain())
theorem le_congr: ∀A B P. mat_le(A, B) → mat_le(P.T*A*P, P.T*B*P)
  proof: structured proof (4 steps)
  axioms: le_def, psd_def
  ledger facts: congr_diff
  theorems: psd_unfold, psd_congr
  verified: yes

```


## Remarks

  * The trusted base includes SymPy's automatic evaluation, e.g.
    `Implies(A, A)` is `True` at once, and the computations of
    `by_cas`, `Theorem.cas` and ledger theorems.
  * `auto` is not complete; a long proof may need larger `depth` and
    `max_nodes`.  A proof step takes about 0.1 seconds.
  * In 0.5.0, an unsound behavior of `remove_trivial_assumptions` was
    fixed: bound variables of a premise were dropped and could be
    unified with outer metavariables.
