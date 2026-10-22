# First-order logic on the minimal core

[Back to README](../README.md)

```python
>>> from sympy import *
>>> from sympy_matrix_tools import *

```


**Experimental.**

From 0.5.0, `sympy_matrix_tools.fol` builds the first-order logic on
the minimal core of the predicate logic above (ForAll, Implies and
And). The other connectives are defined impredicatively as in
higher-order logic, and their rules are proved as theorems:

  * `Exists(Lambda(x, P(x)))` == ∀R. (∀x. P(x) → R) → R
  * `A | B` == ∀R. (A → R) ∧ (B → R) → R
  * `Bot` == ∀R. R, and `neg(A)` == A → Bot
  * The classical logic is given by the axiom `ccontr`: ¬¬A → A.

```python
>>> from sympy_matrix_tools.fol import *
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

Predicates are symbols applied by `PredApply`, so that `P(x, y)` is
written as `PredApply(P, x, y)`. A theorem is proved by the tactic
`auto`, which searches proofs by backtracking with iterative
deepening. It uses the higher-order pattern unification `ho_unify`
and applies the eliminations of ∃ and ∨ in assumptions (`SAFE`)
without backtracking. `formula_str` and `print_proofstate(...,
pretty=True)` show formulas in the logical notation.

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

The classical drinker paradox is also provable. The witness of the
outer ∃ remains a metavariable, which means that the domain is not
empty, as usual.

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

See `tests/test_fol.py` for more examples and non-theorems that are
rejected.

Remarks:

  * SymPy's `false` and `Not` are not used for ⊥ and ¬, because SymPy
    evaluates `Implies(A, false)` to `~A` at once. SymPy's `Or` is
    used; its automatic simplification is classically sound.
  * The trusted base includes SymPy's automatic evaluation, e.g.
    `Implies(A, A)` is `True` at once.
  * In 0.5.0, an unsound behavior of `remove_trivial_assumptions` was
    fixed: bound variables of a premise were dropped and could be
    unified with outer metavariables.
