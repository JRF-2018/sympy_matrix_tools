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


## Remarks

  * The trusted base includes SymPy's automatic evaluation, e.g.
    `Implies(A, A)` is `True` at once.
  * `auto` is not complete; a long proof may need larger `depth` and
    `max_nodes`.  A proof step takes about 0.1 seconds.
  * In 0.5.0, an unsound behavior of `remove_trivial_assumptions` was
    fixed: bound variables of a premise were dropped and could be
    unified with outer metavariables.
