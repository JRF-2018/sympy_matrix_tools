# Unification, rewriting and Isabelle-style proofstates

[Back to README](../README.md)

This document describes the logic tools, from simple to advanced:

  1. **Rewriting with patterns**: turn equations into rewrite rules
     whose listed symbols are pattern variables, and apply them
     (`rewriterules_from_eqs`, `apply_rewriterules`, ...), and use an
     implication with given facts (`resolve_mp`).
  2. **Higher-order pattern unification** (`ho_unify`), which handles
     bound variables and function variables correctly.
  3. **Helpers**: renaming bound symbols uniquely, and applying a
     function only where some conditions hold (`conditional_apply`).
  4. **Proofstates** in the style of the proof assistant Isabelle, and
     the proof search `auto`.  First-order logic built on them is in
     [fol.md](fol.md).

These tools are experimental.  They are used to apply a theorem (an
equation or an implication with conditions) to a SymPy expression in a
controlled way; for rigorous proofs, use Lean or Isabelle.

```python
>>> from sympy import *
>>> from sympy_matrix_tools import *

```


## Rewriting with patterns

In the following functions, the symbols listed in `variables` are
pattern variables, which match any expression; the other symbols must
match themselves.  They use SymPy's `sympy.unify`.

  * `rewriterules_from_eqs(eqs, variables)` makes rewrite rules from
    equations: `Eq(lhs, rhs)` rewrites `lhs` to `rhs`.
  * `apply_rewriterules(rules, z)` applies the rules to the
    subexpressions of `z` that match them, all at once (by
    `xreplace`); the results are not rewritten again.
  * `rewrite_subterm(z, from_term, to_term, variables)` does the same
    for a single rule given directly.
  * `try_rewriterule(source, target, variables)` returns a function
    that rewrites an expression if it matches `source` and returns it
    unchanged otherwise.
  * `resolve_mp(impl, conds, variables)` applies modus ponens: it
    unifies the premises of the implication `impl` with the facts
    `conds`, instantiates the variables, and returns the implication of
    the remaining premises (or the conclusion if no premise remains).
    Premises that become `True` after the instantiation, such as
    `1 < 3`, are removed.  With `complete=True`, it raises an error if
    some premises remain.

```python
>>> Q = Symbol("Q")
>>> f = Function("f")
>>> b = Predicate("b")
>>> g = Function("g")
>>> c = Predicate("c")
>>> x = Symbol("x")

```

```python
>>> rls = rewriterules_from_eqs([Eq(g(Q), f(Q)), Eq(Q * g(Q), g(Q * Q))],
...                             variables=[Q])
>>> apply_rewriterules(rls, g(x) + x * g(x))
f(x) + g(x**2)
>>> rewrite_subterm(g(x) + x * g(x), Q * g(Q), g(Q * Q), variables=[Q])
g(x) + g(x**2)
>>> resolve_mp(Implies(And(c(x), b(x), x < 3), Implies(b(x + 1), c(x + 2))),
...            [b(1), c(1)], variables=[x])
Implies(Q.b(2), Q.c(3))
>>> resolve_mp(Implies(And(c(x), b(x), x < 3), Implies(b(x + 1), c(x + 2))),
...            [b(1), c(1), b(2)], variables=[x])
Q.c(3)

```

Be careful: a pattern variable matches itself too.  Since `Q` is the
pattern variable below, `2*g(Q)` is rewritten as if `Q` were 2:

```python
>>> rls = rewriterules_from_eqs([Eq(Q * g(Q), g(Q * Q))],
...                             variables=[Q])
>>> apply_rewriterules(rls, 2 * g(Q))
g(2**2)

```

### Example: applying a definition about functions

With `Apply`, a pattern variable can stand for a function (see
[wild_apply.md](wild_apply.md)), and a theorem about all functions
with some property can be applied to a concrete function.  Here
`LinearHomogeneous2(Q)` means that a function `Q` of two variables is
homogeneous of degree 1.  Each `Lambda` is regarded as a function (a
set), as in higher-order logic.

```python
>>> a = Symbol("a")
>>> K = Symbol("K")
>>> N = Symbol("N")
>>> y = Symbol("y")
>>> u = Symbol("u")
>>> v = Symbol("v")

```

```python
>>> LinearHomogeneous2 = Predicate("LinearHomogeneous2")
>>> lhdef = Implies(LinearHomogeneous2(Q),
...                 Eq(Apply(Q, a * K, a * N), a * Apply(Q, K, N)))
>>> lhf = LinearHomogeneous2(Lambda((x, y), f(x, y)))
>>> rls = rewriterules_from_eqs([resolve_mp(lhdef, [lhf], variables=[Q],
...                                         complete=True)],
...                             variables=[a, K, N])
>>> apply_rewriterules(rls, f(3 * u, 3 * v))
3*f(u, v)

```


## Higher-order pattern unification

`ho_unify(s, t)` yields the unifiers of `s` and `t`: dicts from the
`Wild`'s (the metavariables) to terms that make `s` and `t` equal.  It
is used by the proofstates below, and is better than `sympy.unify` in
two points:

  * Bound variables of `Lambda` are compared up to renaming
    (alpha-equivalence), and a metavariable is never bound to a term
    containing a bound variable that would escape.
  * A metavariable applied to distinct bound variables, such as
    `Apply(F, x)` with `F = Wild("F")`, is a pattern in the sense of
    Miller, and is solved by abstraction: `F = Lambda(x, ...)`.

Arguments of commutative operators (`And`, `Or`, `Add`, `Mul`) are
tried in all orders, so there may be several unifiers.

```python
>>> F = Wild("F")
>>> f = Symbol("f")
>>> list(ho_unify(Lambda(x, Apply(F, x)), Lambda(y, Apply(f, y + 1))))
[{F_: Lambda(_x, Apply(f, _x + 1))}]
>>> A = Wild("A")
>>> B = Wild("B")
>>> len(list(ho_unify(A + 2 * B, x + 2 * y)))
1

```


## Renaming bound symbols

`rename_bound_symbols_uniquely(z)` renames the bound symbols in `z`
(of `Lambda`, `Sum` and so on) so that each bound symbol has a name
different from all the other symbols.  This avoids confusion when
terms are combined or printed.

```python
>>> x = Symbol("x")
>>> f = Function("f")
>>> x1 = Symbol("x1")

>>> rename_bound_symbols_uniquely(x + x1 + f(Lambda(x1, x1 + 1)))
x + x1 + f(Lambda(x2, x2 + 1))

```


## conditional_apply

`conditional_apply(z, func, conds)` applies `func` only to the
subexpressions of `z` where the conditions `conds` are known to hold
by the context:

  * in a branch of `Piecewise` or `ITE` whose condition is in `conds`,
    or after earlier branches whose conditions are false (so that
    `Not(c)` holds);
  * in the conclusion of an `Implies` whose premises are in `conds`;
  * inside a binder (`Lambda`, `Sum`, ...) of a symbol in `conds`.

If `conds` is a list, the conditions must be met in this order; if it
is a set, in any order.  `match` restricts `func` to the
subexpressions that match it.

In the following, `func` rewrites `f(x)` to `x`.  It is applied under
the condition `Not(P(x))`: not in the first branch (where `P(x)`
holds), but in the later branches, which SymPy then merges.

```python
>>> x = Symbol("x")
>>> y = Symbol("y")
>>> f = Function("f")
>>> P = Predicate("P")
>>> Q = Predicate("Q")
>>> R = Symbol("R", bool=True)
>>> conditional_apply(Piecewise((f(x) + 1, P(x)), (f(x), Q(x)), (f(x), R)),
...                   try_rewriterule(f(x), x), [Not(P(x))])
Piecewise((f(x) + 1, Q.P(x)), (x, R | Q.Q(x)))

```


## Proofstates (**Experimental**)

This part follows the generic proof assistant Isabelle.  The
connectives of the core are `ForAll` (a predicate applied to a
`Lambda`), `Implies` and `And`; they play the role of the meta-level
of Isabelle.  Other connectives are defined on them in
[fol.md](fol.md).

  * A **proofstate** is an ordinary formula
    `Implies(And(subgoal_1, ..., subgoal_n), FrozenGoal(goal))`: if all
    subgoals are proved, the goal is proved.  `FrozenGoal` is the
    identity, which marks the goal ("`gl*`" when printed).
    `make_proofstate(goal)` makes the initial proofstate, whose only
    subgoal is the goal itself.
  * A proof is a sequence of steps that transform the proofstate.  When
    no subgoal remains, `melt_theorem` removes `FrozenGoal` and returns
    the proved theorem.
  * `Wild` symbols are **metavariables** (unknowns to be found by
    unification), like `?x` of Isabelle.  Metavariables are "lifted"
    over the bound variables of a subgoal, so that they may depend on
    them.  Specify `boolean=True` for a `Wild` or `Symbol` that stands
    for a formula.

The steps:

| Function | Isabelle | What it does |
|---|---|---|
| `resolve_implications(ps, th)` | `rule` | unify a subgoal with the conclusion of the theorem `th`, and replace the subgoal by the premises of `th` (backward reasoning) |
| `eresolve_implications(ps, th, elim_index=...)` | `erule` | like `rule`, but also unify a premise of `th` with an assumption of the subgoal and remove the assumption (elimination) |
| `dresolve_implications(ps, th, elim_index=...)` | `drule` | replace an assumption that matches a premise of `th` by the conclusion of `th` (forward reasoning) |
| `forall_eresolve_implications(ps)` | `allE` | instantiate a `ForAll` assumption with a new metavariable |
| `sresolve_implications(ps, th)` | `subgoal_tac` | add the conclusion of `th` as an assumption, with new subgoals for its premises |
| `remove_trivial_assumptions(ps)` | `assumption` | prove a subgoal by one of its assumptions (also backward through assumptions with premises) |
| `split_goals(ps)` | `conjI` | split subgoals whose conclusions are `And` or `Equivalent` |
| `auto(ps, rules)` | `auto` | search a proof by backtracking with the steps above |

Each `*_implications` function has an `iter_*` version that yields all
the results, and takes `index` (which subgoal) and `num` (which
result); they may be numbers or conditions such as a type or a
predicate (see `formula_match`).  `print_proofstate(ps, pretty=True)`
and `formula_str` show formulas in the logical notation, and
`parse_formula` reads that notation back (see [fol.md](fol.md)).

```python
>>> x = Symbol("x")
>>> f = Function("f")
>>> P = Predicate("P")
>>> Q = Predicate("Q")
>>> p1 = ForAll(Lambda(x, P(x)))
>>> p2 = ForAll(Lambda(x, Implies(P(x), Q(f(x)))))
>>> g3 = ForAll(Lambda(x, Q(f(x))))  # goal

>>> # At first, make proofstate with specifying a goal.
>>> z = make_proofstate(g3)
>>> z
Implies(Q.ForAll(Lambda(x, Q.Q(f(x)))), FrozenGoal(Q.ForAll(Lambda(x1, Q.Q(f(x1))))))
>>> # FrozenGoal is semantically an identity (lambda x: x).
>>> # When printing, it is marked by "*".
>>> print_proofstate(z)  # Each prem is a subgoal.
prem 0: Q.ForAll(Lambda(x, Q.Q(f(x))))
gl*: Q.ForAll(Lambda(x1, Q.Q(f(x1))))
>>> z = resolve_implications(z, p2)
>>> print_proofstate(z)
prem 0: Q.ForAll(Lambda(x, Q.P(x)))
gl*: Q.ForAll(Lambda(x1, Q.Q(f(x1))))
>>> z = resolve_implications(z, p1)
>>> print_proofstate(z)
gl*: Q.ForAll(Lambda(x, Q.Q(f(x))))
>>> # To end the proof, melt the FrozenGoal.
>>> z = melt_theorem(z)
>>> print_proofstate(z)  # "*" was removed from "gl".
ForAll: x
gl: Q.Q(f(x))
>>> z
Q.ForAll(Lambda(x, Q.Q(f(x))))

```

The same proof by `auto`, shown in the logical notation:

```python
>>> path = auto(make_proofstate(g3), [p1, p2], trace=True)
>>> for z in path:
...     print_proofstate(z, pretty=True)
prem 0: ∀x. Q(f(x))
gl*: ∀x1. Q(f(x1))
prem 0: ∀x. P(x)
gl*: ∀x1. Q(f(x1))
gl*: ∀x. Q(f(x))

```

A rule with metavariables for formulas, such as the elimination of
implication:

```python
>>> U = Wild("U", boolean=True)
>>> V = Wild("V", boolean=True)
>>> W = Wild("W", boolean=True)
>>> impE = Implies(And(Implies(U, V), U, Implies(V, W)), W)
>>> impE
Implies(U_ & (Implies(U_, V_)) & (Implies(V_, W_)), W_)

```
