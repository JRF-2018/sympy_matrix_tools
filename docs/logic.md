# Unification, rewriting and Isabelle-style proofstates

[Back to README](../README.md)

```python
>>> from sympy import *
>>> from sympy_matrix_tools import *

```


## Unification and rewriting

`apply_rewriterules`, `rewrite_subterm`, `rewriterules_from_eqs`, `resolve_mp`:

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

A possible terrible mistake, because `Q` is used as a pattern:

```python
>>> rls = rewriterules_from_eqs([Eq(Q * g(Q), g(Q * Q))],
...                             variables=[Q])
>>> apply_rewriterules(rls, 2 * g(Q))
g(2**2)

```

Mixed Example:

Let Sympy behave like Higher Order Logic. Each Lambda() must semantically be a set as function.

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


## Renaming bound symbols

`rename_bound_symbols_uniquely`:

```python
>>> x = Symbol("x")
>>> f = Function("f")
>>> x1 = Symbol("x1")

>>> rename_bound_symbols_uniquely(x + x1 + f(Lambda(x1, x1 + 1)))
x + x1 + f(Lambda(x2, x2 + 1))

```


## conditional_apply

`conditional_apply`, `try_rewriterule`:

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

conditional_apply takes conditions which can also include bound variables.



## Predicate logic and proofstates (**Experimental**)

I introduced ForAll as a predicate.  You can prove predicate logic
through proofstates. Each proofstate is an ordinary theorem
itself. You can resolve a usual theorem that does not use FrozenGoal.

I referred to the generic proof assistant Isabelle.

`resolve_implications`, `remove_trivial_assumptions`:

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

In this framework, Wild symbols behave like meta variables in
Isabelle. You need to specify boolean=True when creating Wild or
Symbol as boolean.

```python
>>> U = Wild("U", boolean=True)
>>> V = Wild("V", boolean=True)
>>> W = Wild("W", boolean=True)
>>> impE = Implies(And(Implies(U, V), U, Implies(V, W)), W)
>>> impE
Implies(U_ & (Implies(U_, V_)) & (Implies(V_, W_)), W_)

```
