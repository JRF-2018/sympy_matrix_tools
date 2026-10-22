# Wild symbols, dummy symbols and application

[Back to README](../README.md)

This document describes two groups of building blocks.  They are
mainly used inside the other parts of this library (the chain rule,
the unification and the logic tools), but they are also useful
directly.

  * **Matrix versions of `Wild` and `Dummy`**: `MatrixWild`,
    `MatrixDummy` and `WildMatrixFunction`.  SymPy's `Wild` and
    `Dummy` are scalars, so they cannot stand for matrices in matrix
    expressions.
  * **Application of an unknown function**: `Apply`, `MatApply` and
    `PredApply`.  They let a plain symbol act as a function, which
    becomes a real application when the symbol is replaced by a
    `Lambda` or a function.  This is what higher-order patterns and
    predicate variables need.

```python
>>> from sympy import *
>>> from sympy_matrix_tools import *
>>> N = Symbol("N")
>>> M1 = MatrixSymbol("M1", N, N)
>>> M2 = MatrixSymbol("M2", N, N)

```


## MatrixWild: a pattern variable for matrices (**Experimental**)

`MatrixWild(name, n, m, exclude=(), properties=())` is a pattern
variable that matches a matrix expression of shape `(n, m)` in `match`
and `replace`, like SymPy's `Wild` for scalars.  The shape may contain
`Wild`'s, which are matched with the shape of the expression.

SymPy's `Wild` does not work for this purpose: `(M2*M1).match(a*M1)`
with `a = Wild("a")` returns `None`, because `a*M1` is a scalar times a
matrix.

```python
>>> a = Wild("a")
>>> print((M2 * M1).match(a * M1))
None
>>> aw = MatrixWild("aw", N, N)
>>> aw + M1
aw_ + M1
>>> (M2 * M1).match(aw * M1)
{aw_: M2}
>>> (M2 * M1).replace(aw * M1, aw * 2)
2*M2

```


## MatrixDummy: a unique matrix symbol (**Experimental**)

`MatrixDummy(name, n, m)` is a matrix symbol that is distinct from any
other symbol even if it has the same name, like SymPy's `Dummy`.  It is
used for bound variables: the variable of a `Lambda` over matrices,
the variable of a `Subs` in the chain rule, and so on.  SymPy's `Dummy`
cannot be used in matrix expressions.

```python
>>> Dummy("d") + M1
Traceback (most recent call last):
  ...
TypeError: Mix of Matrix and Scalar symbols
>>> ad = MatrixDummy("ad", N, N)
>>> ad + M1
_ad + M1
>>> ad == MatrixDummy("ad", N, N)
False

```

`make_dummy(expr, name)` and `make_wild(expr, name)` make a dummy or a
wild symbol of the same kind as `expr`: a `MatrixDummy` or
`MatrixWild` of the same shape for a matrix, and a `Dummy` or `Wild`
with the same assumptions for a scalar.

```python
>>> make_dummy(M1, "X")
_X
>>> make_wild(M1, "X")
X_
>>> make_wild(Symbol("t", real=True), "s").is_real
True

```

To use `MatrixDummy` with `as_dummy` and `dummy_eq` of SymPy, call
`fix_Basic_MatrixSymbol_about_dummy()` (see [fixes.md](fixes.md)).


## WildMatrixFunction: a pattern for applied matrix functions (**Experimental**)

`WildMatrixFunction(name, n, m)` matches any applied `MatrixFunction`
whose values are `n x m` matrices, like SymPy's `WildFunction`.  It is
used to find or replace the unknown functions in an expression.

```python
>>> from sympy.abc import x, y
>>> awf = WildMatrixFunction("awf", N, N)
>>> Mf = MatrixFunction("Mf", N, N)
>>> Mf(x, y).match(awf)
{awf_: Mf(x, y)}
>>> (Mf(x) + M1).atoms(MatrixFunction)
{Mf(x)}

```


## Apply, MatApply and PredApply: application of an unknown function (**Experimental**)

In SymPy, `f(x)` needs `f` to be a function.  But in higher-order
patterns and in logic, we often need a **variable that stands for a
function**, which is applied to arguments and later replaced by a
concrete function.  For example, a theorem about any predicate `P`
contains `P(x)` with a variable `P`.

`Apply(Q, x, ...)` is the application of a symbol `Q` to arguments.  It
stays unevaluated while `Q` is a symbol, and when `Q` is replaced by a
`Lambda` (or a function), it becomes the real application (the
beta-reduction).  There are three versions by the kind of the value:

  * `Apply(Q, args...)`: a scalar value.
  * `MatApply(n, m, Q, args...)`: an `n x m` matrix value.
  * `PredApply(Q, args...)`: a boolean value (a predicate).  Since
    SymPy's `Predicate` takes only one argument, use
    `PredApply(P, x, y)` for a predicate with more arguments.

If the function symbol is a `Wild`, `Apply` matches by abstraction:
the pattern `Apply(F, x)` with `F = Wild("F")` matches `x + 1` with
`F = Lambda(x, x + 1)`.  The unifier `ho_unify` (see
[logic.md](logic.md)) uses this.

```python
>>> Q = Symbol("Q")
>>> MQ = MatrixSymbol("MQ", 3, 3)
>>> f = Function("f")
>>> Mf = MatrixFunction("Mf", 3, 3)
>>> b = Predicate("b")

```

```python
>>> Apply(Q, 3) + 3
3 + Apply(Q, 3)
>>> MatApply(3, 3, Q, 3) + MQ
MatApply(3, 3, Q, 3) + MQ
>>> And(PredApply(Q, 3), b(1))
PredApply(Q, 3) & Q.b(1)
>>> Apply(Q, 3).subs(Q, Lambda(Q, Q + 1))
4
>>> Apply(Q, 3).subs(Q, f)
f(3)
>>> MatApply(3, 3, Q, 3).subs(Q, Lambda(Q, Mf(Q)))
Mf(3)
>>> PredApply(Q, 3).subs(Q, Lambda(Q, b(Q)))
Q.b(3)

```

```python
>>> F = Wild("F")
>>> (x + 1).match(Apply(F, x))
{F_: Lambda(_x, _x + 1)}

```
