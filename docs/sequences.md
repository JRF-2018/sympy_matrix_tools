# Functions for sequences, MatSum and MatProduct

[Back to README](../README.md)

Derivations of dynamic programming and optimal control are full of
symbolic sums over time, such as `Sum(beta**(t - tau)*w(y(t)), (t, tau,
T))`: we take the first term out of the sum, shift the index, collect
sums over the same range, or pull a factor out.  SymPy has `Sum` but
few functions to transform it step by step.  This module provides
them, and also sums and products of matrices.

| Function | Purpose |
|---|---|
| `atoms_list(z, Sum)` | all `Sum`'s in `z` in order (with duplicates) |
| `nth_Sum(z, i)`, `nth_Product(z, i)` | the `i`th `Sum` (`Product`) in `z` |
| `Sum_expand(z)` | split sums of sums and flatten nested `Sum`'s |
| `Sum_swap(z)` | swap the order of two summations |
| `Sum_collect(z, level=1)` | merge sums over the same range into one |
| `Sum_step_forward(z, step=1)` | take terms out of the beginning of a sum |
| `Sum_step_backward(z, step=1)` | take terms out of the end of a sum |
| `Sum_coeff_mul(z, x)` | pull a factor `x` out of a sum |
| `Product_expand`, `Product_step_forward`, `Product_step_backward` | the same for `Product` |
| `MatSum`, `MatProduct` | `Sum` and `Product` of matrices |

These functions do not check convergence or infinite ranges.

```python
>>> from sympy import *
>>> from sympy_matrix_tools import *
>>> f = Function("f")
>>> g = Function("g")
>>> n = Symbol("n", integer=True)
>>> m = Symbol("m", integer=True)
>>> m2 = Symbol("m2", integer=True)
>>> T = Symbol("T", integer=True)
>>> tau = Symbol("tau", integer=True)

```


## Finding, expanding and swapping sums

  * `atoms_list(z, *types)` returns the subexpressions of the given
    types in the preorder of the expression tree.  Unlike SymPy's
    `atoms`, which returns a set, it keeps the order and the
    duplicates, so that a sum can be specified by its index.
  * `nth_Sum(z, i)` returns `atoms_list(z, Sum)[i]`.
  * `Sum_expand(z)` splits a sum of a sum of terms into sums of each
    term, and flattens a nested `Sum` into one `Sum` with several
    limits.
  * `Sum_swap(z, nth=None, first=1)` swaps two adjacent limits of a
    `Sum` with several limits, i.e. the order of summations.  `nth`
    selects the sum, and `first` the pair of limits.

```python
>>> z = Sum(g(m) + Sum(f(n, m), (n, 0, T)), (m, 0, tau)) + Sum(f(n, 0), (n, 0, tau))
>>> atoms_list(z, Sum)
[Sum(g(m) + Sum(f(n, m), (n, 0, T)), (m, 0, tau)), Sum(f(n, m), (n, 0, T)), Sum(f(n, 0), (n, 0, tau))]
>>> nth_Sum(z, 0)
Sum(g(m) + Sum(f(n, m), (n, 0, T)), (m, 0, tau))
>>> Sum_expand(z)
Sum(f(n, 0), (n, 0, tau)) + Sum(g(m), (m, 0, tau)) + Sum(f(n, m), (n, 0, T), (m, 0, tau))
>>> Sum_swap(Sum_expand(z))
Sum(f(n, 0), (n, 0, tau)) + Sum(g(m), (m, 0, tau)) + Sum(f(n, m), (m, 0, tau), (n, 0, T))

```


## Collecting sums

`Sum_collect(z, level=1)` merges the sums in the sum `z` whose
outermost ranges (the last limits) are the same into one sum, with a
new dummy variable `_xi` as the index.  It is the inverse of
`Sum_expand`.

```python
>>> Sum_collect(z)
Sum(f(_xi, 0) + g(_xi) + Sum(f(n, _xi), (n, 0, T)), (_xi, 0, tau))

```

With `level=2`, also the second limits from the outside must have the
same ranges, and the outer variables must have the same names.  The
dummy variable can be renamed by `replace`:

```python
>>> z = Sum(g(m) + Sum(f(n, m), (n, 0, T)), (m, 0, tau)) + Sum(f(m2, 2 * m), (m2, 0, T), (m, 0, tau))
>>> z = Sum_expand(z)
>>> z
Sum(g(m), (m, 0, tau)) + Sum(f(m2, 2*m), (m2, 0, T), (m, 0, tau)) + Sum(f(n, m), (n, 0, T), (m, 0, tau))
>>> z2 = Sum_collect(z, level=2)
>>> z2
Sum(g(m), (m, 0, tau)) + Sum(f(_xi, m) + f(_xi, 2*m), (_xi, 0, T), (m, 0, tau))
>>> xi = list(z2.atoms(Dummy))[0]
>>> z2.replace(xi, n)
Sum(g(m), (m, 0, tau)) + Sum(f(n, m) + f(n, 2*m), (n, 0, T), (m, 0, tau))

```


## Stepping a sum

`Sum_step_forward(z, step=1, begin=None, end=None, minus=False,
where=-1)` takes `step` terms out of the beginning of the sum `z`:
`Sum(f(n), (n, a, b))` becomes `f(a) + ... + f(a + step - 1) +
Sum(f(n), (n, a + step, b))`.  `Sum_step_backward` takes terms out of
the end.  This is the basic step of a Bellman equation: separate the
first period from the rest.

  * A negative `step` extends the range instead, and subtracts the
    added terms.
  * `begin` (or `end`) shifts the index of the remaining sum so that
    it starts at `begin` (or ends at `end`).
  * For a symbolic `step` such as `tau`, the terms cannot be listed, so
    the sum is split into two sums.  Since the sign of a symbolic step
    is unknown, `minus=True` tells that it is negative.
  * `where` selects the limit of a sum with several limits (the last
    one, i.e. the outermost summation, by default).

```python
>>> z = Sum(f(n + 1), (n, 0, T))
>>> Sum_step_forward(z)
f(1) + Sum(f(n + 1), (n, 1, T))
>>> Sum_step_forward(z, step=2)
f(1) + f(2) + Sum(f(n + 1), (n, 2, T))
>>> Sum_step_forward(z, step=2, begin=0)
f(1) + f(2) + Sum(f(n + 3), (n, 0, T - 2))
>>> Sum_step_forward(z, step=2, end=tau)
f(1) + f(2) + Sum(f(T + n - tau + 1), (n, -T + tau + 2, tau))
>>> Sum_step_forward(z, step=-2)
-f(-1) - f(0) + Sum(f(n + 1), (n, -2, T))
>>> Sum_step_backward(z)
f(T + 1) + Sum(f(n + 1), (n, 0, T - 1))
>>> Sum_step_backward(z, step=2)
f(T) + f(T + 1) + Sum(f(n + 1), (n, 0, T - 2))
>>> Sum_step_backward(z, step=-2)
-f(T + 2) - f(T + 3) + Sum(f(n + 1), (n, 0, T + 2))
>>> Sum_step_forward(z, step=tau, minus=False)
Sum(f(n + 1), (n, 0, tau)) + Sum(f(n + 1), (n, tau + 1, T))
>>> Sum_step_forward(z, step=-tau, minus=True)
Sum(-f(n + 1), (n, -tau, -1)) + Sum(f(n + 1), (n, -tau, T))
>>> Sum_step_backward(z, step=tau, minus=False)
Sum(f(n + 1), (n, 0, T - tau - 1)) + Sum(f(n + 1), (n, T - tau, T))
>>> Sum_step_backward(z, step=-tau, minus=True)
Sum(-f(n + 1), (n, T + 1, T + tau)) + Sum(f(n + 1), (n, 0, T + tau))
>>> z = Sum(f(n + m), (n, 0, T), (m, 0, tau))
>>> Sum_step_forward(z, where=-1) # default 'where'.
Sum(f(n), (n, 0, T)) + Sum(f(m + n), (n, 0, T), (m, 1, tau))
>>> Sum_step_forward(z, where=-2)
Sum(f(m) + Sum(f(m + n), (n, 1, T)), (m, 0, tau))

```


## Pulling a factor out

`Sum_coeff_mul(z, x, right=True, nth=0)` rewrites the `nth` sum in `z`
as `x * Sum(term / x, ...)`, which is useful to make a common factor
such as a discount factor `beta` visible.  `x` must not contain the
index of the sum.

```python
>>> w = Function("w")
>>> y = Function("y")
>>> t = Symbol("t", integer=True)
>>> beta = Symbol("beta", real=True)

>>> z = Sum_step_forward(Sum(beta ** (t - tau) * w(y(t)), (t, tau, T)))
>>> z
w(y(tau)) + Sum(beta**(t - tau)*w(y(t)), (t, tau + 1, T))
>>> z = Sum_coeff_mul(z, beta)
>>> z
beta*Sum(beta**(t - tau)*w(y(t))/beta, (t, tau + 1, T)) + w(y(tau))

```


## MatSum and MatProduct (**Experimental**)

A `Sum` of matrices is not a matrix expression in SymPy, so it cannot
be multiplied with other matrices correctly:

```python
>>> Mf = MatrixFunction("Mf", n, n)
>>> z = Sum(Mf(m), (m, 0, 2))
>>> Sum_step_forward(z, step=2)
Mf(0) + Mf(1) + Mf(2)
>>> z = Sum(Mf(m), (m, 0, T))
>>> z
Sum(Mf(m), (m, 0, T))
>>> z.is_Matrix
False

```

`MatSum(term, (i, a, b), ...)` and `MatProduct(term, (i, a, b), ...)`
are the sum and the (ordered) product of matrices, which are matrix
expressions.  The functions above work for them too.

```python
>>> z = MatSum(Mf(m), (m, 0, 2))
>>> z
MatSum(Mf(m), (m, 0, 2))
>>> z.is_Matrix
True
>>> latex(z)
'\\sum_{m=0}^{2} \\operatorname{Mf}{\\left(m \\right)}'
>>> z = MatSum(MatSum(Mf(n + m), (m, 0, T)), (m2, 0, tau))
>>> z
MatSum(Mf(m + n), (m, 0, T), (m2, 0, tau))
>>> Sum_step_forward(z)
MatSum(Mf(m + n), (m, 0, T)) + MatSum(Mf(m + n), (m, 0, T), (m2, 1, tau))
>>> z = MatProduct(Mf(m), (m, 0, T))
>>> Product_step_forward(z)
Mf(0)*MatProduct(Mf(m), (m, 1, T))

```
