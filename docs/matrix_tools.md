# Functions for matrix expressions

[Back to README](../README.md)

SymPy can collect and factor commutative expressions (`collect`,
`coeff`, `factor`), but in a matrix expression the order of factors
matters, and those functions do not work well.  This module provides
the corresponding operations for **non-commutative matrix expressions**:

| Function | Purpose |
|---|---|
| `mat_coeff(X, Y)` | the coefficient `C` of one term `C*Y` (or `Y*C`) of `X` |
| `mat_collect(X, Y)` | collect all terms `C*Y` of `X` into `(sum of C)*Y` |
| `mat_divide(X, Y)` | multiply `X` by `Y**-1` and cancel `Y` where possible |
| `mat_trivial_divide(X)` | cancel products like `(A*B)**-1*A*B` |
| `mat_expand(X)` | `expand()` that also expands scalar coefficients |
| `partial_apply(X, arg, f)` | apply `f` only to the subexpression `arg` |

They are typically used to transform the result of a derivation by
hand: collect the terms of the control variable, divide by a matrix,
and so on.  The functions for `sympy.stats` (`fix_ExpectationMatrix_expand`)
are described in [fixes.md](fixes.md).

```python
>>> from sympy import *
>>> from sympy_matrix_tools import *
>>> N = Symbol("N", integer=True)
>>> x1 = Symbol("x1")
>>> M1 = MatrixSymbol("M1", N, N)
>>> M2 = MatrixSymbol("M2", N, N)

```


## mat_coeff

`mat_coeff(X, Y, right=False, expand_pow=False, nth=0)` looks for a
term of the matrix sum `X` of the form `C*Y` and returns `C`.  With
`right=True` it looks for `Y*C` instead.

  * It returns the coefficient of **one** term, the `nth` matching one
    (counted from 0), not the sum of the coefficients; use
    `mat_collect` for that.  It returns `None` if no term matches.
  * `Y` may have scalar factors, which are divided out: the coefficient
    of `2*M2` in `2*x1*(x1 + 1)*M1*M2` is `x1*(x1 + 1)*M1`.
  * `Y` may be a sum such as `M1 + M2`, which is matched as a factor.
  * `expand_pow` controls how powers are seen.  If it is `True`, a
    power with an integer exponent such as `M1**2` is treated as the
    product `M1*M1`.  An integer `k` splits a power `M**N` as
    `M**(N - k)*M**k`, and a dict gives such a split for each power
    (see the examples at the end).

```python
>>> z = (M1 + M2) * M1 + 2 * x1 * (x1 + 1) * M1 * M2 + M2 * (M1 ** 2) + 3 * M2 * M1
>>> mat_coeff(z, M1)
M1 + M2
>>> mat_coeff(z, 2 * M2)
x1*(x1 + 1)*M1
>>> mat_coeff(z, M1, right=True)
2*x1*(x1 + 1)*M2
>>> mat_coeff(z, M1 + M2, right=True)
M1
>>> mat_coeff(z, M1, nth=1)
3*M2
>>> mat_coeff(z, M1 ** 2)
M2
>>> mat_coeff(z, M1, nth=1, expand_pow=True)
M2*M1

```


## mat_collect

`mat_collect(X, Y, right=False, expand_pow=False, expand_inner=False)`
collects all terms of `X` of the form `C*Y` into `(sum of C)*Y` (or
`Y*(sum of C)` with `right=True`), and leaves the other terms as they
are.  A term `c*Y` with a scalar `c` contributes `c*I`, where `I` is
the identity matrix.  `expand_pow` is the same as in `mat_coeff`, and
`expand_inner=True` applies `mat_expand` to the collected coefficient.

```python
>>> z = M1 * M2 + (M2 + M1) * M1 + (M2 + M1) * (M1 ** 2) + (M2 + M1) * M2 + 3 * M1
>>> mat_collect(z, M1)
(M1 + M2)*M1**2 + (M1 + M2)*M2 + (3*I + M1 + M2)*M1 + M1*M2
>>> mat_collect(z, M1, expand_pow=True)
(M1 + M2)*M2 + ((M1 + M2)*M1 + 3*I + M1 + M2)*M1 + M1*M2

```


## partial_apply

`partial_apply(X, arg, applyf)` replaces the subexpression `arg` of `X`
by `applyf(arg)`.  It is useful when only a part of a large expression
should be expanded or simplified.

```python
>>> z
(M1 + M2)*M1**2 + (M1 + M2)*M1 + (M1 + M2)*M2 + 3*M1 + M1*M2
>>> z.args[1]
(M1 + M2)*M1**2
>>> partial_apply(z, z.args[1], lambda y: y.expand().doit())
(M1**3 + M2*M1**2) + (M1 + M2)*M1 + (M1 + M2)*M2 + 3*M1 + M1*M2

```


## mat_divide

`mat_divide(X, Y, right=False)` multiplies `X` by `Y**-1` from the left
(from the right with `right=True`) and cancels `Y` in each term where
`Y` is the leftmost (rightmost) factor.  Where it cannot be cancelled,
`Y**-1` remains.  `Y` must be square or a scalar, and a scalar factor
of `Y` divides the coefficients.

```python
>>> mat_divide(M1 * M2, M1)
M2
>>> mat_divide(M1 * M2, M1, right=True)
M1*M2*M1**(-1)
>>> mat_divide(M1 * M2, M2)
M2**(-1)*M1*M2
>>> mat_divide(M1 * M2, M2, right=True)
M1
>>> mat_divide(M1 * M2 + 2 * x1 * (M1 + M2) * M2, M1)
(2*x1)*M1**(-1)*(M1 + M2)*M2 + M2
>>> mat_divide(M1 * M2 + 2 * x1 * (M1 + M2) * M2, 3 * x1 * (M1 + M2))
(1/(3*x1))*(M1 + M2)**(-1)*M1*M2 + (2/3)*M2

```


## mat_trivial_divide

`mat_trivial_divide(X)` cancels the products of the form
`(A*B)**-k*A*B` that SymPy leaves as they are, where `A*B` is a
product of non-square matrices.

```python
>>> N2 = Symbol("N2", integer=True)
>>> M3 = MatrixSymbol("M3", N, N2)
>>> M4 = MatrixSymbol("M4", N2, N)
>>> z = M1 + ((M3 * M4) ** -2) * M3 * M4
>>> z
(M3*M4)**(-2)*M3*M4 + M1
>>> mat_trivial_divide(z)
(M3*M4)**(-1) + M1
>>> z = M1 + ((M3 * M4) ** -1) * M3 * M4
>>> mat_trivial_divide(z)
I + M1

```


## Combining the functions

Collect the terms of `x1*M1`, take the coefficient, and divide:

```python
>>> z = M1 * M2 + (x1 ** 2) * (M2 + M1) * M1 + x1 * (M2 + M1) * (M1 ** 2) + (M2 + M1) * M2 + 3 * M1
>>> z2 = mat_collect(z, x1 * M1, expand_pow=True)
>>> z2
x1*(x1*(M1 + M2) + (M1 + M2)*M1)*M1 + (M1 + M2)*M2 + 3*M1 + M1*M2
>>> mat_coeff(z2, x1 * M1)
x1*(M1 + M2) + (M1 + M2)*M1
>>> mat_divide(z2, x1 * M1, right=True)
(3/x1)*I + x1*(M1 + M2) + (M1 + M2)*M1 + 1/x1*(M1 + M2)*M2*M1**(-1) + 1/x1*M1*M2*M1**(-1)

```


## mat_expand and symbolic powers

In SymPy 1.14, `MatMul.expand()` doesn't expand scalar coefficients,
and `x1*x1**(N - 2)` is not combined.  `mat_expand(X)` expands them and
applies `powsimp` to them.  `expand_inner=True` of `mat_collect` uses
`mat_expand`.

The following also shows `expand_pow` with symbolic exponents: `2`
splits `M1**N` as `M1**(N - 2)*M1**2`, and the dict
`{x1**N: 1, M1**N: 2}` splits `x1**N` as `x1**(N - 1)*x1` and `M1**N`
as `M1**(N - 2)*M1**2`.

```python
>>> z = x1 ** N * M1**N + x1 * M1 ** 2
>>> z2 = mat_collect(z, x1 *  M1 ** 2, expand_pow=2)
>>> z2
x1*((x1*x1**(N - 2))*M1**(N - 2) + I)*M1**2
>>> z2.args[1]
(x1*x1**(N - 2))*M1**(N - 2) + I
>>> partial_apply(z2, z2.args[1], lambda y: y.expand().doit())
x1*((x1*x1**(N - 2))*M1**(N - 2) + I)*M1**2
>>> partial_apply(z2, z2.args[1], mat_expand)
x1*(x1**(N - 1)*M1**(N - 2) + I)*M1**2
>>> mat_collect(z, x1 *  M1 ** 2, expand_pow=2, expand_inner=True)
x1*(x1**(N - 1)*M1**(N - 2) + I)*M1**2
>>> mat_collect(z, x1 *  M1 ** 2, expand_pow={x1**N: 1, M1**N: 2})
x1*(x1**(N - 1)*M1**(N - 2) + I)*M1**2

```
