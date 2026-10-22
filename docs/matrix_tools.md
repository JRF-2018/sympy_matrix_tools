# Functions for matrix expressions

Coefficients, collection and division in non-commutative matrix expressions.

[Back to README](../README.md)

```python
>>> from sympy import *
>>> from sympy_matrix_tools import *

```


## Define preliminary symbols

``` python
>>> from sympy import *
>>> N = Symbol("N", integer=True)
>>> x1 = Symbol("x1")
>>> M1 = MatrixSymbol("M1", N, N)
>>> M2 = MatrixSymbol("M2", N, N)

```


## Functions for matrices

`mat_coeff`:

```python
>>> from sympy_matrix_tools import *
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

`mat_collect`:

```python
>>> z = M1 * M2 + (M2 + M1) * M1 + (M2 + M1) * (M1 ** 2) + (M2 + M1) * M2 + 3 * M1
>>> mat_collect(z, M1)
(M1 + M2)*M1**2 + (M1 + M2)*M2 + (3*I + M1 + M2)*M1 + M1*M2
>>> mat_collect(z, M1, expand_pow=True)
(M1 + M2)*M2 + ((M1 + M2)*M1 + 3*I + M1 + M2)*M1 + M1*M2

```

`partial_apply`:

```python
>> z
(M1 + M2)*M1**2 + (M1 + M2)*M1 + (M1 + M2)*M2 + 3*M1 + M1*M2
>> z.args[1]
(M1 + M2)*M1**2
>>> partial_apply(z, z.args[1], lambda y: y.expand().doit())
(M1**3 + M2*M1**2) + (M1 + M2)*M1 + (M1 + M2)*M2 + 3*M1 + M1*M2

```

`mat_divide`:

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

`mat_trivial_divide`:

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

Mixed Example:

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

In SymPy 1.14, `MatMul.expand()` doesn't expand scalar coefficients
and `x1*x1**(N - 2)` is not combined. `mat_expand` expands them and
applies `powsimp` to them. `expand_inner=True` of `mat_collect` uses
`mat_expand`.

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
