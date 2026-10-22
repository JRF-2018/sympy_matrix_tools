# Fixes of SymPy and known SymPy bugs

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


## Fix ExpectationMatrix.expand

This module can fix ExpectationMatrix.expand that currently (SymPy
version 1.14.0) works like below (not tested, because the fix changes
SymPy globally):

``` python
Expectation(A * zeta).expand()             # A*ExpectationMatrix(zeta)
Expectation(zeta.T * A).expand()           # ExpectationMatrix(zeta.T*A)
Expectation(epsilon * Identity(1)).expand() # ExpectationMatrix(epsilon*I)
```

To fix it...

``` python
>>> from sympy.stats import Expectation
>>> from sympy.stats.rv import RandomSymbol, RandomMatrixSymbol
>>> epsilon = RandomSymbol("epsilon")
>>> zeta = RandomMatrixSymbol("zeta", 2, 1)
>>> A = MatrixSymbol("A", 2, 2)

>>> import sympy_matrix_tools
>>> sympy_matrix_tools.fix_ExpectationMatrix_expand()

>>> Expectation(A * zeta).expand()
A*ExpectationMatrix(zeta)
>>> Expectation(zeta.T * A).expand()
ExpectationMatrix(zeta).T*A
>>> Expectation(epsilon * Identity(1)).expand()
Expectation(epsilon)*I

```


## Fix Basic and MatrixSymbol about dummy symbols

If you use logical tools of this package, you should fix Basic and
MatrixSymbol about dummy symbols.

The original SymPy errs like below (not tested, because the fix
changes SymPy globally):

```python
Lambda(M1, M1 + Identity(N)).as_dummy().doit()
# TypeError: Mix of Matrix and Scalar symbols
```

To fix it...

```python
>>> import sympy_matrix_tools
>>> sympy_matrix_tools.fix_Basic_MatrixSymbol_about_dummy()
>>> Lambda(M1, M1 + Identity(N)).as_dummy().doit()
Lambda(_0, I + _0)

```


## Fix Transpose._eval_derivative

Importing sympy_matrix_tools fixes `Transpose._eval_derivative` of
SymPy 1.14 (`fix_Transpose_eval_derivative` is called automatically).
The original forgets to transpose: `Transpose(t * A).diff(t)` returns
`A`, not `A.T`.

```python
>>> t = Symbol("t")
>>> A = MatrixSymbol("A", 3, 1)
>>> Transpose(t * A).diff(t)
A.T

```


## Known bugs of SymPy 1.14

The following bugs were found during the development of 0.4.0 and
later.

| Symptom | Workaround |
|---|---|
| `Transpose(t*A).diff(t)` returns `A` instead of `A.T`. | Fixed automatically on import (`fix_Transpose_eval_derivative`). |
| `convert_array_to_matrix` returns a wrong array, or raises `IndexError`, for some arrays of rank 4. `X.diff(Y)` uses it, e.g. `F(Y*Z).diff(Y)` for a `MatrixFunction` `F` is wrong. | Use `mat_diff`, which returns the array expression when it is not converted to a matrix. |
| `((t*A)**2).diff(t)` raises "Mix of Matrix and Scalar symbols". | Write `X*X` instead of `X**2`. |
| `(t**2*A + t*B).diff(t, 2)` raises "Mix of Matrix and Scalar symbols" (`factor_terms` in `Derivative`). | Write `.diff(t).diff(t)`. |
| `z.diff(x(t))` for an applied `MatrixFunction` `x(t)` raises "'Dummy' object has no attribute 'shape'". | Use `freeze_matrix_function` and `melt_matrix_function`. |
| `MatMul.expand()` does not expand scalar coefficients, and `x1*x1**(N - 2)` is not combined. | Use `mat_expand`. |
| `Implies(A, false)` is evaluated to `~A`, and `Lambda(x, x)` becomes `IdentityFunction` whose `args` is empty. | Handled in the logic tools. |
| `Predicate` takes only one argument. | Use `PredApply(P, x, y)` for predicates with more arguments. |
