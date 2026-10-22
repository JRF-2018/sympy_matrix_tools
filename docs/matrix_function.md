# MatrixFunction and the chain rule

[Back to README](../README.md)

SymPy has undefined scalar functions (`Function("f")`), but no
undefined **matrix-valued** functions.  In a derivation about dynamic
systems, we need such functions: a state vector `x(t)` at time `t`, a
matrix `P(t)` of a Riccati equation, or a function `F(X)` of a matrix
`X` whose form is not specified.  This module provides:

| Name | Purpose |
|---|---|
| `MatrixFunction(name, n, m)` | an undefined function whose values are `n x m` matrices |
| `.diff(t)` | the derivative by a scalar, by the chain rule |
| `mat_diff(expr, X)` | the derivative by a matrix, by the chain rule |
| `derivative_str(expr)` | a short notation of the derivatives: `x'(t)`, `F'(a)`, `DF(a)[d]` |
| `freeze_matrix_function`, `melt_matrix_function` | replace applied functions by matrix symbols and back |

An applied `MatrixFunction` is a matrix expression, so it can be added
and multiplied with other matrices.  To make it concrete, substitute a
`Lambda` for it and call `doit()`.  The derivatives are kept
unevaluated (as `Subs` of `Derivative`) until such a substitution.

```python
>>> from sympy import *
>>> from sympy_matrix_tools import *
>>> N = Symbol("N", integer=True)
>>> x1 = Symbol("x1")
>>> M1 = MatrixSymbol("M1", N + 1, N + 1)
>>> M2 = MatrixSymbol("M2", N + 1, N + 1)

```


## MatrixFunction (**Experimental**)

`MatrixFunction("Mf", n, m)` makes an undefined function `Mf`.  It can
be applied to any number of arguments, scalars or matrices, and the
result is an `n x m` matrix expression.  The shape may be symbolic and
is substituted with the other symbols.

```python
>>> Mf = MatrixFunction("Mf", N + 1, N + 1)
>>> Mf(N) + M1
Mf(N) + M1
>>> Mf(N).shape
(N + 1, N + 1)
>>> (Mf(N) + M1).subs(Mf, Lambda(x1, x1 * M1 + M2)).subs(N, 2).doit()
3*M1 + M2

```


## The derivative by a scalar

`.diff(t)` differentiates by a scalar `t` by the chain rule.  For a
scalar argument `a(t)`, the result is `a'(t)` times the derivative of
the function at `a(t)`, which SymPy writes as
`Subs(Derivative(Mf(xi), xi), xi, a(t))`.  `derivative_str` shows it as
`Mf'(a(t))`.

```python
>>> (Mf(2 * x1 + 1) + M1).diff(x1)
2*Subs(Derivative(Mf(_xi_1), _xi_1), _xi_1, 2*x1 + 1) + 0
>>> derivative_str((Mf(2 * x1 + 1) + M1).diff(x1))
"2*Mf'(2*x1 + 1) + 0"
>>> (Mf(2 * x1 + 1) + M1).diff(x1).subs(Mf, Lambda(x1, x1 * M1 + M2)).doit()
2*M1
>>> derivative_str(Mf(x1 ** 2).diff(x1).diff(x1))
"(4*x1**2)*Mf''(x1**2) + 2*Mf'(x1**2)"

```

When an argument is a matrix `a(t)`, the derivative of the function by
the matrix is a tensor of rank 4, which is not a matrix expression.
So the derivative by a scalar is represented by the **directional
derivative** instead: `d/dt Mf(a(t)) = d/ds Mf(a(t) + s a'(t))` at
`s = 0`, which `derivative_str` shows as `DMf(a(t))[a'(t)]`.

```python
>>> X = MatrixSymbol("X", N + 1, N + 1)
>>> Mf(x1 * M2).diff(x1)
Subs(Derivative(Mf(x1*M2 + _s_1*M2), _s_1), _s_1, 0)
>>> derivative_str(Mf(x1 * M2).diff(x1))
'DMf(x1*M2)[M2]'
>>> Mf(x1 * M2).diff(x1).subs(Mf, Lambda(X, X * M1 * X)).doit()
(2*x1)*M2*M1*M2

```

A typical use is a quadratic form of a state vector `x(t)`.  Write the
derivative of `x(t)` as `x(t).diff(t)`.

```python
>>> t = Symbol("t")
>>> A = MatrixSymbol("A", N, N)
>>> x = MatrixFunction("x", N, 1)
>>> (x(t).T * A * x(t)).diff(t)
x(t).T*A*Derivative(x(t), t) + Derivative(x(t), t).T*A*x(t)
>>> derivative_str((x(t).T * A * x(t)).diff(t))
"x(t).T*A*x'(t) + x'(t).T*A*x(t)"

```


## The derivative by a matrix

`mat_diff(expr, X)` differentiates `expr` by a matrix `X`.  The result
is an array expression whose indices of `X` come first, as in SymPy's
`derive_by_array`; it is converted to a matrix expression when
possible.  The chain rule contracts the indices of each argument of a
`MatrixFunction`.

Use `mat_diff` rather than `.diff(X)`: `.diff` converts the result by
SymPy's `convert_array_to_matrix`, which has bugs about arrays of rank
4 in SymPy 1.14.  In the following, the derivative by `.diff` is
wrong, as a comparison with the direct differentiation shows.

```python
>>> Y = MatrixSymbol("Y", 2, 2)
>>> Z = MatrixSymbol("Z", 2, 2)
>>> B = MatrixSymbol("B", 2, 2)
>>> F = MatrixFunction("F", 2, 2)
>>> mat_diff(F(Y * Z), Y)
ArrayContraction(ArrayTensorProduct(Z, I, I, Subs(Derivative(F(_xi_1), _xi_1), _xi_1, Y*Z)), (0, 5), (1, 7), (3, 6))
>>> direct = derive_by_array((Y * Z * B * Y * Z).as_explicit(), Y.as_explicit())
>>> (Array(mat_diff(F(Y * Z), Y).subs(F, Lambda(Y, Y * B * Y)).doit().as_explicit()) - direct).applyfunc(expand) == Array.zeros(2, 2, 2, 2)
True
>>> (Array(F(Y * Z).diff(Y).subs(F, Lambda(Y, Y * B * Y)).doit().as_explicit()) - direct).applyfunc(expand) == Array.zeros(2, 2, 2, 2)
False

```

For a scalar function of a matrix, the result is a matrix:

```python
>>> mat_diff(Trace(F(Y)), Y).shape
(2, 2)

```

Remarks:

  * Importing sympy_matrix_tools fixes `Transpose._eval_derivative`
    of SymPy 1.14, which forgets to transpose: originally
    `Transpose(t * A).diff(t)` returns `A`, not `A.T`
    (`fix_Transpose_eval_derivative`).
  * Other bugs of SymPy 1.14 that you may meet: `((t * A) **
    2).diff(t)` and `(t ** 2 * A + t * B).diff(t, 2)` fail with "Mix of
    Matrix and Scalar symbols".  Use `X * X` instead of `X ** 2`, and
    `.diff(t).diff(t)` instead of `.diff(t, 2)`.  See
    [fixes.md](fixes.md).
  * The tests compare the chain rule with the direct differentiation
    after substituting concrete functions, for scalar and matrix
    variables, several arguments, nested functions and second
    derivatives.  See also
    [examples/chain_rule.ipynb](../examples/chain_rule.ipynb).


## freeze_matrix_function and melt_matrix_function

SymPy cannot differentiate by an applied function such as `x(t)`
itself (it replaces `x(t)` by a scalar `Dummy` internally).
`freeze_matrix_function(z)` replaces each applied `MatrixFunction` in
`z` by a `MatrixSymbol` with the same name and shape, and returns the
new expression and the substitution.  After differentiating by the
symbol, `melt_matrix_function` restores the functions.

```python
>>> A = MatrixSymbol("A", N, N)
>>> x = MatrixFunction("x", N, 1)
>>> t = Symbol("t", integer=True)

>>> z = x(t).T * A * x(t)
>>> atoms_list(z, MatrixSymbol)
[A]
>>> z.diff(x(t))
Traceback (most recent call last):
  ...
AttributeError: 'Dummy' object has no attribute 'shape'

>>> z2, sigma = freeze_matrix_function(z)
>>> z2
x(t).T*A*x(t)
>>> sigma
{x(t): x(t)}
>>> atoms_list(z2, MatrixSymbol)
[x(t), A, x(t)]
>>> z3 = melt_matrix_function(z2.diff(sigma[x(t)]), sigma)
>>> z3
A*x(t) + A.T*x(t)
>>> atoms_list(z3, MatrixSymbol)
[A, A]

```
