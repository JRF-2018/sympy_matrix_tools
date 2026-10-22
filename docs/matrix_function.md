# MatrixFunction and the chain rule

Undefined matrix-valued functions and their derivatives.

[Back to README](../README.md)

```python
>>> from sympy import *
>>> from sympy_matrix_tools import *

```


## MatrixFunction (**Experimental**)

```python
>>> N = Symbol("N", integer=True)
>>> N2 = Symbol("N2", integer=True)
>>> x1 = Symbol("x1")
>>> M1 = MatrixSymbol("M1", N + 1, N + 1)
>>> M2 = MatrixSymbol("M2", N + 1, N + 1)

```

```python
>>> Mf = MatrixFunction("Mf", N + 1, N + 1)
>>> Mf(N) + M1
Mf(N) + M1
>>> (Mf(N) + M1).subs(Mf, Lambda(x1, x1 * M1 + M2)).subs(N, 2).doit()
3*M1 + M2

```

### Chain rule

From 0.4.0, `MatrixFunction` can be differentiated by the chain rule.

By a scalar:

```python
>>> (Mf(2 * x1 + 1) + M1).diff(x1)
2*Subs(Derivative(Mf(_xi_1), _xi_1), _xi_1, 2*x1 + 1) + 0
>>> (Mf(2 * x1 + 1) + M1).diff(x1).subs(Mf, Lambda(x1, x1 * M1 + M2)).doit()
2*M1
>>> Mf(x1 ** 2).diff(x1).diff(x1)
(4*x1**2)*Subs(Derivative(Mf(_xi_1), (_xi_1, 2)), _xi_1, x1**2) + 2*Subs(Derivative(Mf(_xi_1), _xi_1), _xi_1, x1**2)

```

When an argument is a matrix, the derivative by a scalar is
represented by a directional derivative: d/dx Mf(a(x)) = d/ds Mf(a(x)
+ s a'(x)) at s = 0.

```python
>>> X = MatrixSymbol("X", N + 1, N + 1)
>>> Mf(x1 * M2).diff(x1)
Subs(Derivative(Mf(x1*M2 + _s_1*M2), _s_1), _s_1, 0)
>>> Mf(x1 * M2).diff(x1).subs(Mf, Lambda(X, X * M1 * X)).doit()
(2*x1)*M2*M1*M2

```

```python
>>> t = Symbol("t")
>>> A = MatrixSymbol("A", N, N)
>>> x = MatrixFunction("x", N, 1)
>>> (x(t).T * A * x(t)).diff(t)
x(t).T*A*Derivative(x(t), t) + Derivative(x(t), t).T*A*x(t)

```

By a matrix, the derivative is an array expression whose indices of
the variable come first, as SymPy's `derive_by_array`. Use `mat_diff`
rather than `.diff`, because `.diff` uses `convert_array_to_matrix`,
which has bugs about arrays of rank 4 in SymPy 1.14:

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

Remarks:

  * Importing sympy_matrix_tools fixes `Transpose._eval_derivative`
    of SymPy 1.14, which forgets to transpose: originally
    `Transpose(t * A).diff(t)` returns `A`, not `A.T`
    (`fix_Transpose_eval_derivative`).

  * Other bugs of SymPy 1.14 that you may meet: `((t * A) **
    2).diff(t)` and `(t ** 2 * A + t * B).diff(t, 2)` fail with "Mix of
    Matrix and Scalar symbols". Use `X * X` instead of `X ** 2`, and
    `.diff(t).diff(t)` instead of `.diff(t, 2)`.


`freeze_matrix_function` is still needed to differentiate by `x(t)`
itself:

`freeze_matrix_function`, `melt_matrix_functioin`:

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
