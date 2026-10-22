# sympy_matrix_tools

<!-- Time-stamp: "2026-10-22T16:34:30Z" -->

Some tools for SymPy matrices.

**Status (2026-10):** sympy_matrix_tools is in maintenance mode. Version
0.4.0 supports SymPy 1.14 and implements the chain rule for
`MatrixFunction`. Version 0.5.0 builds the first-order logic on the
minimal core of the logic tools. As of SymPy 1.14, SymPy still has nothing like
`MatrixFunction` (an undefined matrix-valued function), `MatrixWild`,
`MatSum` and so on.

## Download in Google Colab

You need SymPy >= 1.14. With Google Colab, you should download SymPy >= 1.14:

```python
!pip install "sympy>=1.14"
```

and then:

``` python
!pip install git+https://github.com/JRF-2018/sympy_matrix_tools
```

To fix the version of sympy_matrix_tools, specify a tag like below:

``` python
!pip install git+https://github.com/JRF-2018/sympy_matrix_tools@v0.4.0
```

The examples below are tested in SymPy=1.14.0.


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
version 1.14.0) works like below.

``` python
>>> from sympy.stats import Expectation
>>> from sympy.stats.rv import RandomSymbol, RandomMatrixSymbol

>>> epsilon = RandomSymbol("epsilon")
>>> zeta = RandomMatrixSymbol("zeta", 2, 1)
>>> A = MatrixSymbol("A", 2, 2)

>>> Expectation(A * zeta).expand()
A*ExpectationMatrix(zeta)
>>> Expectation(zeta.T * A).expand()
ExpectationMatrix(zeta.T*A)
>>> Expectation(epsilon * Identity(1)).expand()
ExpectationMatrix(epsilon*I)

```

To fix it...

``` python
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

The original SymPy errs like below:

```python
>>> Lambda(M1, M1 + Identity(N)).as_dummy().doit()
Traceback (most recent call last):
  ...
TypeError: Mix of Matrix and Scalar symbols

````

To fix it...

```python
>>> import sympy_matrix_tools
>>> sympy_matrix_tools.fix_Basic_MatrixSymbol_about_dummy()
>>> Lambda(M1, M1 + Identity(N)).as_dummy().doit()
Lambda(_0, I + _0)

```


## Usage of functions for matrices

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


## Usage of MatrixFunction (**Experimental**)

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

  * `sympy_matrix_tools_chain_rule_py.py`, which was a prototype of
    the chain rule written by Gemini, is not used anymore. It assumed
    the wrong index order of `array_derive`.

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


## Usage of functions for sequences

```python
>>> f = Function("f")
>>> g = Function("g")
>>> n = Symbol("n", integer=True)
>>> m = Symbol("m", integer=True)
>>> m2 = Symbol("m2", integer=True)
>>> T = Symbol("T", integer=True)
>>> tau = Symbol("tau", integer=True)

```

`atoms_list`, `nth_Sum`, `Sum_expand`, `Sum_swap`, `Sum_collect`:

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
>>> Sum_collect(z)
Sum(f(_xi, 0) + g(_xi) + Sum(f(n, _xi), (n, 0, T)), (_xi, 0, tau))

```


When collectiong with level >= 2, matching the names of bound symbols is required. A bound symbol can be replaced as follows. But the example is not appropriate, because it has already been collected...

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

Partially expand symbolic series.

`Sum_step_forward`, `Sum_step_backward`:

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

I intentionally omitted infinite checks, etc.

`Sum_coeff_mul`:

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

Mixed Example:

```python
>>> Mf = MatrixFunction("Mf", n, n)
>>> z = Sum(Mf(m), (m, 0, 2))
>>> Sum_step_forward(z, step=2)
Mf(0) + Mf(1) + Mf(2)
>>> z = Sum(Mf(m), (m, 0, T))
>>> z
Sum(Mf(m), (m, 0, T))
>>> # However, z is not correctly typed.
>>> z.is_Matrix
False

```

If you want something like above, you should use `MatSum`.


## Usage of MatSum and MatProduct (**Experimental**)

```python
>>> Mf = MatrixFunction("Mf", n, n)
>>> z = MatSum(Mf(m), (m, 0, 2))
>>> z
MatSum(Mf(m), (m, 0, 2))
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


## Usage of MatrixWild and MatrixDummy (**Experimental**)

`MatrixWild`, `MatrixDummy`, `WildMatrixFunction`:

```python
>>> N = Symbol("N")
>>> M1 = MatrixSymbol("M1", N, N)
>>> M2 = MatrixSymbol("M2", N, N)

```

```python
>>> aw = MatrixWild("aw", N, N)
>>> aw + M1
aw_ + M1
>>> (M2 * M1).replace(aw * M1, aw * 2)
2*M2
>>> ad = MatrixDummy("ad", N, N)
>>> ad + M1
_ad + M1
>>> from sympy.abc import x, y
>>> awf = WildMatrixFunction("awf", N, N)
>>> Mf = MatrixFunction("Mf", N, N)
>>> Mf(x, y).match(awf)
{awf_: Mf(x, y)}

```


## Usage of Apply and MatApply and PredApply (**Experimental**)

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


## Usage of functions about unification and reasoning

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


## Usage of renaming bound symbols

`rename_bound_symbols_uniquely`:

```python
>>> x = Symbol("x")
>>> f = Function("f")
>>> x1 = Symbol("x1")

>>> rename_bound_symbols_uniquely(x + x1 + f(Lambda(x1, x1 + 1)))
x + x1 + f(Lambda(x2, x2 + 1))

```


## Usage of conditional_apply

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



## Usage of predicate logic (**Experimental**)

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


## First-order logic on the minimal core (**Experimental**)

From 0.5.0, `sympy_matrix_tools.fol` builds the first-order logic on
the minimal core of the predicate logic above (ForAll, Implies and
And). The other connectives are defined impredicatively as in
higher-order logic, and their rules are proved as theorems:

  * `Exists(Lambda(x, P(x)))` == ∀R. (∀x. P(x) → R) → R
  * `A | B` == ∀R. (A → R) ∧ (B → R) → R
  * `Bot` == ∀R. R, and `neg(A)` == A → Bot
  * The classical logic is given by the axiom `ccontr`: ¬¬A → A.

```python
>>> from sympy_matrix_tools.fol import *
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

Predicates are symbols applied by `PredApply`, so that `P(x, y)` is
written as `PredApply(P, x, y)`. A theorem is proved by the tactic
`auto`, which searches proofs by backtracking with iterative
deepening. It uses the higher-order pattern unification `ho_unify`
and applies the eliminations of ∃ and ∨ in assumptions (`SAFE`)
without backtracking. `formula_str` and `print_proofstate(...,
pretty=True)` show formulas in the logical notation.

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

The classical drinker paradox is also provable. The witness of the
outer ∃ remains a metavariable, which means that the domain is not
empty, as usual.

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

See `tests/test_fol.py` for more examples and non-theorems that are
rejected.

Remarks:

  * SymPy's `false` and `Not` are not used for ⊥ and ¬, because SymPy
    evaluates `Implies(A, false)` to `~A` at once. SymPy's `Or` is
    used; its automatic simplification is classically sound.
  * The trusted base includes SymPy's automatic evaluation, e.g.
    `Implies(A, A)` is `True` at once.
  * In 0.5.0, an unsound behavior of `remove_trivial_assumptions` was
    fixed: bound variables of a premise were dropped and could be
    unified with outer metavariables.


## Tests (for developers)

You need SymPy >= 1.14 and pytest.

```sh
$ python -m pytest -s --doctest-modules
$ python -m doctest README.md -v

```

`pytest.ini` limits the targets of pytest to `sympy_matrix_tools`
and `tests`.

The tests of the chain rule (`tests/test_matrix_function.py`)
substitute concrete functions for `MatrixFunction`'s and compare the
results with the direct differentiation by `derive_by_array`.


## Remark

The primary purpose of sympy_matrix_tools is for my personal use and experimental use. Of course, I would appreciate it if other people could use it, but since MatrixFunction etc. may appear at any time in the original SymPy, when using it, it may be better to specify the version of sympy_matrix_tools and SymPy. Or you can think it's enough to make a notebook like me and leave the output from time to time.

Normally, I should write my hopes in SymPy's Issues and have them incorporated, but I can't speak English, I don't have much technical skills, and I can't handle it that much.

If I am honored to have Sympy developers want to (partially) integrate these codes with the original Sympy, then feel free to do so.

The update history is logged below in Japanese (Sorry).

《sympy_matrix_tools を作った。 - JRF のひとこと》  
http://jrf.cocolog-nifty.com/statuses/2022/06/post-2c0e50.html

An actual usage example is on the following repository (in Japanese).
Its notebook is fixed to sympy_matrix_tools v0.4.0.

《JRF-2018/economy_control: 村田安雄『動的経済システムの最適制御』の検算＆シミュレーション。》  
https://github.com/JRF-2018/economy_control


## Author

JRF ( http://jrf.cocolog-nifty.com/statuses )

The author is Japanese.  I would be grateful if you could ask in Japanese.


## Acknowledgments

Version 0.4.0 (support for SymPy 1.14, the chain rule for
`MatrixFunction`, and this README) and version 0.5.0 (the first-order
logic, the higher-order pattern unification and the auto tactic) were
developed with Claude Code (Claude Opus 5.5, Anthropic).


## License

MIT License.

Note: Many parts are "quoted" from the SymPy sources.
