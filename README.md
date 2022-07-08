# sympy_matrix_tools

<!-- Time-stamp: "2022-07-08T17:21:39Z" -->

Some tools for SymPy matrices.

## Download in Google Colab

You need SymPy >= 1.10.1. With Google Colab, you should download SymPy >= 1.10.1:

```python
!pip install "sympy>=1.10.1"
```

and then:

``` python
!pip install git+https://github.com/JRF-2018/sympy_matrix_tools
```


## Fix coeff

This module can fix a MatMul.args_cnc error like below...

``` python
>>> from sympy import *
>>> N = Symbol("N", integer=True)
>>> x1 = Symbol("x1")
>>> M1 = MatrixSymbol("M1", N, N)
>>> M2 = MatrixSymbol("M2", N, N)
>>> z = (M2 + 2 * (M2 + Identity(N)) * M1 + Identity(N))
>>> z.coeff(M1)
Traceback (most recent call last):
  ...
AttributeError: 'list' object has no attribute 'difference'

```

To fix it...

``` python
>>> import sympy_matrix_tools
>>> sympy_matrix_tools.fix_MatMul_args_cnc()
>>> z.coeff(M1)
2*I + 2*M2

```


## Fix ExpectationMatrix.expand

This module can fix ExpectationMatrix.expand that currently (SymPy
version 1.10.1) works like below.

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
AttributeError: 'Symbol' object has no attribute 'as_coeff_mmul'

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
M1**3 + M2*M1**2 + (M1 + M2)*M1 + (M1 + M2)*M2 + 3*M1 + M1*M2

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
>>> z
(M3*M4)**(-1)*M3*M4 + M1
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

```python
>>> z = x1 ** N * M1**N + x1 * M1 ** 2
>>> z2 = mat_collect(z, x1 *  M1 ** 2, expand_pow=2)
>>> z2
x1*((x1*x1**(N - 2))*M1**(N - 2) + I)*M1**2
>>> z2.args[1]
(x1*x1**(N - 2))*M1**(N - 2) + I
>>> partial_apply(z2, z2.args[1], lambda y: y.expand().doit())
x1*((x1**N/x1)*M1**(N - 2) + I)*M1**2
>>> mat_collect(z, x1 *  M1 ** 2, expand_pow=2, expand_inner=True)
x1*((x1**N/x1)*M1**(N - 2) + I)*M1**2
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
>>> # The following example seems to be able to differentiate well, but it doesn't really work.
>>> (Mf(2 * x1 + 1) + M1).diff(x1)
2*Subs(Derivative(Mf(_xi_1), _xi_1), _xi_1, 2*x1 + 1) + 0
>>> (Mf(2 * x1 + 1) + M1).diff(x1) + ZeroMatrix(N+1,N+1)
Traceback (most recent call last):
  ...
AttributeError: 'Mul' object has no attribute 'shape'
>>> # However, it may work in practice.
>>> (Mf(2 * x1 + 1) + M1).diff(x1).subs(Mf, Lambda(x1, x1 * M1 + M2)).doit()
2*M1

```

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


## Tests (for developers)

```sh
$ python -m pytest -s --doctest-modules
$ python -m doctest README.md -v

```


## Remark

The primary purpose of sympy_matrix_tools is for my personal use and experimental use. Of course, I would appreciate it if other people could use it, but since MatrixFunction etc. may appear at any time in the original SymPy, when using it, it may be better to specify the version of sympy_matrix_tools and SymPy. Or you can think it's enough to make a notebook like me and leave the output from time to time.

Normally, I should write my hopes in SymPy's Issues and have them incorporated, but I can't speak English, I don't have much technical skills, and I can't handle it that much.

If I am honored to have Sympy developers want to (partially) integrate these codes with the original Sympy, then feel free to do so.

The update history is logged below in Japanese (Sorry).

《sympy_matrix_tools を作った。 - JRF のひとこと》  
http://jrf.cocolog-nifty.com/statuses/2022/06/post-2c0e50.html

An actual usage example is on the following repository (in Japanese).

《JRF-2018/economy_control: 村田安雄『動的経済システムの最適制御』の検算＆シミュレーション。》  
https://github.com/JRF-2018/economy_control


## Author

JRF ( http://jrf.cocolog-nifty.com/statuses )

The author is Japanese.  I would be grateful if you could ask in Japanese.


## License

MIT License.

Note: Many parts are "quoted" from the SymPy sources.
