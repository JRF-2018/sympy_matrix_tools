# Functions for sequences, MatSum and MatProduct

[Back to README](../README.md)

```python
>>> from sympy import *
>>> from sympy_matrix_tools import *

```


## Functions for sequences

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


## MatSum and MatProduct (**Experimental**)

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
