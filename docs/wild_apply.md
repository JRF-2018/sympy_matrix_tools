# MatrixWild, MatrixDummy, Apply, MatApply and PredApply

[Back to README](../README.md)

```python
>>> from sympy import *
>>> from sympy_matrix_tools import *

```


## MatrixWild and MatrixDummy (**Experimental**)

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


## Apply, MatApply and PredApply (**Experimental**)

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
