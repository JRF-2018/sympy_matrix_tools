# sympy_matrix_tools

<!-- Time-stamp: "2022-04-29T05:15:24Z" -->

Some tools for sympy matrices.

## Download in Google Colab

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
  File "<stdin>", line 1, in <module>
  File "C:\Users\USER1\Anaconda3\lib\site-packages\sympy\core\expr.py", line 156
1, in coeff
    rv = d.coeff(x, right=right, _first=False)
  File "C:\Users\USER1\Anaconda3\lib\site-packages\sympy\core\expr.py", line 157
8, in coeff
    resid = margs.difference(xargs)
AttributeError: 'list' object has no attribute 'difference'
```

To fix it...

``` python
>>> import sympy_matrix_tools
>>> sympy_matrix_tools.fix_MatMul_args_cnc()
>>> z.coeff(M1)
2*I + 2*M2
```


## Usage of other functions

```python
>>> from sympy_matrix_tools import *
>>> z =   (M1 + M2) * M1 + 2 * x1 * (x1 + 1) * M1 * M2 + M2 * (M1 ** 2) + 3 * M2 * M1
>>> mat_coeff(z, M1)
M1 + M2
>>> mat_coeff(z, 2 * M2)
x1*(x1 + 1)*M1
>>> mat_coeff(z, M1, right=True)
x1*(x1 + 1)*2*M2
>>> mat_coeff(z, M1 + M2, right=True)
M1
>>> mat_coeff(z, M1, nth=1)
3*M2
>>> mat_coeff(z, M1 ** 2)
M2
>>> mat_coeff(z, M1, nth=1, expand_pow=True)
M2*M1
```

```python
>>> z = M1 * M2 + (M2 + M1) * M1 + (M2 + M1) * (M1 ** 2) + (M2 + M1) * M2 + 3 * M1
>>> mat_collect(z, M1)
(M1 + M2 + 3*I)*M1 + (M1 + M2)*M1**2 + (M1 + M2)*M2 + M1*M2
>>> mat_collect(z, M1, expand_pow=True)
(M1 + M2 + (M1 + M2)*M1 + 3*I)*M1 + (M1 + M2)*M2 + M1*M2
```

```python
>>> partial_apply(z, z.args[1], lambda y: y.expand().doit())
M1**3 + M2*M1**2 + (M1 + M2)*M1 + (M1 + M2)*M2 + 3*M1 + M1*M2
```

```python
>>> mat_divide(M1 * M2, M1)
M2
>>> mat_divide(M1 * M2, M1, right=True)
(M1*M2)*M1**(-1)
>>> mat_divide(M1 * M2, M2)
M2**(-1)*(M1*M2)
>>> mat_divide(M1 * M2, M2, right=True)
M1
>>> mat_divide(M1 * M2 + 2 * x1 * (M1 + M2) * M2, M1)
M1**(-1)*((2*x1)*(M1 + M2)*M2) + M2
>>> mat_divide(M1 * M2 + 2 * x1 * (M1 + M2) * M2, 3 * x1 * (M1 + M2))
((1/(3*x1))*(M1 + M2)**(-1))*(M1*M2) + (2/3)*M2
```

```python
>>> z = M1 * M2 + (x1 ** 2) * (M2 + M1) * M1 + x1 * (M2 + M1) * (M1 ** 2) + (M2 + M1) * M2 + 3 * M1
>>> z2 = mat_collect(z, x1 * M1, expand_pow=True)
>>> z2
(x1*(M1 + M2) + (M1 + M2)*M1)*(x1*M1) + (M1 + M2)*M2 + 3*M1 + M1*M2
>>> mat_coeff(z2, x1*M1)

Mixed Example:

```python
>>> z = M1 * M2 + (x1 ** 2) * (M2 + M1) * M1 + x1 * (M2 + M1) * (M1 ** 2) + (M2 + M1) * M2 + 3 * M1
>>> z2 = mat_collect(z, x1 * M1, expand_pow=True)
>>> z2
(x1*(M1 + M2) + (M1 + M2)*M1)*(x1*M1) + (M1 + M2)*M2 + 3*M1 + M1*M2
>>> mat_coeff(z2, x1 * M1)
x1*(M1 + M2) + (M1 + M2)*M1
>>> mat_divide(z2, x1 * M1, right=True)
((M1 + M2)*M2)*(1/x1*M1**(-1)) + (M1*M2)*(1/x1*M1**(-1)) + x1*(M1 + M2) + (M1 + M2)*M1 + 3/x1
```

## Author

JRF ( http://jrf.cocolog-nifty.com/statuses )


## License

MIT License.
