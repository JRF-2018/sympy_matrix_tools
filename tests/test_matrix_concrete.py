# -*- coding: utf-8 -*-
__version__ = '0.1.0' # Time-stamp: <2022-05-04T12:46:52Z>
## Language: Japanese/UTF-8

import pytest
from sympy import MatrixSymbol, Symbol, latex
from sympy_matrix_tools import *

f = Function("f")
g = Function("g")
n = Symbol("n", integer=True)
m = Symbol("m", integer=True)
m2 = Symbol("m2", integer=True)
T = Symbol("T", integer=True)
tau = Symbol("tau", integer=True)
Mf = MatrixFunction("Mf", n, n)

def test_MatSum ():
    z = MatSum(Mf(m), (m, 0, 2))
    assert \
        z \
        == MatSum(Mf(m), (m, 0, 2))
    assert \
        latex(z) \
        == '\\sum_{m=0}^{2} \\operatorname{Mf}{\\left(m \\right)}'
    assert \
        Sum_step_forward(z, step=2) \
        == Mf(0) + Mf(1) + Mf(2)
    z = MatSum(Mf(m), (m, 0, T))
    assert \
        Sum_step_forward(z) \
        == Mf(0) + MatSum(Mf(m), (m, 1, T))
    z = MatSum(MatSum(Mf(n + m), (m, 0, T)), (n, 0, tau))
    assert \
        z \
        == MatSum(Mf(m + n), (m, 0, T), (n, 0, tau))


def test_MatProduct ():
    z = MatProduct(Mf(m), (m, 0, T))
    assert \
        Product_step_forward(z) \
        == Mf(0)*MatProduct(Mf(m), (m, 1, T))
