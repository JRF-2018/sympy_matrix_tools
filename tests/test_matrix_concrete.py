# -*- coding: utf-8 -*-
__version__ = '0.1.3' # Time-stamp: <2022-05-05T10:09:32Z>

import pytest
from sympy import MatrixSymbol, Symbol, Matrix, \
    latex, Adjoint, Trace, adjoint, trace
from sympy_matrix_tools import *

f = Function("f")
g = Function("g")
n = Symbol("n", integer=True)
m = Symbol("m", integer=True)
m2 = Symbol("m2", integer=True)
T = Symbol("T", integer=True)
tau = Symbol("tau", integer=True)
M1 = MatrixSymbol("M1", n, n)
Mf = MatrixFunction("Mf", n, n)


def test_MatSum_1 ():
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


def test_MatSum_2 ():
    z = MatSum(Mf(m), (m, 0, T))
    assert \
        z.T \
        == MatSum(Mf(m).T, (m, 0, T))
    assert \
        adjoint(z) \
        == MatSum(Adjoint(Mf(m)), (m, 0, T))
    assert \
        str(z.subs(n, 2).as_explicit()) \
        == 'Matrix([[Sum((Mf(m))[0, 0], (m, 0, T)), Sum((Mf(m))[0, 1], (m, 0, T))], [Sum((Mf(m))[1, 0], (m, 0, T)), Sum((Mf(m))[1, 1], (m, 0, T))]])'
    assert \
        trace(z.subs(n, 2)) \
        == Sum(Trace(Mf(m).subs(n, 2)), (m, 0, T))
    assert \
        z.doit() \
        == MatSum(Mf(m), (m, 0, T))


def test_MatProduct_1 ():
    z = MatProduct(Mf(m), (m, 0, T))
    assert \
        Product_step_forward(z) \
        == Mf(0)*MatProduct(Mf(m), (m, 1, T))


def test_MatProduct_2 ():
    z = MatProduct(Mf(m), (m, 0, T))
    assert \
        z.T \
        == MatProduct(Mf(m), (m, 0, T)).T
    assert \
        adjoint(z) \
        == adjoint(MatProduct(Mf(m), (m, 0, T)))
    assert \
        z.doit() \
        == MatProduct(Mf(m), (m, 0, T))
    assert \
        str(z.subs(n, 2).as_explicit()) \
        == 'Matrix([[MatProduct(Mf(m), (m, 0, T))[0, 0], MatProduct(Mf(m), (m, 0, T))[0, 1]], [MatProduct(Mf(m), (m, 0, T))[1, 0], MatProduct(Mf(m), (m, 0, T))[1, 1]]])'


def test_Misc_1 ():
    assert \
        MatSum(m * M1, (m, 0, 2)).doit() \
        == 3*M1
    assert \
        MatSum(Mf(m), (m, 0, 2)).doit() \
        == Mf(0) + Mf(1) + Mf(2)
    assert \
        MatProduct(Mf(m), (m, 0, 2)).doit() \
        == Mf(0)*Mf(1)*Mf(2)
