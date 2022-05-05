# -*- coding: utf-8 -*-
__version__ = '0.1.3' # Time-stamp: <2022-05-05T10:09:53Z>

import pytest
from sympy import MatrixSymbol, Symbol, MatMul, Rational, Identity, Dummy
from sympy_matrix_tools import *

f = Function("f")
g = Function("g")
n = Symbol("n", integer=True)
m = Symbol("m", integer=True)
m2 = Symbol("m2", integer=True)
T = Symbol("T", integer=True)
tau = Symbol("tau", integer=True)


def test_seq_tools_1 ():
    z = Sum(g(m) + Sum(f(n, m), (n, 0, T)), (m, 0, tau)) + Sum(f(n, 0), (n, 0, tau))
    assert \
        z \
        == Sum(g(m) + Sum(f(n, m), (n, 0, T)), (m, 0, tau)) + Sum(f(n, 0), (n, 0, tau))
    assert \
        atoms_list(z, Sum) \
        == [Sum(g(m) + Sum(f(n, m), (n, 0, T)), (m, 0, tau)), Sum(f(n, m), (n, 0, T)), Sum(f(n, 0), (n, 0, tau))]
    assert \
        nth_Sum(z, 0) \
        == Sum(g(m) + Sum(f(n, m), (n, 0, T)), (m, 0, tau))
    assert \
        Sum_expand(z) \
        == Sum(f(n, 0), (n, 0, tau)) + Sum(g(m), (m, 0, tau)) + Sum(f(n, m), (n, 0, T), (m, 0, tau))
    assert \
        Sum_swap(Sum_expand(z)) \
        == Sum(f(n, 0), (n, 0, tau)) + Sum(g(m), (m, 0, tau)) + Sum(f(n, m), (m, 0, tau), (n, 0, T))

def test_Sum_collect ():
    z = Sum(g(m) + Sum(f(n, m), (n, 0, T)), (m, 0, tau)) + Sum(f(n, 0), (n, 0, tau))
    z2 = Sum_collect(z)
    xi = list(z2.atoms(Dummy))[0]
    assert \
        z2.replace(xi, n) \
        == Sum(f(n, 0) + g(n) + Sum(f(n, n), (n, 0, T)), (n, 0, tau))
    z = Sum(g(m) + Sum(f(n, m), (n, 0, T)), (m, 0, tau)) + Sum(f(m2, 2 * m), (m2, 0, T), (m, 0, tau))
    z = Sum_expand(z)
    assert \
        z \
        == Sum(g(m), (m, 0, tau)) + Sum(f(m2, 2*m), (m2, 0, T), (m, 0, tau)) + Sum(f(n, m), (n, 0, T), (m, 0, tau))
    z2 = Sum_collect(z, level=2)
    xi = list(z2.atoms(Dummy))[0]
    assert \
        z2.replace(xi, n) \
        == Sum(g(m), (m, 0, tau)) + Sum(f(n, m) + f(n, 2*m), (n, 0, T), (m, 0, tau))


def test_Sum_step_1 ():
    z = Sum(f(n + 1), (n, 0, T))
    assert \
        Sum_step_forward(z) \
        == f(1) + Sum(f(n + 1), (n, 1, T))
    assert \
        Sum_step_forward(z, step=2) \
        == f(1) + f(2) + Sum(f(n + 1), (n, 2, T))
    assert \
        Sum_step_forward(z, step=2, begin=0) \
        == f(1) + f(2) + Sum(f(n + 3), (n, 0, T - 2))
    assert \
        Sum_step_forward(z, step=2, end=tau) \
        == f(1) + f(2) + Sum(f(T + n - tau + 1), (n, -T + tau + 2, tau))
    assert \
        Sum_step_forward(z, step=-2) \
        == -f(-1) - f(0) + Sum(f(n + 1), (n, -2, T))
    assert \
        Sum_step_backward(z) \
        == f(T + 1) + Sum(f(n + 1), (n, 0, T - 1))
    assert \
        Sum_step_backward(z, step=2) \
        == f(T) + f(T + 1) + Sum(f(n + 1), (n, 0, T - 2))
    assert \
        Sum_step_backward(z, step=-2) \
        == -f(T + 2) - f(T + 3) + Sum(f(n + 1), (n, 0, T + 2))
    assert \
        Sum_step_forward(z, step=tau, minus=False) \
        == Sum(f(n + 1), (n, 0, tau)) + Sum(f(n + 1), (n, tau + 1, T))
    assert \
        Sum_step_forward(z, step=-tau, minus=True) \
        == Sum(-f(n + 1), (n, -tau, -1)) + Sum(f(n + 1), (n, -tau, T))
    assert \
        Sum_step_backward(z, step=tau, minus=False) \
        == Sum(f(n + 1), (n, 0, T - tau - 1)) + Sum(f(n + 1), (n, T - tau, T))
    assert \
        Sum_step_backward(z, step=-tau, minus=True) \
        == Sum(-f(n + 1), (n, T + 1, T + tau)) + Sum(f(n + 1), (n, 0, T + tau))
    z = Sum(f(n + m), (n, 0, T), (m, 0, tau))
    assert \
        Sum_step_forward(z, where=-1) \
        == Sum(f(n), (n, 0, T)) + Sum(f(m + n), (n, 0, T), (m, 1, tau))

    assert \
        Sum_step_forward(z, where=-2) \
        == Sum(f(m) + Sum(f(m + n), (n, 1, T)), (m, 0, tau))


def test_Sum_step_2 ():
    z = Sum(f(n * 2) + f(m + 1), (n, 0, T), (m, 0, tau))
    assert \
        Sum_step_forward(z) \
        == Sum(f(1) + f(2*n), (n, 0, T)) + Sum(f(2*n) + f(m + 1), (n, 0, T), (m, 1, tau))
    assert \
        Sum_step_forward(z, step=2) \
        == Sum(f(1) + f(2*n), (n, 0, T)) + Sum(f(2) + f(2*n), (n, 0, T)) + Sum(f(2*n) + f(m + 1), (n, 0, T), (m, 2, tau))
    assert \
        Sum_step_forward(z, step=-2) \
        == -Sum(f(-1) + f(2*n), (n, 0, T)) - Sum(f(0) + f(2*n), (n, 0, T)) + Sum(f(2*n) + f(m + 1), (n, 0, T), (m, -2, tau))
    assert \
        Sum_step_backward(z) \
        == Sum(f(2*n) + f(tau + 1), (n, 0, T)) + Sum(f(2*n) + f(m + 1), (n, 0, T), (m, 0, tau - 1))
    assert \
        Sum_step_backward(z, step=2) \
        == Sum(f(2*n) + f(tau), (n, 0, T)) + Sum(f(2*n) + f(tau + 1), (n, 0, T)) + Sum(f(2*n) + f(m + 1), (n, 0, T), (m, 0, tau - 2))
    assert \
        Sum_step_backward(z, step=-2) \
        == -Sum(f(2*n) + f(tau + 2), (n, 0, T)) - Sum(f(2*n) + f(tau + 3), (n, 0, T)) + Sum(f(2*n) + f(m + 1), (n, 0, T), (m, 0, tau + 2))
    assert \
        Sum_step_forward(z, where=-2) \
        == Sum(f(0) + f(m + 1) + Sum(f(2*n) + f(m + 1), (n, 1, T)), (m, 0, tau))

    assert \
        Sum_step_backward(z, where=-2) \
        == Sum(f(2*T) + f(m + 1) + Sum(f(2*n) + f(m + 1), (n, 0, T - 1)), (m, 0, tau))


def test_Sum_step_3 ():
    z = Sum(f(n * 2) + f(m + 1), (n, 0, T))
    assert \
        Sum_step_forward(z, step=tau, minus=False) \
        == Sum(f(2*n) + f(m + 1), (n, 0, tau)) + Sum(f(2*n) + f(m + 1), (n, tau + 1, T))
    assert \
        Sum_step_forward(z, step=-tau, minus=True) \
        == Sum(-f(2*n) - f(m + 1), (n, -tau, -1)) + Sum(f(2*n) + f(m + 1), (n, -tau, T))
    assert \
        Sum_step_backward(z, step=tau, minus=False) \
        == Sum(f(2*n) + f(m + 1), (n, 0, T - tau - 1)) + Sum(f(2*n) + f(m + 1), (n, T - tau, T))
    assert \
        Sum_step_backward(z, step=-tau, minus=True) \
        == Sum(-f(2*n) - f(m + 1), (n, T + 1, T + tau)) + Sum(f(2*n) + f(m + 1), (n, 0, T + tau))


def test_Product_step_1 ():
    z = Product(f(n * 2) + f(m + 1), (n, 0, T), (m, 0, tau))
    assert \
        Product_step_forward(z) \
        == Product(f(1) + f(2*n), (n, 0, T))*Product(f(2*n) + f(m + 1), (n, 0, T), (m, 1, tau))
    assert \
        Product_step_backward(z) \
        == Product(f(2*n) + f(tau + 1), (n, 0, T))*Product(f(2*n) + f(m + 1), (n, 0, T), (m, 0, tau - 1))
    assert \
        Product_step_forward(z, where=-2) \
        == Product((f(0) + f(m + 1))*Product(f(2*n) + f(m + 1), (n, 1, T)), (m, 0, tau))
    assert \
        Product_step_backward(z, where=-2) \
        == Product((f(2*T) + f(m + 1))*Product(f(2*n) + f(m + 1), (n, 0, T - 1)), (m, 0, tau))


def test_Product_step_2 ():
    z = Product(f(n * 2) + f(m + 1), (n, 0, T))
    assert \
        Product_step_forward(z, step=tau, minus=False) \
        == Product(f(2*n) + f(m + 1), (n, 0, tau))*Product(f(2*n) + f(m + 1), (n, tau + 1, T))
    assert \
        Product_step_forward(z, step=-tau, minus=True) \
        == Product(1/(f(2*n) + f(m + 1)), (n, -tau, -1))*Product(f(2*n) + f(m + 1), (n, -tau, T))
    assert \
        Product_step_backward(z, step=tau, minus=False) \
        == Product(f(2*n) + f(m + 1), (n, 0, T - tau - 1))*Product(f(2*n) + f(m + 1), (n, T - tau, T))
    assert \
        Product_step_backward(z, step=-tau, minus=True) \
        == Product(1/(f(2*n) + f(m + 1)), (n, T + 1, T + tau))*Product(f(2*n) + f(m + 1), (n, 0, T + tau))


def test_Sum_MatrixFunction ():
    Mf = MatrixFunction("Mf", n, n)
    z = Sum(Mf(m), (m, 0, 2))
    assert \
        z \
        == Sum(Mf(m), (m, 0, 2))
    assert \
        Sum_step_forward(z, step=2) \
        == Mf(0) + Mf(1) + Mf(2)
    assert Sum_step_forward(z, step=2).is_Matrix
    z = Sum(Mf(m), (m, 0, T))
    assert \
        z \
        == Sum(Mf(m), (m, 0, T))
    # # The next is failed.
    # assert z.is_Matrix
    z2 = Sum_step_forward(z)
    # # The next is also failed
    # assert \
    #     z2 \
    #     == Mf(0) + Sum(Mf(m), (m, 1, T))
    # print([z2.args[0], z2.args[1]])
    # print([z2.args[0].is_Matrix, z2.args[1].is_Matrix])
