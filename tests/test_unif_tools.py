# -*- coding: utf-8 -*-
__version__ = '0.2.0' # Time-stamp: <2022-06-19T01:28:38Z>

import pytest
from sympy import MatrixSymbol, Symbol, Function, Predicate, Lambda, And, Implies
from sympy_matrix_tools import *

Q = Symbol("Q")
f = Function("f")
b = Predicate("b")
g = Function("g")
c = Predicate("c")
x = Symbol("x")
a = Symbol("a")
K = Symbol("K")
N = Symbol("N")
y = Symbol("y")
u = Symbol("u")
v = Symbol("v")

def test_unif_tools ():
    rls = rewriterules_from_eqs([Eq(g(Q), f(Q)), Eq(Q * g(Q), g(Q * Q))], variables=[Q])
    assert \
        apply_rewriterules(rls, g(x) + x * g(x)) \
        == f(x) + g(x**2)
    assert \
        rewrite_subterm(g(x) + x * g(x), Q * g(Q), g(Q * Q), variables=[Q]) \
        == g(x) + g(x**2)
    assert \
        resolve_mp(Implies(And(c(x), b(x), x < 3), Implies(b(x + 1), c(x + 2))),
                   [b(1), c(1)], variables=[x]) \
        == Implies(b(2), c(3))
    assert \
        resolve_mp(Implies(And(c(x), b(x), x < 3), Implies(b(x + 1), c(x + 2))),
                   [b(1), c(1), b(2)], variables=[x]) \
        == c(3)


def test_linear_homogeneous ():
    LinearHomogeneous2 = Predicate("LinearHomogeneous2")
    lhdef = Implies(LinearHomogeneous2(Q),
                    Eq(Apply(Q, a * K, a * N), a * Apply(Q, K, N)))
    lhf = LinearHomogeneous2(Lambda((x, y), f(x, y)))
    rls = rewriterules_from_eqs([resolve_mp(lhdef, [lhf], variables=[Q],
                                            complete=True)],
                                variables=[a, K, N])
    assert \
        apply_rewriterules(rls, f(3 * u, 3 * v)) \
        == 3 * f(u, v)
