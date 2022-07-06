# -*- coding: utf-8 -*-
__version__ = '0.2.4' # Time-stamp: <2022-07-05T09:02:03Z>

import pytest
from sympy import MatrixSymbol, Symbol, Function, Predicate, Lambda, \
    And, Implies, ITE, Piecewise, Sum
from sympy.unify import rewriterule
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



def test_rename_bound_symbols_uniquely ():
    x = Symbol("x")
    y = Symbol("y")
    f = Function("f")
    w = Wild("w")

    x1 = Symbol("x@[1,2]@1")
    assert \
        str(clean_addressed_symbols(x + x1 + f(Lambda(x1, x1 + 1)))) \
        == "x + x1 + f(Lambda(x1, x1 + 1))"
    assert \
        str(rename_bound_symbols_uniquely(x + x1 + f(Lambda(x1, x1 + 1)))) \
        == "x + x1 + f(Lambda(x2, x2 + 1))"
    
    x1 = Symbol("@[1,2]@1")
    assert \
        str(clean_addressed_symbols(x + x1 + f(Lambda(x1, x1 + 1)))) \
        == "0 + x + f(Lambda(0, 0 + 1))"
    assert \
        str(rename_bound_symbols_uniquely(x + x1 + f(Lambda(x1, x1 + 1)))) \
        == "0 + x + f(Lambda(1, 1 + 1))"

    x1 = Symbol("x@[1,2]@1")
    assert \
        str(rename_bound_symbols_uniquely([x + x1 + f(Lambda(x1, x1 + 1)) + w,
                                           w + x1])) \
        == "[x + x1 + w_ + f(Lambda(x2, x2 + 1)), x + w1_]"
    assert \
        str(rename_bound_symbols_uniquely([[x + x1 + f(Lambda(x1, x1 + 1)) + f(Lambda(x1, 2 * x1 + 1)) + w,
            w + x1], [x + x1 + f(Lambda(x1, x1 + 1)) + w]])) \
        == "[[x + x1 + w_ + f(Lambda(x2, x2 + 1)) + f(Lambda(x3, 2*x3 + 1)), x + w_], [x + x1 + w1_ + f(Lambda(x4, x4 + 1))]]"

    assert \
        str(rename_bound_symbols_uniquely(x + Sum(f(x), (x, 0, x)))) \
        == "x + Sum(f(x1), (x1, 0, x))"

    x1 = Symbol("x1")
    assert \
        str(rename_bound_symbols_uniquely(x + x1 + f(Lambda(x1, x1 + 1)))) \
        == "x + x1 + f(Lambda(x2, x2 + 1))"


def test_conditional_apply ():
    x = Symbol("x")
    y = Symbol("y")
    f = Function("f")
    P = Predicate("P")
    Q = Predicate("Q")
    R = Symbol("R", bool=True)

    assert \
        conditional_apply(ITE(P(x), Q(x), R),
                          lambda q: next(rewriterule(R, P(x))(q)),
                          [Not(P(x))]) \
        == ITE(P(x), Q(x), P(x))
    assert \
        conditional_apply(Piecewise((f(x) + 1, P(x)), (f(x), Q(x)), (f(x), R)),
                          lambda q: next(rewriterule(f(x), x)(q)), [Not(P(x))]) \
        == Piecewise((f(x) + 1, P(x)), (x, R | Q(x)))
    assert \
        conditional_apply(Implies(And(P(x), Q(x)), R),
                          lambda q: next(rewriterule(R, P(x))(q)), [P(x)]) \
        == Implies(P(x) & Q(x), P(x))
    assert \
        conditional_apply(Lambda(x, f(x)),
                          lambda q: next(rewriterule(f(x), x + 1)(q)), [x]) \
        == Lambda(x, x + 1)
    assert \
        conditional_apply(Sum(f(x), (x, 0, y)),
                          lambda q: next(rewriterule(f(x), x + 1)(q)), [x]) \
        == Sum(x + 1, (x, 0, y))
    assert \
        conditional_apply(Sum(f(x), (x, 0, y)),
                          try_rewriterule(f(x), x + 1), [x]) \
        == Sum(x + 1, (x, 0, y))
