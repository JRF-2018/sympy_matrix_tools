# -*- coding: utf-8 -*-
__version__ = '0.2.4' # Time-stamp: <2022-07-06T04:23:09Z>

import pytest
from sympy import MatrixSymbol, Symbol, Function, Lambda, Wild
from sympy_matrix_tools import *


def test_and_implication_normal_form ():
    x = Symbol("x")
    y = Symbol("y")
    f = Function("f")
    P = Predicate("P")
    Q = Predicate("Q")
    R = Symbol("R", bool=True)

    z = and_implication_normal_form(ForAll(Lambda(x, Implies(P(x), ForAll(Lambda(y, Implies(Q(y), R)))))))

    assert \
        z \
        == ForAll(Lambda((x, y), Implies(P(x) & Q(y), R)))

    assert \
        and_implication_normal_form(z) \
        == z
    
    assert \
        str(make_outer_forall_wild(z)) \
        == "Implies(Q.P(x_) & Q.Q(y_), R)"

    assert \
        get_numbered_premise(z, 1) \
        == Q(y)

    assert \
        str(make_bound_symbols_wild(z)) \
        == "Q.ForAll(Lambda((x_, y_), Implies(Q.P(x_) & Q.Q(y_), R)))"

    assert \
        bound_symbols_for_free(z) \
        == {R: {x, y}}


def test_resolve_implications (capfd):
    x = Symbol("x")
    y = Symbol("y")
    f = Function("f")
    P = Predicate("P")
    Q = Predicate("Q")
    p1 = ForAll(Lambda(x, P(x)))
    p2 = ForAll(Lambda(x, Implies(P(x), Q(f(x)))))
    g3 = ForAll(Lambda(x, Q(f(x))))

    z = resolve_implications(g3, p2, goal=True)
    assert \
        str(z) \
        == "Q.ForAll(Lambda(x, Implies(Q.P(x), Q.Q(f(x)))))"

    z = resolve_implications(z, p1)
    assert \
        str(z) \
        == "Q.ForAll(Lambda(x, Q.Q(f(x))))"

    g4 = Implies(And(p1, p2), g3)
    g4 = and_implication_normal_form(g4)

    R = Wild("R")
    V = Wild("V", boolean=True)
    x1 = Wild("x1")
    p5 = ForAll(Lambda(x1, Implies(And(ForAll(Lambda(x, PredApply(R, x))),
                                       Implies(PredApply(R, x1), V)), V)))
    p5 = and_implication_normal_form(p5)
    z = resolve_implications(g4, p5, goal=True)
    print_proofstate(z)
    out, err = capfd.readouterr()
    assert \
        out == \
        "ForAll: x1\n" + \
        "prem 0: Q.ForAll(Lambda(x2, Q.P(x2)))\n" + \
        "prem 1: Q.ForAll(Lambda(x3, Implies(Q.P(x3), Q.Q(f(x3)))))\n" + \
        "prem 2: Implies(PredApply(R_, x_, x1) & Q.ForAll(Lambda(x5, Q.P(x5))) & Q.ForAll(Lambda(x6, Implies(Q.P(x6), Q.Q(f(x6))))), Q.Q(f(x1)))\n" + \
        "prem 3: Q.ForAll(Lambda(x4, Implies(Q.ForAll(Lambda(x7, Q.P(x7))) & Q.ForAll(Lambda(x8, Implies(Q.P(x8), Q.Q(f(x8))))), PredApply(R_, x4, x1))))\n" + \
        "gl: Q.Q(f(x1))\n"
    z = try_remove_trivial_assumptions(z, 3)
    z = resolve_implications(z, p5, index=2)
    z = try_remove_trivial_assumptions(z, 3, 2)
    z = remove_trivial_assumptions(z, index=2, num=1)
    z = remove_trivial_assumptions(z, index=2, num=0)
    print_proofstate(z)
    out, err = capfd.readouterr()
    assert \
        out == \
        "ForAll: x\n" + \
        "prem 0: Q.ForAll(Lambda(x1, Q.P(x1)))\n" + \
        "prem 1: Q.ForAll(Lambda(x2, Implies(Q.P(x2), Q.Q(f(x2)))))\n" + \
        "gl: Q.Q(f(x))\n"

    z = eresolve_implications(g4, p5, goal=True, elim_index=1)
    z = eresolve_implications(z, p5, index=2, elim_index=1)
    W = Wild("W", boolean=True)
    U = Wild("U", boolean=True)
    impE = Implies(And(Implies(U, V), U, Implies(V, W)), W)
    z = eresolve_implications(z, impE, index=2, elim_index=1)
    z = try_remove_trivial_assumptions(z)
    print_proofstate(z)
    out, err = capfd.readouterr()
    assert \
        out == \
        "ForAll: x\n" + \
        "prem 0: Q.ForAll(Lambda(x1, Q.P(x1)))\n" + \
        "prem 1: Q.ForAll(Lambda(x2, Implies(Q.P(x2), Q.Q(f(x2)))))\n" + \
        "gl: Q.Q(f(x))\n"

    z = eresolve_implications(g4, p5, goal=True, elim_index=1)
    z = eresolve_implications(z, p5, index=2, elim_index=1)
    W = Wild("W", boolean=True)
    U = Wild("U", boolean=True)
    impE = Implies(And(Implies(U, V), U, Implies(V, W)), W)
    z = dresolve_implications(z, impE, index=2, elim_index=1)
    z = try_remove_trivial_assumptions(z, 3)
    z = try_remove_trivial_assumptions(z)
    print_proofstate(z)
    out, err = capfd.readouterr()
    assert \
        out == \
        "ForAll: x\n" + \
        "prem 0: Q.ForAll(Lambda(x1, Q.P(x1)))\n" + \
        "prem 1: Q.ForAll(Lambda(x2, Implies(Q.P(x2), Q.Q(f(x2)))))\n" + \
        "gl: Q.Q(f(x))\n"
