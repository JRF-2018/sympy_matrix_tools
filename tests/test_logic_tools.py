# -*- coding: utf-8 -*-
__version__ = '0.3.5' # Time-stamp: <2022-07-24T19:15:52Z>

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

    z = make_proofstate(g3)
    assert \
        str(z) \
        == "Implies(Q.ForAll(Lambda(x, Q.Q(f(x)))), FrozenGoal(Q.ForAll(Lambda(x1, Q.Q(f(x1))))))"
    print_proofstate(z)
    out, err = capfd.readouterr()
    assert \
        out == \
        "prem 0: Q.ForAll(Lambda(x, Q.Q(f(x))))\n" + \
        "gl*: Q.ForAll(Lambda(x1, Q.Q(f(x1))))\n"
    z = resolve_implications(z, p2)
    z = resolve_implications(z, p1)
    z = melt_theorem(z)
    print_proofstate(z)
    out, err = capfd.readouterr()
    assert \
        out == \
        "ForAll: x\n" + \
        "gl: Q.Q(f(x))\n"
    assert \
        z \
        == ForAll(Lambda(x, Q(f(x))))


    # *X6Y0*
    g4 = Implies(And(p1, p2), g3)
    g4 = and_implication_normal_form(g4)
    R = Wild("R")
    V = Wild("V", boolean=True)
    x1 = Wild("x1")
    p5 = ForAll(Lambda(x1, Implies(And(ForAll(Lambda(x, PredApply(R, x))),
                                       Implies(PredApply(R, x1), V)), V)))
    p5 = and_implication_normal_form(p5)
    z = make_proofstate(g4)
    # *X6Y1
    z = resolve_implications(z, p5)
    # *X6Y2*
    z = remove_trivial_assumptions(z, num=("prem", ("prem", [("not", ForAll), P])))
    # *X6Y3*
    z = resolve_implications(z, p5)
    # *X6Y4*
    z = remove_trivial_assumptions(z, num=("prem", ("prem", [("not", ForAll), Implies])))
    # *X6Y5*
    z = remove_trivial_assumptions(z, num=("prem", ("prem", [("not", Implies)])))
    # *X6Y6*
    z = remove_trivial_assumptions(z)
    # *X6Y7*
    z = melt_theorem(z)
    print_proofstate(z)
    out, err = capfd.readouterr()
    assert \
        out == \
        "ForAll: x\n" + \
        "prem 0: Q.ForAll(Lambda(x1, Q.P(x1)))\n" + \
        "prem 1: Q.ForAll(Lambda(x2, Implies(Q.P(x2), Q.Q(f(x2)))))\n" + \
        "gl: Q.Q(f(x))\n"


    # *X7*
    z = make_proofstate(g4)
    # *X7Y1*
    z = eresolve_implications(z, p5, elim_index=ForAll)
    # *X7Y2*
    z = eresolve_implications(z, p5, elim_index=ForAll)
    # *X7Y3*
    W = Wild("W", boolean=True)
    U = Wild("U", boolean=True)
    impE = Implies(And(Implies(U, V), U, Implies(V, W)), W)
    z = eresolve_implications(z, impE, elim_index=["U", "V"])
    # *X7Y4*
    z = try_remove_trivial_assumptions(z)
    # *X7Y5*
    z = melt_theorem(z)
    print_proofstate(z)
    out, err = capfd.readouterr()
    assert \
        out == \
        "ForAll: x\n" + \
        "prem 0: Q.ForAll(Lambda(x1, Q.P(x1)))\n" + \
        "prem 1: Q.ForAll(Lambda(x2, Implies(Q.P(x2), Q.Q(f(x2)))))\n" + \
        "gl: Q.Q(f(x))\n"
    
    # *X8*
    z = make_proofstate(g4)
    # *X8Y2*
    z = eresolve_implications(z, p5, elim_index=ForAll)
    # *X8Y3*
    z = eresolve_implications(z, p5, elim_index=ForAll)
    # *X8Y4*
    W = Wild("W", boolean=True)
    U = Wild("U", boolean=True)
    impE = Implies(And(Implies(U, V), U, Implies(V, W)), W)
    z = dresolve_implications(z, impE, elim_index=[Implies(U, V), "U", "V"])
    # *X8Y5*
    z = try_remove_trivial_assumptions(z, index=[("prem", PredApply), 0])
    z = try_remove_trivial_assumptions(z)
    # *X8Y6*
    z = melt_theorem(z)
    print_proofstate(z)
    out, err = capfd.readouterr()
    assert \
        out == \
        "ForAll: x\n" + \
        "prem 0: Q.ForAll(Lambda(x1, Q.P(x1)))\n" + \
        "prem 1: Q.ForAll(Lambda(x2, Implies(Q.P(x2), Q.Q(f(x2)))))\n" + \
        "gl: Q.Q(f(x))\n"
    
    # *X9*
    z = make_proofstate(g4)
    # *X9*Y1*
    z = forall_eresolve_implications(z, forall_index=0)
    # *X9*Y2*
    z = forall_eresolve_implications(z)
    # *X9*Y3*
    z = try_remove_trivial_assumptions(z, index=("premnum", [">", 2]))
    # *X9Y4*
    z = melt_theorem(z)
    print_proofstate(z)
    out, err = capfd.readouterr()
    assert \
        out == \
        "ForAll: x\n" + \
        "prem 0: Q.ForAll(Lambda(x1, Q.P(x1)))\n" + \
        "prem 1: Q.ForAll(Lambda(x2, Implies(Q.P(x2), Q.Q(f(x2)))))\n" + \
        "gl: Q.Q(f(x))\n"
    
    # *X10*
    z = make_proofstate(g4)
    z = sresolve_implications(z, p5, index=0)
    print_proofstate(z)
    out, err = capfd.readouterr()
    assert \
        out == \
        "prem 0: Q.ForAll(Lambda((x1, x2), Implies(Q.ForAll(Lambda(x6, Q.P(x6))) & Q.ForAll(Lambda(x7, Implies(Q.P(x7), Q.Q(f(x7))))), PredApply(R_, x2, x1))))\n" + \
        "prem 1: Q.ForAll(Lambda(x3, Implies(PredApply(V_, x3) & Q.ForAll(Lambda(x8, Q.P(x8))) & Q.ForAll(Lambda(x9, Implies(Q.P(x9), Q.Q(f(x9))))), Q.Q(f(x3)))))\n" + \
        "prem 2: Q.ForAll(Lambda(x4, Implies(PredApply(R_, Apply(x_, x4), x4) & Q.ForAll(Lambda(x10, Q.P(x10))) & Q.ForAll(Lambda(x11, Implies(Q.P(x11), Q.Q(f(x11))))), PredApply(V_, x4))))\n" + \
        "gl*: Q.ForAll(Lambda(x5, Implies(Q.ForAll(Lambda(x12, Q.P(x12))) & Q.ForAll(Lambda(x13, Implies(Q.P(x13), Q.Q(f(x13))))), Q.Q(f(x5)))))\n"

    # *X11*
    prems, z = make_proofstate(g4, split_prems=True)
    z = resolve_implications(z, prems[1])
    z = resolve_implications(z, prems[0])
    z = melt_theorem(z)
    print_proofstate(z)
    out, err = capfd.readouterr()
    assert \
        out == \
        "ForAll: x\n" + \
        "prem 0: Q.ForAll(Lambda(x1, Q.P(x1)))\n" + \
        "prem 1: Q.ForAll(Lambda(x2, Implies(Q.P(x2), Q.Q(f(x2)))))\n" + \
        "gl: Q.Q(f(x))\n"


    # *X12*
    prems, z = make_proofstate(impE, split_prems=True)
    z = resolve_implications(z, prems[2])
    z = resolve_implications(z, prems[1])
    z = resolve_implications(z, prems[0])
    z = melt_theorem(z)
    print_proofstate(z)
    out, err = capfd.readouterr()
    assert \
        out == \
        "prem 0: U_\n" + \
        "prem 1: Implies(U_, V_)\n" + \
        "prem 2: Implies(V_, W_)\n" + \
        "gl: W_\n"

    # *X15*
    U = Dummy("U", boolean=True)
    V = Dummy("V", boolean=True)
    W = Dummy("W", boolean=True)
    g = Implies(U, W)
    p1 = Implies(U, V)
    p2 = Implies(V, W)
    z = make_proofstate(g)
    z = resolve_implications(z, p2)
    z = eresolve_implications(z, p1)
    assert \
        str(melt_theorem(z)) \
        == "Implies(_U, _W)"
