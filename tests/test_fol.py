# -*- coding: utf-8 -*-
__version__ = '0.5.0' # Time-stamp: <2026-10-22T16:34:30Z>

import pytest
from sympy import Symbol, Lambda, Implies, And, Or, Equivalent, symbols
from sympy_matrix_tools import PredApply, ForAll
from sympy_matrix_tools.fol import *

x, y = symbols("x y")
P = Symbol("P")
Q = Symbol("Q")
R = Symbol("R")
D = Symbol("D")
A, B = symbols("A B", boolean=True)


def p (f, *args):
    return PredApply(f, *args)


def E (v, body):
    return Exists(Lambda(v, body))


def F (v, body):
    return ForAll(Lambda(v, body))


def prove (name, goal, *rules, **kwargs):
    kwargs.setdefault('depth', 8)
    kwargs.setdefault('max_nodes', 600)
    kwargs.setdefault('safe', SAFE)
    th = Theorem(name, goal, by_auto(*rules, **kwargs))
    th.verify()
    return th


def test_derived_rules ():
    for th in THEOREMS:
        assert th.verify()


def test_first_order_theorems ():
    prove("syllogism",
          Implies(And(F(x, Implies(p(P, x), p(Q, x))),
                      F(x, Implies(p(Q, x), p(R, x)))),
                  F(x, Implies(p(P, x), p(R, x)))))
    prove("forall_swap",
          Implies(F(x, F(y, p(P, x, y))), F(y, F(x, p(P, x, y)))))
    prove("ex_all_swap",
          Implies(E(x, F(y, p(P, x, y))), F(y, E(x, p(P, x, y)))), exI)
    prove("ex_disj_distrib1",
          Implies(E(x, Or(p(P, x), p(Q, x))),
                  Or(E(x, p(P, x)), E(x, p(Q, x)))), exI, disjI1)
    prove("ex_disj_distrib2",
          Implies(Or(E(x, p(P, x)), E(x, p(Q, x))),
                  E(x, Or(p(P, x), p(Q, x)))), exI, disjI1)
    prove("all_conj",
          Equivalent(F(x, And(p(P, x), p(Q, x))),
                     And(F(x, p(P, x)), F(x, p(Q, x)))))
    prove("not_ex_all_not",
          Implies(neg(E(x, p(P, x))), F(x, neg(p(P, x)))), exI)
    prove("all_not_not_ex",
          Implies(F(x, neg(p(P, x))), neg(E(x, p(P, x)))))


def test_classical_theorems ():
    prove("de_morgan_or",
          Implies(neg(Or(A, B)), And(neg(A), neg(B))), disjI1, disjI2)
    prove("peirce", Implies(Implies(Implies(A, B), A), A), ccontr, botE)
    prove("not_all_ex_not",
          Implies(neg(F(x, p(P, x))), E(x, neg(p(P, x)))),
          ccontr, exI, botE)
    # The witness of the outer Exists remains a metavariable, which
    # means that the domain is not empty, as usual.
    prove("drinker", E(x, Implies(p(D, x), F(y, p(D, y)))),
          excluded_middle, exI, botE, ccontr, depth=12, max_nodes=4000)


def test_non_theorems ():
    with pytest.raises(ValueError):
        prove("all_ex_swap",
              Implies(F(y, E(x, p(R, x, y))), E(x, F(y, p(R, x, y)))),
              exI, ccontr, botE, max_nodes=60)
    with pytest.raises(ValueError):
        # Without the classical axiom, the drinker paradox is not
        # provable.
        prove("drinker", E(x, Implies(p(D, x), F(y, p(D, y)))),
              exI, botE, max_nodes=60)
    with pytest.raises(ValueError):
        prove("or_to_and", Implies(Or(A, B), A),
              disjI1, ccontr, botE, max_nodes=60)
