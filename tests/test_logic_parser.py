# -*- coding: utf-8 -*-
__version__ = '0.8.0' # Time-stamp: <2026-10-22T18:27:21Z>

import pytest
from sympy import Symbol, Wild, Lambda, Implies, And, Or, Equivalent, \
    Predicate, Function, symbols
from sympy_matrix_tools import PredApply, Apply, ForAll, formula_str, \
    make_proofstate, auto, and_implication_normal_form
from sympy_matrix_tools.fol import Exists, Bot, neg, AXIOMS, THEOREMS, \
    Theorem, by_auto, exI, SAFE
from sympy_matrix_tools.logic_parser import parse_formula, FormulaSyntaxError

x, y = symbols("x y")
P = Symbol("P")
Q = Symbol("Q")
A, B, C = symbols("A B C", boolean=True)


def test_basic ():
    assert parse_formula("A") == A
    assert parse_formula("P(x)") == PredApply(P, x)
    assert parse_formula("P(x, y + 1)") == PredApply(P, x, y + 1)
    assert parse_formula("A → B → C") == Implies(A, Implies(B, C))
    assert parse_formula("A ∨ B ∧ C") == Or(A, And(B, C))
    assert parse_formula("¬A ∧ B") == And(neg(A), B)
    assert parse_formula("A ↔ B") == Equivalent(A, B)
    assert parse_formula("⊥") == Bot
    assert parse_formula("∀x y. P(x, y)") == \
        ForAll(Lambda((x, y), PredApply(P, x, y)))
    assert parse_formula("∃x y. P(x, y)") == \
        Exists(Lambda(x, Exists(Lambda(y, PredApply(P, x, y)))))
    # A binder extends as far to the right as possible.
    assert parse_formula("∀x. P(x) → Q(x)") == \
        ForAll(Lambda(x, Implies(PredApply(P, x), PredApply(Q, x))))
    assert parse_formula("A → ∀x. P(x)") == \
        Implies(A, ForAll(Lambda(x, PredApply(P, x))))


def test_ascii ():
    assert parse_formula("ALL x. P(x) -> EX y. Q(y) | ~A & Bot") == \
        parse_formula("∀x. P(x) → ∃y. Q(y) ∨ ¬A ∧ ⊥")
    assert parse_formula("A <-> B") == Equivalent(A, B)


def test_metavariables_and_sympy ():
    f = Function("f")
    R = Wild("R", boolean=True)
    F = Wild("F")
    assert parse_formula("?R → ?R") == Implies(R, R)
    assert parse_formula("?F(x)") == PredApply(F, x)
    assert parse_formula("P(?x(x, y), f(y))") == \
        PredApply(P, Apply(Wild("x"), x, y), f(y))
    assert parse_formula("∀x. [x < 3] → P(x)") == \
        ForAll(Lambda(x, Implies(x < 3, PredApply(P, x))))


def test_symbols_and_predicates ():
    t = Symbol("t", real=True)
    b = Predicate("b")
    r = parse_formula("b(t) → Q(t)", symbols={"t": t},
                      predicates={"b": b})
    assert r == Implies(b(t), PredApply(Q, t))


def test_round_trip ():
    for th in AXIOMS + THEOREMS:
        s = formula_str(th.statement)
        r = and_implication_normal_form(parse_formula(s))
        assert r == th.statement, s


def test_round_trip_proofstates ():
    th = Theorem("ex_all_swap",
                 parse_formula("(∃x. ∀y. P(x, y)) → (∀y. ∃x. P(x, y))"),
                 by_auto(exI, safe=SAFE))
    assert th.verify()
    path = auto(make_proofstate(th.statement), [exI.statement],
                safe=[(r.statement, i) for r, i in SAFE], trace=True)
    for z in path:
        s = formula_str(z)
        assert formula_str(parse_formula(s)) == s


def test_errors ():
    for s in ["A →", "∀. A", "(A", "P(x", "A B"]:
        with pytest.raises(FormulaSyntaxError):
            parse_formula(s)
