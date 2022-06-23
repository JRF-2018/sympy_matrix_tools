# -*- coding: utf-8 -*-
__version__ = '0.2.1' # Time-stamp: <2022-06-23T12:48:16Z>

import pytest
from sympy import MatrixSymbol, Symbol, Function, Predicate, \
    Lambda, And, Implies
from sympy_matrix_tools import *


N = Symbol("N")
M1 = MatrixSymbol("M1", N, N)
M2 = MatrixSymbol("M2", N, N)


def test_matrix_wild ():
    #fix_AssocOp__matches_commutative()
    aw = MatrixWild("aw", N, N)

    assert \
        aw.subs(N, 3).shape \
        == (3, 3)

    assert \
        str(aw + M1) \
        == "aw_ + M1"

    assert \
        (M2 * M1).replace(aw * M1, aw * 2) \
        == 2 * M2

    assert \
        (M2 + M1).replace(aw, aw * 2).doit() \
        == 2*(2*M1 + 2*M2)

    assert \
        (M2 + M1).replace(aw + M1, aw * 2).doit() \
        == M2

    assert \
        str(MatrixElement(aw, 0, 0)) \
        == "aw_[0, 0]"

    assert \
        (aw.subs(N, N) == aw) \
        == True

    assert \
        (aw.subs(N, 3) == aw) \
        == False


def test_matrix_dummy ():
    ad = MatrixDummy("ad", N, N)
    ad2 = MatrixDummy("ad", 3, 3)

    assert \
        ad.subs(N, 3).shape \
        == (3, 3)

    assert \
        str(ad + M1) \
        == "_ad + M1"

    assert \
        str(MatrixElement(ad, 0, 0)) \
        == "_ad[0, 0]"

    assert \
        str(Eq(ad.subs(N, 3), ad2).simplify()) \
        == "Eq(_ad, _ad)"

    assert \
        Eq(ad.subs(N, 3), ad.subs(N, 3)).simplify() \
        == True

    assert \
        (ad.subs(N, N) == ad) \
        == True

    assert \
        (ad.subs(N, 3) == ad) \
        == False


def test_wild_matrix_function ():
    from sympy.abc import x, y
    awf = WildMatrixFunction("awf", N, N)
    Mf = MatrixFunction("Mf", N, N)

    assert x.match(awf) \
        == None

    assert \
        str(awf.match(awf)) \
        == "{awf_: awf_}"
    
    with pytest.raises(TypeError, match=r".*match*"):
        Mf.match(awf)

    assert \
        str(awf.match(awf)) \
        == "{awf_: awf_}"

    assert \
        str(Mf(x).match(awf)) \
        == "{awf_: Mf(x)}"

    assert \
        str(Mf(x, y).match(awf)) \
        == "{awf_: Mf(x, y)}"

    awf = WildMatrixFunction('awf', N, N, nargs=2)

    assert \
        awf.nargs == {2}

    assert \
        Mf(x).match(awf) \
        == None

    assert \
        str(Mf(x, y).match(awf)) \
        == "{awf_: Mf(x, y)}"

    awf = WildMatrixFunction('awf', N, N, nargs=(1, 2))

    assert \
        awf.nargs == {1, 2}

    assert \
        str(Mf(x).match(awf)) \
        == "{awf_: Mf(x)}"

    assert \
        str(Mf(x, y).match(awf)) \
        == "{awf_: Mf(x, y)}"

    assert \
        Mf(x, y, 1).match(awf) \
        == None

    awf2 = WildMatrixFunction('awf', 3, 3, nargs=(1, 2))

    assert \
        (awf == awf2) \
        == False

    awf2 = WildMatrixFunction('awf', N, N, nargs=(1, 2))

    assert \
        (awf == awf2) \
        == True


def test_mixed_1 ():
    Mf = MatrixFunction("Mf", N, N)
    Mg = MatrixFunction("Mg", N, N)
    Q = MatrixWild("Q", N, N)

    assert \
        rewrite_subterm(M1 + M2 * Mf(M2), Q * Mf(Q), Mg(2 * Q), variables=[Q]) \
        == Mg(2 * M2) + M1


def test_matches ():
    a = Wild("a")
    M1 = MatrixDummy("M1", a, a, 1)
    M1b = MatrixDummy("M1", 3, 3, 1)
    assert \
        str(M1.matches(M1b)) \
        == "{a_: 3}"

    M1 = MatrixWild("M1", a, a)
    M1b = MatrixWild("M1", 3, 3)
    assert \
        str(M1.matches(M1b)) \
        == "{a_: 3, M1_: M1_}"


def test_xreplace ():
    a = Wild("a")
    A = MatrixWild("A", a, a)
    assert \
        A.xreplace({a: 3}).shape \
        == (3, 3)
    D = MatrixDummy("D", a, a)
    assert \
        D.xreplace({a: 3}).shape \
        == (3, 3)


def test_as_dummy ():
    fix_Basic_MatrixSymbol_about_dummy()
    N = Symbol("N")
    x = MatrixSymbol("x", N, N)
    y = MatrixSymbol("y", N, N)
    f = MatrixFunction("f", N, N)
    g = MatrixFunction("g", N, N)
    q = Lambda(x, y + f(x) + x)

    assert \
        str(q.as_dummy().doit(deep=True)) \
        == "Lambda(_0, f(_0) + _0 + y)"

    assert \
        str(x.as_dummy()) \
        == "_x"
