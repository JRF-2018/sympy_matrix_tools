# -*- coding: utf-8 -*-
__version__ = '0.3.6' # Time-stamp: <2022-10-11T17:02:52Z>

import pytest
import sympy
from packaging.version import parse as parse_version
from sympy import MatrixSymbol, Symbol, ZeroMatrix, Identity, Matrix, Dummy, Derivative, Subs, Lambda
from sympy_matrix_tools import *

N = Symbol("N", integer=True)
N2 = Symbol("N2", integer=True)
x1 = Symbol("x1")
M1 = MatrixSymbol("M1", N + 1, N + 1)
M2 = MatrixSymbol("M2", N + 1, N + 1)
I = Identity(N + 1)
Z = ZeroMatrix(N + 1, N + 1)
f = Function("f")


def test_matrix_function ():
    Mf = MatrixFunction("Mf", N + 1, N + 1)
    assert \
      Mf(N) + M1 \
      == Mf(N) + M1

    # # The following example seems to be able to differentiate well, but it doesn't really work.
    # _xi_1 = Dummy('_xi_1')
    # assert \
    #   (Mf(2 * x1 + 1) + M1).diff(x1).doit() \
    #   == 2*Subs(Derivative(Mf(_xi_1), _xi_1), _xi_1, 2*x1 + 1)
    if parse_version(sympy.__version__) < parse_version('1.11'):
        with pytest.raises(AttributeError, match=r".no attribute.*shape.*"):
            (Mf(2 * x1 + 1) + M1).diff(x1) + ZeroMatrix(N+1,N+1)
    if parse_version(sympy.__version__) >= parse_version('1.11'):
        with pytest.raises(TypeError, match=r".*Mix of Matrix and Scalar.*"):
            (Mf(2 * x1 + 1) + M1).diff(x1) + ZeroMatrix(N+1,N+1)
        
    # # However, it may work in practice.
    if parse_version(sympy.__version__) < parse_version('1.11'):
        assert \
            (Mf(2 * x1 + 1) + M1).diff(x1).subs(Mf, Lambda(x1, x1 * M1 + M2)).doit() \
            == 2*M1
    if parse_version(sympy.__version__) >= parse_version('1.11'):
        with pytest.raises(TypeError, match=r".*Mix of Matrix and Scalar.*"):
            (Mf(2 * x1 + 1) + M1).diff(x1).subs(Mf, Lambda(x1, x1 * M1 + M2)).doit() \
                == 2*M1
    
    assert \
      (Mf(N) + M1).subs(Mf, Lambda(x1, x1 * M1 + M2)).subs(N, 2).doit() \
      == (3*M1 + M2).subs(N, 2)

    # *
    assert \
        Mf(M2).subs({N:2}).shape \
        == (3, 3)
    assert \
        Mf(M2).subs({N:2}).rows \
        == 3
    assert \
        (Mf(M2) + M1).subs(Mf, Lambda(M2, M1 + M2)).doit() \
        == 2*M1 + M2
    Mf_2 = MatrixFunction("Mf", 3, 3)
    assert Mf != Mf_2
    assert \
        Mf(2 * N + 1, x1).subs(N, 2).subs({x1: 3}).as_explicit() \
        == Matrix([[(Mf_2(5, 3))[0, 0], (Mf_2(5, 3))[0, 1], (Mf_2(5, 3))[0, 2]], [(Mf_2(5, 3))[1, 0
], (Mf_2(5, 3))[1, 1], (Mf_2(5, 3))[1, 2]], [(Mf_2(5, 3))[2, 0], (Mf_2(5, 3))[2, 1], (Mf_2
(5, 3))[2, 2]]])
    if parse_version(sympy.__version__) < parse_version('1.11'):
        z = (Mf(N) + M1).subs(N, 2).subs(Mf_2, Lambda(x1, x1 * M1 + M2)).subs(N, 2).doit()
    if parse_version(sympy.__version__) >= parse_version('1.11'):
        z = (Mf(N) + M1).subs(N, 2).subs(Mf_2, Lambda(x1, x1 * M1 + M2).subs(N, 2)).doit()
    assert \
        z \
        == (M1 + 2*M1 + M2).subs(N,2)
    assert \
        z.shape \
        == (2 + 1, 2 + 1)

    M3 = MatrixSymbol("M3", 3, 3)
    x3 = MatrixSymbol("x3", 3, 1)
    Mf3 = MatrixFunction("Mf3", 3, 1)
    assert \
        (Mf3(x3).T * M3 * Mf3(x3)).subs(Mf3(x3), x3).diff(x3).subs(x3, Mf3(x3)) \
        == M3*Mf3(x3) + M3.T*Mf3(x3)


# def test_matrix_function_0 ():
#     Mf = MatrixFunction0("Mf")

#     with pytest.raises(TypeError, match=r".*MatrixFunction needs.*"):
#         Mf(x1)

#     assert \
#         Mf(N + 1, N + 1, N) + M1 \
#         == Mf(N + 1, N + 1, N) + M1
#     assert \
#         Mf(N + 1, N + 1, 2 * N + 1, x1).subs(N, 2).subs({x1: 3}).as_explicit() \
#         == Matrix([[(Mf(3, 3, 5, 3))[0, 0], (Mf(3, 3, 5, 3))[0, 1], (Mf(3, 3, 5, 3))[0, 2]]
# , [(Mf(3, 3, 5, 3))[1, 0], (Mf(3, 3, 5, 3))[1, 1], (Mf(3, 3, 5, 3))[1, 2]], [(Mf
# (3, 3, 5, 3))[2, 0], (Mf(3, 3, 5, 3))[2, 1], (Mf(3, 3, 5, 3))[2, 2]]])

#     # # The following example doesn't work.
#     # _xi_3 = Dummy("_xi_3")
#     # assert \
#     #     (Mf(N + 1, N + 1, 2 * x1 + 1) + M1).diff(x1) \
#     #     == 2*Subs(Derivative(Mf(N + 1, N + 1, _xi_3), _xi_3), _xi_3, 2*x1 + 1) + 0

#     assert \
#         (Mf(N + 1, N + 1, N) + M1).subs(Mf, Lambda((N, N2, x1), (x1 * M1 + M2).subs(N, 2))).subs(N, 2).doit() \
#         == (M1 + 2*M1 + M2).subs(N, 2)
#     assert \
#         (Mf(N + 1, N + 1, M2) + M1).subs(N, 2).subs(Mf, Lambda((N, N2, M2), (M1 + M2).subs(N, 2))).doit() \
#         == (2*M1 + M2).subs(N, 2)
#     assert \
#         Mf(N + 1, N + 1, M2).subs({N:2}).shape \
#         == (3, 3)


def test_freeze_matrix_function ():
    A = MatrixSymbol("A", N, N)
    x = MatrixFunction("x", N, 1)
    t = Symbol("t", integer=True)

    z = x(t).T * A  * x(t)
    with pytest.raises(AttributeError, match=r"'Dummy' object has no attribute 'shape'"):
        z.diff(x(t))

    z2, sigma = freeze_matrix_function(z)
    assert \
        str(z2) \
        == "x(t).T*A*x(t)"
    assert \
        str(sigma) \
        == "{x(t): x(t)}"
    assert \
        str(atoms_list(z2, MatrixSymbol)) \
        == "[x(t), A, x(t)]"

    z3 = melt_matrix_function(z2.diff(sigma[x(t)]), sigma)
    assert \
        z3 \
        == A*x(t) + A.T*x(t)
    assert \
        str(atoms_list(z3, MatrixSymbol)) \
        == "[A, A]"


def test_wild_shape ():
    a = Wild("a")
    x = Symbol("x")
    Mf = MatrixFunction("Mf", a, a)
    Mfb = MatrixFunction("Mf", N, N)
    assert \
        str(Mf(x).matches(Mfb(x))) \
        == "{a_: N}"


def test_diff ():
    m = Symbol('m')
    W = MatrixSymbol("W", m, m)
    u = MatrixSymbol("u", m, 1)
    phi = MatrixFunction("phi", m, 1)
    r = phi(u)
    ex = u - W * r
    
    with pytest.raises(NotImplementedError):
        ex.diff(u)
