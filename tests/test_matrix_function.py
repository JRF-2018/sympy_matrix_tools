# -*- coding: utf-8 -*-
__version__ = '0.4.0' # Time-stamp: <2026-10-22T14:29:04Z>

import pytest
from sympy import MatrixSymbol, Symbol, ZeroMatrix, Identity, Matrix, Dummy, Derivative, Subs, Lambda, \
    Array, Trace, Transpose, derive_by_array, simplify
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

    assert \
      str((Mf(2 * x1 + 1) + M1).diff(x1) + ZeroMatrix(N+1,N+1)) \
      == "2*Subs(Derivative(Mf(_xi_1), _xi_1), _xi_1, 2*x1 + 1)"
    assert \
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
    
    # convert_array_to_matrix of SymPy 1.14 fails.
    with pytest.raises(IndexError):
        ex.diff(u)
    assert \
        str(mat_diff(ex, u)) \
        == "ArrayAdd(PermuteDims(ArrayContraction(ArrayTensorProduct(-1, Subs(Derivative(phi(_xi_1), _xi_1), _xi_1, u), W), (2, 5)), (2 3)), PermuteDims(ArrayTensorProduct(I, I), (3)(1 2)))"


def _explicit (e):
    if hasattr(e, 'as_explicit'):
        e = e.as_explicit()
    return Array(e)


def _check_chain_rule (expr, x, sigma, n=1):
    """Compare the n-th derivative of expr by x with the derivative of
    expr where MatrixFunction's are substituted by sigma."""
    lhs = expr
    for i in range(n):
        lhs = mat_diff(lhs, x)
    for f, lam in sigma:
        lhs = lhs.subs(f, lam)
    lhs = _explicit(lhs.doit())
    rhs = expr
    for f, lam in sigma:
        rhs = rhs.subs(f, lam)
    rhs = _explicit(rhs.doit())
    for i in range(n):
        if x.is_Matrix:
            rhs = derive_by_array(rhs, _explicit(x))
        else:
            rhs = rhs.diff(x)
    lhs = Array(lhs)
    rhs = Array(rhs)
    return lhs.shape == rhs.shape \
        and all(simplify(a - b) == 0 for a, b
                in zip(lhs.reshape(len(lhs)), rhs.reshape(len(rhs))))


def test_chain_rule ():
    X = MatrixSymbol('X', 2, 2)
    Y = MatrixSymbol('Y', 2, 2)
    A = MatrixSymbol('A', 2, 2)
    D = MatrixSymbol('D', 2, 2)
    D2 = MatrixSymbol('D2', 2, 2)
    u = MatrixSymbol('u', 2, 1)
    v = MatrixSymbol('v', 2, 1)
    t = Symbol('t')
    s = Symbol('s')
    F = MatrixFunction('F', 2, 2)
    G = MatrixFunction('G', 2, 2)
    F2 = MatrixFunction('F2', 2, 2)
    phi = MatrixFunction('phi', 2, 1)
    LF = Lambda(D, D * A * D)
    LG = Lambda(D, D.T * D + A)
    LF2 = Lambda((D, D2), D * D2 + A * D2)
    LFs = Lambda(s, s ** 2 * A + s * Identity(2))
    Lphi = Lambda(v, v * (v.T * v))

    # by a matrix.
    assert _check_chain_rule(F(X), X, [(F, LF)])
    assert _check_chain_rule(F(X * Y), X, [(F, LF)])
    assert _check_chain_rule(F(Y * X.T), X, [(F, LF)])
    assert _check_chain_rule(F2(X, X.T), X, [(F2, LF2)])
    assert _check_chain_rule(F(G(X)), X, [(F, LF), (G, LG)])
    assert _check_chain_rule(A * F(X) * Y, X, [(F, LF)])
    assert _check_chain_rule(F(Trace(X)), X, [(F, LFs)])
    assert _check_chain_rule(u - X * phi(u), u, [(phi, Lphi)])

    # by a scalar.
    assert _check_chain_rule(F(t), t, [(F, LFs)])
    assert _check_chain_rule(F(2 * t + 1) + X, t, [(F, LFs)])
    assert _check_chain_rule(F(t ** 2), t, [(F, LFs)])
    assert _check_chain_rule(F(t * X), t, [(F, LF)])
    assert _check_chain_rule(F(G(t * X)), t, [(F, LF), (G, LG)])
    # (t**2*A + t*I).diff(t, 2) fails in SymPy 1.14, so LFs can't be used.
    assert _check_chain_rule(F(t ** 2), t, [(F, Lambda(s, s ** 3 * A))], n=2)
    assert _check_chain_rule(F(t * X), t, [(F, LF)], n=2)


def test_fix_Transpose_eval_derivative ():
    t = Symbol('t')
    A = MatrixSymbol("A", 3, 1)
    x = MatrixFunction("x", 3, 1)
    assert Transpose(t * A).diff(t) == A.T
    assert \
        str((x(t).T * x(t)).diff(t)) \
        == "x(t).T*Derivative(x(t), t) + Derivative(x(t), t).T*x(t)"
