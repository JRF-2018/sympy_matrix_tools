# -*- coding: utf-8 -*-
__version__ = '0.1.10' # Time-stamp: <2022-05-20T02:00:04Z>

import pytest
from sympy import MatrixSymbol, Symbol, MatMul, Rational, Identity
from sympy_matrix_tools import *

N = Symbol("N", integer=True)
N2 = Symbol("N2", integer=True)
x1 = Symbol("x1")
M1 = MatrixSymbol("M1", N, N)
M2 = MatrixSymbol("M2", N, N)
I = Identity(N)


def test_fix_MatMul_args_cnc ():
    z = M2 + 2 * (M2 + Identity(N)) * M1 + Identity(N)
    with pytest.raises(AttributeError, match=r".*list.*difference.*"):
        z.coeff(M1)
    fix_MatMul_args_cnc()
    assert z.coeff(M1) == 2*I + 2*M2


def test_fix_ExpectationMatrix_expand ():
    from sympy.stats import Expectation, ExpectationMatrix
    from sympy.stats.rv import RandomSymbol, RandomMatrixSymbol
    
    zeta = RandomMatrixSymbol("zeta", 2, 1)
    A = MatrixSymbol("A", 2, 2)
    assert \
        Expectation(zeta.T * A).expand() \
        == ExpectationMatrix(zeta.T * A)
    fix_ExpectationMatrix_expand()
    assert \
        Expectation(zeta.T * A).expand() \
        == ExpectationMatrix(zeta).T * A


def test_collect_1 ():
    assert \
        mat_collect(M1 * M2 + M2 * M2, M2, expand_pow=True) \
        == (M1 + M2)*M2
    assert \
        mat_collect(M1 * M2 + M2 * M2, M2) \
        == M2**2 + M1*M2
    assert \
        mat_collect(M1 * M2 + M2 * M2 + 3 * M2, M2, expand_pow=True) \
        == (3*I + M1 + M2)*M2
    assert \
        mat_collect(M1 * M2 + x1 * M2 ** 2 + (M1 - x1 * M2) *M2, M2, expand_pow=True) \
        == 2*M1*M2
    assert mat_collect(M1 * M2 + x1 * M2 ** 2 + (M1 - x1 * M2) *M2 + M2, M2, expand_pow=True) \
        == (I + 2*M1)*M2
    assert \
        mat_collect(mat_collect(x1 * (M1**(-1))**2 + M2 * (M1 ** -1) ** 2 + M1 * M2 + 3 * M2, M1 ** -1, expand_pow=True), M2) \
        == (3*I + M1)*M2 + (x1*M1**(-1) + M2*M1**(-1))*M1**(-1)
    z2 = mat_collect(x1 ** N * M1**(-N) + x1 * (M1 ** -2), M1 ** -1, expand_pow=-1)
    assert \
        z2 \
        == (x1*M1**(-1) + (x1**(N - 1)/x1)*M1**(-N - 1))*M1**(-1)
    z2 = mat_collect(x1 ** N * M1**N + x1 * (M1 ** N2), M1 ** N2, expand_pow=N2, expand_inner=True)
    assert \
        z2 \
        == (x1*I + x1**N*M1**(N - N2))*M1**N2
    assert \
        mat_collect(M1 * M2 + M2 ** 2, M2, expand_pow={M2 ** 2: 2}) \
        == (M1 + M2)*M2


def test_coeff ():
    z = (M1 + M2) * M1 + 2 * x1 * (x1 + 1) * M1 * M2 + M2 * (M1 ** 2) + 3 * M2 * M1
    assert \
        mat_coeff(z, M1) \
        == M1 + M2
    assert \
        mat_coeff(z, 2 * M2) \
        == MatMul(x1, (x1 + 1), M1)
    assert \
        mat_coeff(z, M1, right=True) \
        == MatMul(2, x1, (x1 + 1), M2)
    assert \
        mat_coeff(z, M1 + M2, right=True) \
        == M1
    assert \
        mat_coeff(z, M1, nth=1) \
        == MatMul(3, M2)
    assert \
        mat_coeff(z, M1 ** 2) \
        == M2
    assert \
        mat_coeff(z, M1, nth=1, expand_pow=True) \
        == M2*M1
    # *
    assert \
        mat_coeff(-M1, M1) \
        == -1
    assert \
        mat_coeff(z, 2 * M1, right=True) \
        == MatMul(x1, (x1 + 1), M2)
    assert \
        mat_coeff(2 * x1 * M1, x1 * M1) \
        == 2


def test_collect_2 ():
    z = M1 * M2 + (M2 + M1) * M1 + (M2 + M1) * (M1 ** 2) + (M2 + M1) * M2 + 3 * M1
    z2 = mat_collect(z, M1)
    assert \
        z2 \
        == (M1 + M2)*M1**2 + (M1 + M2)*M2 + (3*I + M1 + M2)*M1 + M1*M2
    assert z2.expand().doit() == z.expand().doit()
    z2 = mat_collect(z, M1, expand_pow=True)
    assert \
        z2 \
        == (M1 + M2)*M2 + ((M1 + M2)*M1 + 3*I + M1 + M2)*M1 + M1*M2
    assert z2.expand().doit() == z.expand().doit()
    # *
    z2 = mat_collect(z, M2)
    assert \
        z2 \
        == (2*M1 + M2)*M2 + (M1 + M2)*M1**2 + (M1 + M2)*M1 + 3*M1
    assert z2.expand().doit() == z.expand().doit()


def test_partial_apply ():
    z = M1 * M2 + (M2 + M1) * M1 + (M2 + M1) * (M1 ** 2) + (M2 + M1) * M2 + 3 * M1
    assert \
        z.args[1] \
        == (M1 + M2)*M1**2
    assert \
        partial_apply(z, z.args[1], lambda y: y.expand().doit()).doit() \
        == M1**3 + M2*M1**2 + (M1 + M2)*M1 + (M1 + M2)*M2 + 3*M1 + M1*M2
    # *
    assert \
        z.args[1] + z.args[2] \
        == (M1 + M2)*M1**2 + (M1 + M2)*M1
    assert \
        partial_apply(z, z.args[1] + z.args[2], lambda y: y.expand().doit()).doit() \
        == (M1 + M2)*M2 + M1**2 + M1**3 + 3*M1 + M1*M2 + M2*M1**2 + M2*M1
    z2 = mat_collect(z, M2)
    assert \
        partial_apply(z2, z2.args[0], lambda y: y.expand().doit()).doit() \
        == (2*M1 + M2)*M2 + (M1 + M2)*M1**2 + (M1 + M2)*M1 + 3*M1


def test_divide ():
    assert \
        mat_divide(M1 * M2, M1) \
        == M2
    assert \
        mat_divide(M1 * M2, M1, right=True) \
        == M1*M2*M1**(-1)
    assert \
        mat_divide(M1 * M2, M2) \
        == M2**(-1)*M1*M2
    assert \
        mat_divide(M1 * M2, M2, right=True) \
        == M1
    assert \
        mat_divide(M1 * M2 + 2 * x1 * (M1 + M2) * M2, M1) \
        == (2*x1)*M1**(-1)*(M1 + M2)*M2 + M2
    assert \
        mat_divide(M1 * M2 + 2 * x1 * (M1 + M2) * M2, 3 * x1 * (M1 + M2)) \
        == (1/(3*x1))*(M1 + M2)**(-1)*M1*M2 + (Rational(2, 3))*M2
    # *
    assert \
        mat_divide(M1 * M2 + 2 * x1 * (M1 + M2) * M2 + M1 ** 2, M1) \
        == (2*x1)*M1**(-1)*(M1 + M2)*M2 + M1 + M2
    assert \
        mat_divide(M1 * M2 + 2 * x1 * (M1 + M2) * M2 + M1 ** 2, 3 * x1 * (M1 + M2)) \
        == (1/(3*x1))*(M1 + M2)**(-1)*M1**2 + (1/(3*x1))*(M1 + M2)**(-1)*M1*M2 + (Rational(2, 3))*M2
    assert \
        mat_divide(M1 * M2 + 2 * x1 * (M1 + M2) * M2 + M1 ** 2 + M1, M1) \
        == (2*x1)*M1**(-1)*(M1 + M2)*M2 + I + M1 + M2
    assert \
        mat_divide(M1 * M2 + 2 * x1 * (M1 + M2) * M2 + M1 ** 2 + M1 + M1 ** N, M1, right=True) \
        == (2*x1)*(M1 + M2)*M2*M1**(-1) + I + M1**(N - 1) + M1 + M1*M2*M1**(-1)


def test_mixed_1 ():
    z = M1 * M2 + (x1 ** 2) * (M2 + M1) * M1 + x1 * (M2 + M1) * (M1 ** 2) + (M2 + M1) * M2 + 3 * M1
    z2 = mat_collect(z, x1 * M1, expand_pow=True)
    assert \
        z2 \
        == x1*(x1*(M1 + M2) + (M1 + M2)*M1)*M1 + (M1 + M2)*M2 + 3*M1 + M1*M2
    assert z2.expand().doit() == z.expand().doit()
    assert \
      mat_coeff(z2, x1 * M1) \
      == x1*(M1 + M2) + (M1 + M2)*M1
    assert \
      mat_divide(z2, x1 * M1, right=True) \
      == (3/x1)*I + x1*(M1 + M2) + (M1 + M2)*M1 + 1/x1*(M1 + M2)*M2*M1**(-1) + 1/x1*M1*M2*M1**(-1)
    # *
    assert \
        mat_coeff(z2, x1) \
        == (x1*(M1 + M2) + (M1 + M2)*M1)*M1
    assert \
        mat_divide(z2, x1) \
        == (3/x1)*M1 + (x1*(M1 + M2) + (M1 + M2)*M1)*M1 + 1/x1*(M1 + M2)*M2 + 1/x1*M1*M2
    assert \
        mat_collect(z, x1, expand_pow=True) \
        == x1*(x1*(M1 + M2)*M1 + (M1 + M2)*M1**2) + (M1 + M2)*M2 + 3*M1 + M1*M2


def test_mixed_2 ():
    z = x1 ** N * M1**N + x1 * M1 ** 2
    z2 = mat_collect(z, x1 *  M1 ** 2, expand_pow=2)
    assert \
        z2 \
        == x1*((x1*x1**(N - 2))*M1**(N - 2) + I)*M1**2
    assert \
        z2.args[1] \
        == (x1*x1**(N - 2))*M1**(N - 2) + I
    assert \
        partial_apply(z2, z2.args[1], lambda y: y.expand().doit()).doit() \
        == x1*((x1**N/x1)*M1**(N - 2) + I)*M1**2
    assert \
        mat_collect(z, x1 *  M1 ** 2, expand_pow=2, expand_inner=True) \
        == x1*((x1**N/x1)*M1**(N - 2) + I)*M1**2
    assert \
        mat_collect(z, x1 *  M1 ** 2, expand_pow={x1**N: 1, M1**N: 2}) \
        == x1*(x1**(N - 1)*M1**(N - 2) + I)*M1**2

    
def test_mat_trivial_divide ():
    N2 = Symbol("N2", integer=True)
    M3 = MatrixSymbol("M3", N, N2)
    M4 = MatrixSymbol("M4", N2, N)
    
    z = M1 + ((M3 * M4) ** -2) * M3 * M4
    assert \
        mat_trivial_divide(z) \
        == (M3*M4)**(-1) + M1
    z = M1 + ((M3 * M4) ** -1) * M3 * M4
    assert \
        mat_trivial_divide(z) \
        == M1 + Identity(N)
