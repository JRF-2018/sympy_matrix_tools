# -*- coding: utf-8 -*-
__version__ = '0.7.0' # Time-stamp: <2026-10-22T17:31:49Z>

from sympy import symbols, MatrixSymbol, Eq, ZeroMatrix, Symbol
from sympy_matrix_tools.derivation import Derivation, Grad, cas, rule

n, m = symbols("n m", integer=True, positive=True)
A = MatrixSymbol("A", n, n)
B = MatrixSymbol("B", n, m)
W = MatrixSymbol("W", n, n)
y0 = MatrixSymbol("y0", n, 1)
z = MatrixSymbol("z", n, 1)
p = MatrixSymbol("p", n, 1)
x = MatrixSymbol("x", m, 1)
K = MatrixSymbol("K", m, m)
yx = A*y0 + B*x + z
J = yx.T*W*yx/2 + p.T*yx


def test_wrong_steps_are_found ():
    d = Derivation("wrong steps")
    d.hypothesis(Eq(W.T, W), "W_sym")
    # Asserted without justification, and the factor 1/2 is wrong.
    d.calc("grad_wrong", Grad(J, x),
           B.T*W*yx/2 + B.T*p)
    # A wrong claim that it is the result of doit.
    d.calc("grad_wrong_cas", Grad(J, x),
           (B.T*W.T*yx + B.T*p/2, cas("doit")))
    assert not d.check(trials=2)
    st1 = d.facts["grad_wrong"].steps[0]
    st2 = d.facts["grad_wrong_cas"].steps[0]
    assert st1.symbolic is None
    assert st1.numeric == (0, 2)
    assert st1.counterexample is not None
    assert st2.symbolic is False
    assert st2.numeric == (0, 2)
    s = d.summary()
    assert s["asserted"] == 1
    assert s["symbolic_failures"] == 1
    assert s["numeric_failures"] == 2


def test_open_obligation ():
    d = Derivation("open obligation")
    d.define(K, B.T*W*B, "K_def")
    d.hypothesis(Eq(K*x + B.T*p, ZeroMatrix(m, 1)), "eq")
    d.solve_linear("x_sol", "eq", x)
    assert d.check(trials=2)
    assert d.summary()["open_obligations"] == 1
    assert "obligation Q.invertible(K): OPEN" in d.report()


def test_scalar_steps ():
    t = Symbol("t")
    a = Symbol("a")
    d = Derivation("scalars")
    d.calc("square", (t + a)**2,
           (t**2 + 2*a*t + a**2, cas("expand")),
           t**2 + 2*a*t + a**2 + 1)
    assert not d.check(trials=2)
    steps = d.facts["square"].steps
    assert steps[0].symbolic is True and steps[0].numeric == (2, 2)
    assert steps[1].symbolic is None and steps[1].numeric == (0, 2)


def test_matrix_functions ():
    from sympy import Function, Symbol
    from sympy_matrix_tools import MatrixFunction
    t = Symbol("t")
    X = MatrixSymbol("X", n, n)
    xf = MatrixFunction("x", n, 1)
    F = MatrixFunction("F", n, n)
    G = MatrixFunction("G", n, n)
    f = Function("f")
    d = Derivation("functions")
    d.hypothesis(Eq(W.T, W), "W_sym")
    d.calc("dq", (xf(t).T*W*xf(t)).diff(t),
           xf(t).T*W*xf(t).diff(t) + xf(t).diff(t).T*W*xf(t),
           (2*xf(t).diff(t).T*W*xf(t), rule("W_sym")),
           # Wrong: the factor 2 is missing.
           xf(t).diff(t).T*W*xf(t))
    d.calc("dF", F(t**2).diff(t),
           2*t*F(t**2).diff(t).args[1],
           # Wrong: the factor t is missing.
           2*F(t**2).diff(t).args[1])
    d.calc("dG", G(t*X).diff(t), G(t*X).diff(t))
    d.calc("ff", f(t)**2 + f(t).diff(t),
           (f(t)*f(t) + f(t).diff(t), cas("expand")))
    assert not d.check(trials=2)
    dq = d.facts["dq"].steps
    assert [st.numeric for st in dq] == [(2, 2), (2, 2), (0, 2)]
    dF = d.facts["dF"].steps
    assert [st.numeric for st in dF] == [(2, 2), (0, 2)]
    assert d.facts["dG"].steps[0].numeric == (2, 2)
    assert d.facts["ff"].steps[0].numeric == (2, 2)
    r = d.report()
    assert "x'(t)" in r
    assert "F'(t**2)" in r
    assert "DG(t*X)[X]" in r

