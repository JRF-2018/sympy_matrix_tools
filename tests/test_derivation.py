# -*- coding: utf-8 -*-
__version__ = '0.6.0' # Time-stamp: <2026-10-22T16:54:53Z>

import os
import runpy
from sympy import symbols, MatrixSymbol, Eq, ZeroMatrix, Symbol
from sympy_matrix_tools.derivation import Derivation, Grad, cas

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


def test_lq_control_example ():
    path = os.path.join(os.path.dirname(__file__), "..", "examples",
                        "lq_control.py")
    d = runpy.run_path(path)["d"]
    assert d.check(trials=3)
    s = d.summary()
    assert s["symbolic_failures"] == 0
    assert s["numeric_failures"] == 0
    assert s["open_obligations"] == 0
    assert s["asserted"] == 0
    assert "x = -K**(-1)*B.T*(W*(A*y0 + z) + p)" in d.report()
    assert d._repr_latex_().startswith(r"$\begin{array}")


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
