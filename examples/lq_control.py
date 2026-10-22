# -*- coding: utf-8 -*-
__version__ = '0.6.0' # Time-stamp: <2026-10-22T16:54:53Z>

"""The optimal control of a one-period linear-quadratic problem,
checked by a Derivation ledger.

The state y is given by y = A*y0 + B*x + z, where x is the control.
The objective J = 1/2 * y.T*Q*y + p.T*y is minimized with respect to
x.  We derive the gradient, the first-order condition and the optimal
control x = -K**-1 * B.T*(Q*(A*y0 + z) + p) with K = B.T*Q*B.

This is a simplified version of a derivation in
https://github.com/JRF-2018/economy_control (economy_control_03_2b).
"""

from sympy import symbols, MatrixSymbol, Eq, ZeroMatrix, Q
from sympy_matrix_tools.derivation import Derivation, Grad, cas, rule

n, m = symbols("n m", integer=True, positive=True)
A = MatrixSymbol("A", n, n)
B = MatrixSymbol("B", n, m)
W = MatrixSymbol("W", n, n)     # the weight matrix (Q in the book)
y0 = MatrixSymbol("y0", n, 1)
z = MatrixSymbol("z", n, 1)
p = MatrixSymbol("p", n, 1)
x = MatrixSymbol("x", m, 1)
y = MatrixSymbol("y", n, 1)
K = MatrixSymbol("K", m, m)

d = Derivation("One-period LQ control")
d.hypothesis(Eq(W.T, W), "W_sym")
d.hypothesis(Q.invertible(K), "K_inv")
d.define(y, A*y0 + B*x + z, "state")
d.define(K, B.T*W*B, "K_def")

# Grad needs the objective written explicitly in x.
yx = A*y0 + B*x + z
J = yx.T*W*yx/2 + p.T*yx
d.alias(J, "J")
d.calc("J_by_y", J,
       (y.T*W*y/2 + p.T*y, rule("state")))
d.calc("grad", Grad(J, x),
       (B.T*W*yx/2 + B.T*W.T*yx/2 + B.T*p, cas("doit")),
       (B.T*W*yx + B.T*p, rule("W_sym")),
       (K*x + B.T*(W*(A*y0 + z) + p), rule("K_def")))

d.hypothesis(Eq(Grad(J, x), ZeroMatrix(m, 1), evaluate=False), "foc")
d.calc("foc_linear", K*x + B.T*(W*(A*y0 + z) + p),
       (Grad(J, x), rule("grad")),
       (ZeroMatrix(m, 1), rule("foc")))
d.solve_linear("x_opt", "foc_linear", x)

if __name__ == '__main__':
    d.check(trials=3)
    print(d.report())
