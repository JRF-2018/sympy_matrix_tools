# -*- coding: utf-8 -*-
__version__ = '0.2.0' # Time-stamp: <2022-06-19T01:28:07Z>

import pytest
from sympy import MatrixSymbol, Symbol, Function, Predicate, Lambda, And, latex
from sympy_matrix_tools import *

Q = Symbol("Q")
MQ = MatrixSymbol("MQ", 3, 3)
f = Function("f")
Mf = MatrixFunction("Mf", 3, 3)
b = Predicate("b")

def test_apply ():
    assert \
        str(Apply(Q, 3) + 3) \
        == "3 + Apply(Q, 3)"
    assert \
        str(MatApply(3, 3, Q, 3) + MQ) \
        == "MatApply(3, 3, Q, 3) + MQ"
    assert \
        str(And(PredApply(Q, 3), b(1))) \
        == "PredApply(Q, 3) & Q.b(1)"
    assert \
        latex(PredApply(Q, 3)) \
        == r"\operatorname{Q}_{\text{PredApply}}(Q, 3)"
    assert \
        repr(PredApply(Q, 3)) \
        == "PredApply(Q, 3)"
    assert \
        Apply(Q, 3).subs(Q, Lambda(Q, Q + 1)) \
        == 4
    assert \
        Apply(Q, 3).subs(Q, f) \
        == f(3)
    assert \
        MatApply(3, 3, Q, 3).subs(Q, Lambda(Q, Mf(Q))) \
        == Mf(3)
    assert \
        MatApply(3, 3, Q, 3).subs(Q, Mf) \
        == Mf(3)
    assert \
        PredApply(Q, 3).subs(Q, Lambda(Q, b(Q))) \
        == b(3)
    assert \
        PredApply(Q, 3).subs(Q, b) \
        == b(3)
