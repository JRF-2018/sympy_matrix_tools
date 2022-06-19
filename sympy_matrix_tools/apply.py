# -*- coding: utf-8 -*-
__version__ = '0.2.0' # Time-stamp: <2022-06-19T07:31:41Z>

from sympy import Expr, MatrixExpr, AppliedPredicate, Lambda, sympify
from sympy.logic.boolalg import Boolean
from sympy.printing.str import StrPrinter
from sympy.printing.repr import ReprPrinter
from sympy.printing.latex import LatexPrinter, latex_escape

class Apply (Expr):
    def __new__ (cls, lmd, *args, **kwargs):
        lmd = sympify(lmd)
        if isinstance(lmd, Lambda) or callable(lmd):
            return lmd(*args, **kwargs)
        return Expr.__new__(cls, lmd, *args, **kwargs)

    @property
    def function (self):
        return self.args[0]

    @property
    def arguments (self):
        return self.args[1:]

    def doit (self, **kwargs):
        deep = kwargs.get('deep', True)
        if deep:
            args = [arg.doit(**kwargs) for arg in self.args]
        else:
            args = self.args
        return self.func(*self.args)


class MatApply(Apply, MatrixExpr):
    is_Matrix = True
    
    def __new__ (cls, m, n, lmd, *args, **kwargs):
        lmd = sympify(lmd)
        if isinstance(lmd, Lambda) or callable(lmd):
            return lmd(*args, **kwargs)
        return MatrixExpr.__new__(cls, m, n, lmd, *args, **kwargs)

    @property
    def shape (self):
        return (self.args[0], self.args[1])

    @property
    def function (self):
        return self.args[2]

    @property
    def arguments (self):
        return self.args[3:]

    def could_extract_minus_sign (self):
        return False


class PredApply(Apply, AppliedPredicate):
    def __new__ (cls, lmd, *args, **kwargs):
        lmd = sympify(lmd)
        if isinstance(lmd, Lambda) or callable(lmd):
            return lmd(*args, **kwargs)
        args = map(sympify, args)
        return Boolean.__new__(cls, lmd, *args, **kwargs)
    
    @property
    def function (self):
        return self.args[0]

    @property
    def arguments (self):
        return self.args[1:]


def _LatexPrinter__print_PredApply (self, expr):
    pred = expr.__class__.__name__
    args = expr.args
    pred_latex = r"\operatorname{{Q}}_{{\text{{{}}}}}"\
        .format(latex_escape(pred))
    args_latex = ', '.join([self._print(a) for a in args])
    return '%s(%s)' % (pred_latex, args_latex)

LatexPrinter._print_PredApply = _LatexPrinter__print_PredApply

def _ReprPrinter__print_PredApply (self, expr):
    args = expr._args
    return "%s(%s)" % (expr.__class__.__name__, self.reprify(args, ", "))

ReprPrinter._print_PredApply = _ReprPrinter__print_PredApply

def _StrPrinter__print_PredApply (self, expr):
    return '%s(%s)' % (
        expr.__class__.__name__, self.stringify(expr.args, ", "))

StrPrinter._print_PredApply = _StrPrinter__print_PredApply
