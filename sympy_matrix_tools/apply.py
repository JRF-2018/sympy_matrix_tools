# -*- coding: utf-8 -*-
__version__ = '0.2.1' # Time-stamp: <2022-06-22T17:20:42Z>

from sympy import Expr, MatrixExpr, AppliedPredicate, Lambda,\
    sympify, Wild, Dummy
from sympy.logic.boolalg import Boolean
from sympy.printing.str import StrPrinter
from sympy.printing.repr import ReprPrinter
from sympy.printing.latex import LatexPrinter, latex_escape
from .matrix_symbol import MatrixDummy


def make_dummy (expr):
    if expr.is_Matrix:
        return MatrixDummy("X", *expr.shape)
    return Dummy("x")


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


    def matches (self, expr, repl_dict=None, old=False):
        d = super().matches(expr, repl_dict=repl_dict, old=old)
        if d is not None:
            return d
        if not isinstance(self.function, Wild):
            return None

        if repl_dict is None:
            repl_dict = dict()
        else:
            repl_dict = repl_dict.copy()

        d = repl_dict
        
        s = {}
        xs = []
        for a in self.arguments:
            x = make_dummy(a)
            xs.append(x)
            if a not in s:
                s[a] = x

        d[self.function] = Lambda(tuple(xs), expr.subs(s))
        return d


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

    def matches (self, expr, repl_dict=None, old=False):
        d = super().matches(expr, repl_dict=repl_dict, old=old)
        if d is None:
            return None

        for arg, other_arg in zip(self.shape, expr.shape):
            if arg == other_arg:
                continue
            if arg.is_Relational:
                try:
                    d = arg.xreplace(d).matches(other_arg, d, old=old)
                except TypeError:
                    d = None
            else:
                    d = arg.xreplace(d).matches(other_arg, d, old=old)
            if d is None:
                return None
        return d


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
