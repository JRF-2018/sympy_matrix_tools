# -*- coding: utf-8 -*-
__version__ = '0.1.1' # Time-stamp: <2022-05-04T12:28:51Z>
## Language: Japanese/UTF-8

import sympy
from sympy import Mul, Add, Sum, Product, MatrixExpr, StrPrinter, \
    trace, transpose, adjoint
from sympy.matrices.expressions.matexpr import MatrixElement


def _StrPrinter_print_MatSum (self, expr):
    def _xab_tostr(xab):
        if len(xab) == 1:
            return self._print(xab[0])
        else:
            return self._print((xab[0],) + tuple(xab[1:]))
    L = ', '.join([_xab_tostr(l) for l in expr.limits])
    return 'MatSum(%s, %s)' % (self._print(expr.function), L)
StrPrinter._print_MatSum = _StrPrinter_print_MatSum


def _StrPrinter_print_MatProduct (self, expr):
    def _xab_tostr(xab):
        if len(xab) == 1:
            return self._print(xab[0])
        else:
            return self._print((xab[0],) + tuple(xab[1:]))
    L = ', '.join([_xab_tostr(l) for l in expr.limits])
    return 'MatSum(%s, %s)' % (self._print(expr.function), L)
StrPrinter._print_MatSum = _StrPrinter_print_MatSum


class MatSum (Sum, MatrixExpr):
    @property
    def is_Matrix (self):
        return self.args[0].is_Matrix

    @property
    def shape (self):
        return self.args[0].shape

    def could_extract_minus_sign(self):
        return self.args[0].could_extract_minus_sign()

    def _entry(self, i, j, **kwargs):
        return Sum(self.args[0]._entry(i, j, **kwargs))

    def _eval_transpose(self):
        return self.func(*([transpose(self.args[0])] + self.args[1:])).doit()

    def _eval_adjoint(self):
        return self.func(*([adjoint(self.args[0])] + self.args[1:])).doit()

    def _eval_trace(self):
        return Sum(*([trace(self.args[0])] + self.args[1:])).doit()

    def doit(self, **kwargs):
        deep = kwargs.get('deep', True)
        if deep:
            args = [arg.doit(**kwargs) for arg in self.args]
        else:
            args = self.args
        return self.func(*args)


class MatProduct (Product, MatrixExpr):
    @property
    def is_Matrix (self):
        return self.args[0].is_Matrix

    @property
    def shape (self):
        return self.args[0].shape

    def _entry(self, i, j, **kwargs):
        return MatrixElement(self, i, j)

    def could_extract_minus_sign(self):
        x = self.args[0].could_extract_minus_sign()
        if not x:
            return False
        return Product(x, self.args[1:]).could_extract_minus_sign()

    def doit(self, **kwargs):
        deep = kwargs.get('deep', True)
        if deep:
            args = [arg.doit(**kwargs) for arg in self.args]
        else:
            args = self.args
        return self.func(*args)


