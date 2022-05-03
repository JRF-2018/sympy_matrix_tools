# -*- coding: utf-8 -*-
__version__ = '0.1.0' # Time-stamp: <2022-05-03T05:13:27Z>
## Language: Japanese/UTF-8

from sympy import Function, MatrixExpr, sympify
from sympy.matrices.expressions.matexpr import MatrixElement, _LeftRightArgs
from sympy.core.function import UndefinedFunction, Application


class MatrixFunction0 (Function, MatrixExpr):
    def __new__ (cls, *args, **options):
        if cls is MatrixFunction0:
            return UndefinedMatrixFunction0(*args, **options)
        args = list(args)
        if len(args) < 2:
            raise TypeError("MatrixFunction needs rows and cols as arguments.")
        n, m = args.pop(0), args.pop(0)

        n, m = sympify(n), sympify(m)

        cls._check_dim(m)
        cls._check_dim(n)

        obj = super().__new__(cls, n, m, *args, **options)
        return obj

    @property
    def shape (self):
        return self.args[0], self.args[1]

    def _entry(self, i, j, **kwargs):
        return MatrixElement(self, i, j)


class AppliedMatrixUndef0 (MatrixFunction0):
    is_number = False

    def __new__(cls, *args, **options):
        args = list(map(sympify, args))
        u = [a.name for a in args if isinstance(a, UndefinedFunction)
             or isinstance(a, UndefinedMatrixFunction0)]
        if u:
            raise TypeError('Invalid argument: expecting an expression, not UndefinedFunction%s: %s' % (
                's'*(len(u) > 1), ', '.join(u)))
        obj = super().__new__(cls, *args, **options)
        return obj

    def _eval_as_leading_term(self, x, logx=None, cdir=0):
        return self

    @property
    def _diff_wrt(self):
        return True


class UndefinedMatrixFunction0(UndefinedFunction):
    def __new__(mcl, name, bases=(AppliedMatrixUndef0,), __dict__=None, **kwargs):
        obj = super().__new__(mcl, name, bases, __dict__, **kwargs)
        return obj


class MatrixFunction (Function, MatrixExpr):
    def __new__ (cls, *args, **options):
        if cls is MatrixFunction:
            return UndefinedMatrixFunction(*args, **options)
        obj = super().__new__(cls, *args, **options)
        return obj


class AppliedMatrixUndef (MatrixFunction):
    is_number = False

    def __new__ (cls, *args, **options):
        args = list(map(sympify, args))
        u = [a.name for a in args if isinstance(a, UndefinedFunction)
             or isinstance(a, UndefinedMatrixFunction)]
        if u:
            raise TypeError('Invalid argument: expecting an expression, not UndefinedFunction%s: %s' % (
                's'*(len(u) > 1), ', '.join(u)))
        obj = super().__new__(cls, *args, **options)
        obj.shape = (cls.shape[0], cls.shape[1])
        return obj

    def _eval_as_leading_term (self, x, logx=None, cdir=0):
        return self

    @property
    def _diff_wrt (self):
        return True

    def _entry (self, i, j, **kwargs):
        return MatrixElement(self, i, j)

    @property
    def free_symbols (self):
        return super().free_symbols \
            | self.shape[0].free_symbols \
            | self.shape[1].free_symbols

    def _eval_subs(self, old, new):
        n, m = self.shape
        n_new = n._subs(old, new)
        m_new = m._subs(old, new)
        modified = False
        if n != n_new or m != m_new:
            modified = True
            n = n_new or n
            m = m_new or m
            u = UndefinedMatrixFunction(self.name, m_new, n_new)
            c = u(*[i._subs(old, new) for i in self.args])
        else:
            c = self
        if old == self.func:
            old = c.func
        res = super(AppliedMatrixUndef, c)._eval_subs(old, new)
        if res is not None:
            return res
        elif modified:
            return c


class UndefinedMatrixFunction(UndefinedFunction):
    def __new__(mcl, name, n, m, bases=(AppliedMatrixUndef,), __dict__=None, **kwargs):
        n, m = sympify(n), sympify(m)

        MatrixExpr._check_dim(n)
        MatrixExpr._check_dim(m)

        obj = super().__new__(mcl, name, bases, __dict__, **kwargs)
        obj.shape = (n, m)
        return obj

    _kwargs = {}  # type: tDict[str, Optional[bool]]

    def __hash__(self):
        return hash((self.class_key(), frozenset(self._kwargs.items())))

    def __eq__(self, other):
        return super().__eq__(other) and \
            self.shape == other.shape

    def __ne__(self, other):
        return not self == other
