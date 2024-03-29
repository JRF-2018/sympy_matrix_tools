# -*- coding: utf-8 -*-
__version__ = '0.3.6' # Time-stamp: <2022-10-12T10:21:46Z>

from sympy import Function, MatrixExpr, sympify, MatrixSymbol, Expr
from sympy.matrices.expressions.matexpr import MatrixElement, _LeftRightArgs
from sympy.core.function import UndefinedFunction, Application, AppliedUndef
from packaging.version import parse as parse_version

# class MatrixFunction0 (Function, MatrixExpr):
#     def __new__ (cls, *args, **options):
#         if cls is MatrixFunction0:
#             return UndefinedMatrixFunction0(*args, **options)
#         args = list(args)
#         if len(args) < 2:
#             raise TypeError("MatrixFunction needs rows and cols as arguments.")
#         n, m = args.pop(0), args.pop(0)

#         n, m = sympify(n), sympify(m)

#         cls._check_dim(m)
#         cls._check_dim(n)

#         obj = super().__new__(cls, n, m, *args, **options)
#         return obj

#     @property
#     def shape (self):
#         return self.args[0], self.args[1]

#     def _entry (self, i, j, **kwargs):
#         return MatrixElement(self, i, j)


# class AppliedMatrixUndef0 (MatrixFunction0):
#     is_number = False

#     def __new__ (cls, *args, **options):
#         args = list(map(sympify, args))
#         u = [a.name for a in args if isinstance(a, UndefinedFunction)
#              or isinstance(a, UndefinedMatrixFunction0)]
#         if u:
#             raise TypeError('Invalid argument: expecting an expression, not UndefinedFunction%s: %s' % (
#                 's'*(len(u) > 1), ', '.join(u)))
#         obj = super().__new__(cls, *args, **options)
#         return obj

#     def _eval_as_leading_term (self, x, logx=None, cdir=0):
#         return self

#     @property
#     def _diff_wrt (self):
#         return True


# class UndefinedMatrixFunction0 (UndefinedFunction):
#     def __new__ (mcl, name, bases=(AppliedMatrixUndef0,), __dict__=None, **kwargs):
#         obj = super().__new__(mcl, name, bases, __dict__, **kwargs)
#         return obj


class MatrixFunction (Function, MatrixExpr):
    def __new__ (cls, *args, **options):
        if cls is MatrixFunction:
            return UndefinedMatrixFunction(*args, **options)
        obj = super().__new__(cls, *args, **options)
        return obj

    @property
    def free_symbols (self):
        return super().free_symbols \
            | self.shape[0].free_symbols \
            | self.shape[1].free_symbols

    def _eval_subs (self, old, new):
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
        res = super(MatrixFunction, c)._eval_subs(old, new)
        if res is not None:
            return res
        elif modified:
            return c

    def _xreplace (self, rule):
        n, m = self.shape
        n_new, n_changed = n._xreplace(rule)
        m_new, m_changed = m._xreplace(rule)
        modified = False
        if n_changed or m_changed:
            modified = True
            n = n_new or n
            m = m_new or m
            u = UndefinedMatrixFunction(self.name, m_new, n_new)
            c = u(*[i._xreplace(rule)[0] for i in self.args])
        else:
            c = self
        r, r_changed = super(MatrixFunction, c)._xreplace(rule)
        return r, r_changed or modified

    def matches (self, expr, repl_dict=None, old=False):
        if self.is_Matrix and not expr.is_Matrix:
            return None
        
        if repl_dict is None:
            repl_dict = dict()
        else:
            repl_dict = repl_dict.copy()

        d = repl_dict
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

        return super(MatrixFunction, self.xreplace(d))\
            .matches(expr, d, old=old)


class AppliedMatrixUndef (MatrixFunction, AppliedUndef):
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


class UndefinedMatrixFunction (UndefinedFunction):
    def __init__ (cls, name, n, m, **kwargs):
        super().__init__(cls, name, **kwargs)
        
    def __new__ (mcl, name, n, m, bases=(AppliedMatrixUndef,), __dict__=None, **kwargs):
        n, m = sympify(n), sympify(m)

        MatrixExpr._check_dim(n)
        MatrixExpr._check_dim(m)

        obj = super().__new__(mcl, name, bases, __dict__, **kwargs)
        obj.shape = (n, m)
        return obj

    _kwargs = {}  # type: tDict[str, Optional[bool]]

    def __hash__ (self):
        return hash((self.class_key(), frozenset(self._kwargs.items()),
                     ('shape', self.shape)))

    def __eq__ (self, other):
        return super().__eq__(other) and \
            self.shape == other.shape

    def __ne__ (self, other):
        return not self == other


try:
    from sympy.tensor.array.expressions.arrayexpr_derivatives \
        import array_derive

    @array_derive.register(MatrixFunction)
    def _ (expr: MatrixFunction, x: Expr):
        raise NotImplementedError()

    # from sympy import Subs
    # from sympy.core.function import ArgumentIndexError
    # from sympy.tensor.array.array_derivatives \
    #     import ArrayDerivative
    # from sympy.tensor.array.expressions.array_expressions import (
    #     get_rank, _array_tensor_product, _array_diagonal, _array_add,
    #     _array_contraction, _permute_dims, _ArrayExpr)
    # from .matrix_symbol import make_dummy

    # class ArraySubs (Subs, _ArrayExpr):
    #     shape = property(lambda self: self.args[0].shape)
    #     # @property
    #     # def shape (self):
    #     #     return (self.args[0].shape[-2], self.args[0].shape[-1])

    # # xxx = [(0, 4), (1, 3)]
    # xxx = [(2, 4), (1, 3)]

    # @array_derive.register(AppliedMatrixUndef)
    # def _ (expr: AppliedMatrixUndef, x: Expr):
    #     def fdiff (self, argindex=1):
    #         if not (1 <= argindex <= len(self.args)):
    #             raise ArgumentIndexError(self, argindex)
    #         ix = argindex - 1
    #         A = self.args[ix]
    #         D = make_dummy(A, 'xi_%i' % argindex, dummy_index=hash(A))
    #         args = self.args[:ix] + (D,) + self.args[ix + 1:]
    #         return ArraySubs(ArrayDerivative(self.func(*args), D,
    #                                          evaluate=False), D, A)

    #     # I don't know correct algorithm and am trying plausible one.
    #     i = 0
    #     l = []
    #     b = get_rank(x)
    #     c = get_rank(expr)
    #     for a in expr.args:
    #         i += 1
    #         da = array_derive(a, x)
    #         if da.is_zero:
    #             continue
    #         try:
    #             df = fdiff(expr, i)
    #         except ArgumentIndexError:
    #             df = Function.fdiff(self, i)
    #         # da = _array_diagonal(da, *[(i, b + i) for i in range(b)])
    #         df = _array_diagonal(df, *[(i, c + i) for i in range(c)])
    #         tp = _array_tensor_product(da, df)
    #         # tp = _array_contraction(tp, *xxx)
    #         # tp = _array_diagonal(tp, *[(b + i, b + c + i) for i in range(c)])
    #         tp = _array_diagonal(tp, *[(i, b + i) for i in range(b)])
    #         l.append(tp)
    #     # print(tp.shape, get_rank(tp))
    #     tp = _array_add(*l)
    #     # diag_indices = [(i, b + i, b + c + i) for i in range(b + c)]
    #     # tp = _array_diagonal(tp, *diag_indices)
    #     # print(tp.shape, get_rank(tp))
    #     return tp

except ImportError:
    pass


def freeze_matrix_function (z, exclude=None):
    if exclude is None:
        exclude = set()
    l = []
    s = sorted(list(z.atoms(MatrixFunction)), key=lambda x: len(str(x)),
               reverse=True)
    for x in s:
        if not (x.func in exclude or x in exclude):
            l.append((x, MatrixSymbol(str(x), x.shape[0], x.shape[1])))
    d = dict(l)
    return z.subs(d), d


def melt_matrix_function (z, sigma):
    d = dict([(v, n) for n, v in sigma.items()])
    return z.subs(d)
