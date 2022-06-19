# -*- coding: utf-8 -*-
__version__ = '0.2.0' # Time-stamp: <2022-06-19T07:33:53Z>

from sympy import MatrixSymbol, Wild, Dummy, sympify, MatrixExpr, WildFunction
from sympy.core.symbol import Str
from sympy.core.operations import AssocOp
from sympy.matrices.matrices import MatrixKind
from sympy.printing.str import StrPrinter
from sympy.printing.repr import ReprPrinter
from .matrix_function import MatrixFunction, AppliedMatrixUndef


## fix_AssocOp__matches_commutative was needed for MatrixWild.  However, it seems not to be needed now.

_orig_AssocOp__matches_commutative = None
def _new_AssocOp__matches_commutative (self, expr, repl_dict=None, old=False):
    try:
        return _orig_AssocOp__matches_commutative(self, expr, repl_dict, old)
    except TypeError:
        return None

def fix_AssocOp__matches_commutative ():
    global _orig_AssocOp__matches_commutative
    if _orig_AssocOp__matches_commutative is None:
        _orig_AssocOp__matches_commutative = AssocOp._matches_commutative
    AssocOp._matches_commutative = _new_AssocOp__matches_commutative


class MatrixDummy (Dummy, MatrixSymbol):
    is_Dummy = True
    kind: MatrixKind = MatrixKind()

    def __new__ (cls, name=None, n=None, m=None, dummy_index=None, **assumptions):
        n, m = sympify(n), sympify(m)
        MatrixExpr._check_dim(n)
        MatrixExpr._check_dim(m)
        obj = Dummy.__new__(cls, name=name, dummy_index=dummy_index, **assumptions)
        obj._shape = (n, m)
        return obj

    @property
    def shape (self):
        return self._shape

    @property
    def name (self):
        return self._name

    @name.setter
    def name (self, c):
        self._name = c
    
    def _hashable_content (self):
        return super()._hashable_content() \
            + self.shape

    def _eval_subs (self, old, new):
        n, m = self.shape
        n_new = n._subs(old, new)
        m_new = m._subs(old, new)
        modified = False
        c = self
        if n != n_new or m != m_new:
            modified = True
            n = n_new or n
            m = m_new or m
            c = self.func(name=c.name, n=n, m=m, dummy_index=c.dummy_index,
                          **c.assumptions0)
        res = super(Dummy, c)._eval_subs(old, new)
        if res is not None:
            return res
        elif modified:
            return c


class MatrixWild (Wild, MatrixSymbol):
    is_Wild = True
    kind: MatrixKind = MatrixKind()

    def __new__ (cls, name, n, m, exclude=(), properties=(), **assumptions):
        n, m = sympify(n), sympify(m)
        MatrixExpr._check_dim(n)
        MatrixExpr._check_dim(m)
        obj = Wild.__new__(cls, name, exclude=exclude, properties=properties, **assumptions)
        obj._shape = (n, m)
        # obj._args = (Str(obj.name), n, m)
        return obj

    # @property
    # def args (self):
    #     return self._args

    @property
    def shape (self):
        return self._shape

    @property
    def name (self):
        return self._name

    @name.setter
    def name (self, c):
        self._name = c
    
    def _hashable_content (self):
        return super()._hashable_content() \
            + self.shape

    def _eval_subs (self, old, new):
        n, m = self.shape
        n_new = n._subs(old, new)
        m_new = m._subs(old, new)
        modified = False
        c = self
        if n != n_new or m != m_new:
            modified = True
            n = n_new or n
            m = m_new or m
            c = self.func(c.name, n, m,
                          exclude=c.exclude, properties=c.properties,
                          **c.assumptions0)
        res = super(Wild, c)._eval_subs(old, new)
        if res is not None:
            return res
        elif modified:
            return c

    def matches (self, expr, repl_dict=None, old=False):
        if self.is_Matrix and not expr.is_Matrix:
            return None
        if self.shape != expr.shape:
            return None
        if any(expr.has(x) for x in self.exclude):
            return None
        if not all(f(expr) for f in self.properties):
            return None
        if repl_dict is None:
            repl_dict = dict()
        else:
            repl_dict = repl_dict.copy()
        repl_dict[self] = expr
        return repl_dict


class WildMatrixFunction(WildFunction, MatrixFunction):  # type: ignore
    def __init__ (cls, name, n, m, **assumptions):
        WildFunction.__init__(cls, name, **assumptions)
        MatrixFunction.__init__(cls)

    @property
    def shape (self):
        return (sympify(self.args[1]), sympify(self.args[2]))

    def matches (self, expr, repl_dict=None, old=False):
        if isinstance(expr, (AppliedMatrixUndef, MatrixFunction)):
            if len(expr.args) not in self.nargs:
                return None
        elif isinstance(expr, type) \
             and issubclass(expr, (AppliedMatrixUndef, MatrixFunction)):
            pass
        else:
            return None

        try:
            eshape = expr.shape
        except:
            return None

        if self.shape != eshape:
            return None

        if repl_dict is None:
            repl_dict = dict()
        else:
            repl_dict = repl_dict.copy()

        repl_dict[self] = expr
        return repl_dict


def _ReprPrinter__print_WildMatrixFunction (self, expr):
    return self._print_WildFunction(expr)

ReprPrinter._print_WildMatrixFunction = _ReprPrinter__print_WildMatrixFunction

def _StrPrinter__print_WildMatrixFunction (self, expr):
    return self._print_WildFunction(expr)

StrPrinter._print_WildMatrixFunction = _StrPrinter__print_WildMatrixFunction
