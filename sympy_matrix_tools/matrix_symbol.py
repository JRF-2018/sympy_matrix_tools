# -*- coding: utf-8 -*-
__version__ = '0.2.1' # Time-stamp: <2022-06-23T12:40:36Z>

from sympy import MatrixSymbol, Wild, Dummy, sympify, \
    MatrixExpr, WildFunction, Symbol
from sympy.core.basic import Basic, Atom
from sympy.core.cache import cacheit
from sympy.core.symbol import Str
from sympy.core.operations import AssocOp
from sympy.matrices.matrices import MatrixKind
from sympy.printing.str import StrPrinter
from sympy.printing.repr import ReprPrinter
from .matrix_function import MatrixFunction, AppliedMatrixUndef


def _new_Basic_dummy_eq (self, other, symbol=None):
    s = self.as_dummy()
    o = _sympify(other)
    o = o.as_dummy()

    dummy_symbols = [i for i in s.free_symbols if i.is_Dummy]

    if len(dummy_symbols) == 1:
        dummy = dummy_symbols.pop()
    else:
        return s == o

    if symbol is None:
        symbols = o.free_symbols

        if len(symbols) == 1:
            symbol = symbols.pop()
        else:
            return s == o

    tmp = dummy.as_dummy()

    return s.xreplace({dummy: tmp}) == o.xreplace({symbol: tmp})

def _new_Basic_as_dummy (self):
    def can(x):
        # mask free that shadow bound
        free = x.free_symbols
        bound = set(x.bound_symbols)
        d = {i: i.as_dummy() for i in bound & free}
        x = x.subs(d)
        # replace bound with canonical names
        x = x.xreplace(x.canonical_variables)
        # return after undoing masking
        return x.xreplace({v: k for k, v in d.items()})
    if not self.has(Symbol):
        return self
    return self.replace(
        lambda x: hasattr(x, 'bound_symbols'),
        can,
        simultaneous=False)

def numbered_name(prefix='x', start=0, exclude=()):
    exclude = set(exclude or [])
    while True:
        name = '%s%s' % (prefix, start)
        if name not in exclude:
            yield name
        start += 1

def _new_Basic_canonical_variables (self):
    if not hasattr(self, 'bound_symbols'):
        return {}
    dums = numbered_name('_')
    reps = {}
    bound = self.bound_symbols
    names = {i.name for i in self.free_symbols - set(bound)}
    for b in bound:
        dn = next(dums)
        if b.is_Symbol:
            while dn in names:
                dn = next(dums)
        if b.is_Matrix:
            d = MatrixSymbol(dn, *b.shape)
        else:
            d = Symbol(dn)
        reps[b] = d
    return reps

def _new_MatrixSymbol_as_dummy (self):
    return MatrixDummy(self.name, *self.shape)

def fix_Basic_MatrixSymbol_about_dummy ():
    Basic.dummy_eq = _new_Basic_dummy_eq
    Basic.as_dummy = _new_Basic_as_dummy
    Basic.canonical_variables = property(_new_Basic_canonical_variables)
    MatrixSymbol.as_dummy = _new_MatrixSymbol_as_dummy


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
    is_commutative = False
    kind: MatrixKind = MatrixKind()

    def __new__ (cls, name=None, n=None, m=None, dummy_index=None, **assumptions):
        n, m = sympify(n), sympify(m)
        MatrixExpr._check_dim(n)
        MatrixExpr._check_dim(m)
        obj = Dummy.__new__(cls, name=str(name), dummy_index=dummy_index, **assumptions)
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
    
    @property
    def args (self):
        return (sympify(self.name), *self.shape,
                sympify(self.dummy_index))

    def matches(self, expr, repl_dict=None, old=False):
        return super(MatrixSymbol, self)\
            .matches(expr, repl_dict=repl_dict, old=old)

    def xreplace (self, rule):
        return super(Atom, self).xreplace(rule)

    def _hashable_content (self):
        return super()._hashable_content() \
            + self.shape

    def as_dummy (self):
        return MatrixDummy(self.name, *self.shape)


class MatrixWild (Wild, MatrixSymbol):
    is_Wild = True
    is_commutative = False
    kind: MatrixKind = MatrixKind()

    def __new__ (cls, name, n, m, exclude=(), properties=(), **assumptions):
        n, m = sympify(n), sympify(m)
        MatrixExpr._check_dim(n)
        MatrixExpr._check_dim(m)
        exclude = tuple([sympify(x) for x in exclude])
        properties = tuple(properties)
        cls._sanitize(assumptions, cls)
        obj = MatrixWild.__xnew__(cls, str(name), n, m, exclude=exclude,
                                  properties=properties, **assumptions)
        return obj
    
    @staticmethod
    @cacheit
    def __xnew__(cls, name, n, m, exclude, properties, **assumptions):
        obj = Symbol.__xnew__(cls, name, **assumptions)
        obj.exclude = exclude
        obj.properties = properties
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
    
    @property
    def args (self):
        return (sympify(self.name), *self.shape)

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
        res = super(MatrixWild, c)._eval_subs(old, new)
        if res is not None:
            return res
        elif modified:
            return c

    def xreplace (self, rule):
        return super(Atom, self).xreplace(rule)

    def _xreplace(self, rule):
        n, m = self.shape
        n_new, n_changed = n._xreplace(rule)
        m_new, m_changed = m._xreplace(rule)
        modified = False
        c = self
        if n_changed or m_changed:
            modified = True
            n = n_new or n
            m = m_new or m
            c = self.func(c.name, n, m,
                          exclude=c.exclude, properties=c.properties,
                          **c.assumptions0)
        r, r_changed = super(MatrixWild, c)._xreplace(rule)
        return r, r_changed or modified

    def matches (self, expr, repl_dict=None, old=False):
        if self.is_Matrix and not expr.is_Matrix:
            return None
        if any(expr.has(x) for x in self.exclude):
            return None
        if not all(f(expr) for f in self.properties):
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

        d[self] = expr
        return d

    def as_dummy (self):
        return MatrixDummy(self.name, *self.shape)


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
        # elif isinstance(expr, type) \
        #      and issubclass(expr, (AppliedMatrixUndef, MatrixFunction)):
        #     pass
        else:
            return None

        try:
            eshape = expr.shape
        except:
            return None

        if repl_dict is None:
            repl_dict = dict()
        else:
            repl_dict = repl_dict.copy()

        d = repl_dict
        for arg, other_arg in zip(self.shape, eshape):
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

        d[self] = expr
        return d


def _ReprPrinter__print_WildMatrixFunction (self, expr):
    return self._print_WildFunction(expr)

ReprPrinter._print_WildMatrixFunction = _ReprPrinter__print_WildMatrixFunction

def _StrPrinter__print_WildMatrixFunction (self, expr):
    return self._print_WildFunction(expr)

StrPrinter._print_WildMatrixFunction = _StrPrinter__print_WildMatrixFunction
