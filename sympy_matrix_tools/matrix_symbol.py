# -*- coding: utf-8 -*-
__version__ = '0.2.3' # Time-stamp: <2022-06-24T23:31:21Z>

from sympy import MatrixSymbol, Wild, Dummy, sympify, \
    MatrixExpr, WildFunction, Symbol
from sympy.core.basic import Basic, Atom
from sympy.core.operations import AssocOp
from sympy.core.assumptions import StdFactKB
from sympy.matrices.matrices import MatrixKind
from sympy.printing.str import StrPrinter
from sympy.printing.repr import ReprPrinter
from .matrix_function import MatrixFunction, AppliedMatrixUndef


def _new_Basic_dummy_eq (self, other, symbol=None):
    s = self.as_dummy()
    o = sympify(other) # changed
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

    tmp = dummy.as_dummy() # changed

    return s.xreplace({dummy: tmp}) == o.xreplace({symbol: tmp})

def _new_Basic_as_dummy (self):
    def can(x):
        # mask free that shadow bound
        free = x.free_symbols
        bound = set(x.bound_symbols)
        d = {i: i.as_dummy() for i in bound & free} # changed
        x = x.subs(d)
        # replace bound with canonical names
        x = x.xreplace(x.canonical_variables)
        # return after undoing masking
        return x.xreplace({v: k for k, v in d.items()})
    if not self.has(Symbol) and not self.has(MatrixSymbol): # changed
        return self
    return self.replace(
        lambda x: hasattr(x, 'bound_symbols'),
        can,
        simultaneous=False)

def numbered_name (prefix='x', start=0, exclude=()):
    exclude = set(exclude or [])
    while True:
        name = '%s%s' % (prefix, start)
        if name not in exclude:
            yield name
        start += 1

def _new_Basic_canonical_variables (self):
    if not hasattr(self, 'bound_symbols'):
        return {}
    bound = self.bound_symbols
    names = {i.name for i in self.free_symbols - set(bound)}
    dums = numbered_name('_', exclude=names)
    reps = {}
    for b in bound:
        dn = next(dums)
        if hasattr(b, 'as_symbol'):
            d = b.as_symbol(dn)
        else:
            d = Symbol(dn)
        reps[b] = d
    return reps

def _new_MatrixSymbol_as_dummy (self):
    return MatrixDummy(self.name, *self.shape)

def _new_MatrixSymbol_as_symbol (self, name):
    return MatrixSymbol(name, *self.shape)

def fix_Basic_MatrixSymbol_about_dummy ():
    Basic.dummy_eq = _new_Basic_dummy_eq
    Basic.as_dummy = _new_Basic_as_dummy
    Basic.canonical_variables = property(_new_Basic_canonical_variables)
    MatrixSymbol.as_dummy = _new_MatrixSymbol_as_dummy
    MatrixSymbol.as_symbol = _new_MatrixSymbol_as_symbol


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


def make_dummy (expr, name="x", dummy_index=None):
    if expr.is_Matrix:
        return MatrixDummy(name, *expr.shape, dummy_index)
    return Dummy(name)


class MatrixDummy (Dummy, MatrixSymbol):
    is_Dummy = True
    is_commutative = False
    kind: MatrixKind = MatrixKind()

    def __new__ (cls, name=None, n=None, m=None, dummy_index=None, **assumptions):
        if dummy_index is not None:
            assert name is not None, "If you specify a dummy_index, you must also provide a name"

        if name is None:
            name = "Dummy_" + str(Dummy._count)

        if dummy_index is None:
            dummy_index = Dummy._base_dummy_index + Dummy._count
            Dummy._count += 1

        cls._sanitize(assumptions, cls)
        obj = MatrixSymbol.__new__(cls, name, n, m)
        obj.dummy_index = dummy_index

        tmp_asm_copy = assumptions.copy()
        obj._assumptions = StdFactKB(assumptions)
        obj._assumptions._generator = tmp_asm_copy

        return obj
    
    def __getnewargs_ex__(self):
        return ((self.name, *self.shape, self.dummy_index), self.assumptions0)

    @property
    def name(self):
        return self.args[0].name

    @property
    def args (self):
        return super().args + (sympify(self.dummy_index),)

    def matches (self, expr, repl_dict=None, old=False):
        return Basic\
            .matches(self, expr, repl_dict=repl_dict, old=old)

    def xreplace (self, rule):
        return super(Atom, self).xreplace(rule)

    def _hashable_content (self):
        return MatrixSymbol._hashable_content(self) \
            + (self.dummy_index,) \
            + tuple(sorted(self.assumptions0.items()))

    def as_dummy (self):
        return MatrixDummy(self.name, *self.shape)


class MatrixWild (Wild, MatrixSymbol):
    is_Wild = True
    is_commutative = False
    kind: MatrixKind = MatrixKind()

    def __new__(cls, name, n, m, exclude=(), properties=(), **assumptions):
        obj = MatrixSymbol.__new__(cls, name, n, m)
        obj.exclude = tuple([_sympify(x) for x in exclude])
        obj.properties = tuple(properties)

        cls._sanitize(assumptions, cls)
        tmp_asm_copy = assumptions.copy()
        obj._assumptions = StdFactKB(assumptions)
        obj._assumptions._generator = tmp_asm_copy
        
        return obj

    def __getnewargs_ex__(self):
        return ((self.name, *self.shape, self.exclude, self.properties),
                self.assumptions0)

    @property
    def name(self):
        return self.args[0].name

    def _hashable_content (self):
        return MatrixSymbol._hashable_content(self) \
            + (self.exclude, self.properties) \
            + tuple(sorted(self.assumptions0.items()))

    def xreplace (self, rule):
        return super(Atom, self).xreplace(rule)

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


class WildMatrixFunction (WildFunction, MatrixFunction):  # type: ignore
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
