# -*- coding: utf-8 -*-
__version__ = '0.4.0' # Time-stamp: <2026-10-22T14:29:04Z>

from sympy import Function, MatrixExpr, sympify, MatrixSymbol, Expr, \
    Derivative, Subs, Dummy, ZeroMatrix, MatAdd, MatMul, Symbol, Tuple, \
    Mul, Transpose, transpose
from sympy.matrices.expressions.matexpr import MatrixElement, _LeftRightArgs
from sympy.core.function import UndefinedFunction, Application, AppliedUndef, \
    ArgumentIndexError
from sympy.utilities.iterables import is_sequence, has_dups
from sympy.tensor.array.array_derivatives import ArrayDerivative
from sympy.tensor.array.expressions.array_expressions import _ArrayExpr, \
    ZeroArray, get_rank, get_shape, _array_tensor_product, \
    _array_contraction, _array_add
from sympy.tensor.array.expressions.arrayexpr_derivatives import array_derive
from sympy.tensor.array.expressions.from_matrix_to_array \
    import convert_matrix_to_array

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

    def fdiff (self, argindex=1):
        return _matrix_function_fdiff(self, argindex)

    def _eval_derivative (self, x):
        return _matrix_function_eval_derivative(self, x)


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


# Chain rule for MatrixFunction.
#
# The convention of SymPy's array_derive(expr, x) is that the indices
# of x come first and the indices of expr follow, i.e. the shape of
# the result is x.shape + expr.shape.  So the chain rule for
# F(a_1, ..., a_n) is
#
#   dF/dx = sum_i contract(da_i/dx (x) dF/da_i)
#
# where da_i/dx has the indices [x, a_i], dF/da_i has [a_i, F], and
# the indices of a_i in the middle are contracted.
#
# sympy_matrix_tools_chain_rule_py.py (made by Gemini) assumed the
# opposite order, and contracted wrong indices.

def _fix_matrix_mul (e):
    # Derivative.doit may make a Mul of matrices like Mul(2, X, A, X).
    return e.replace(lambda x: isinstance(x, Mul)
                     and any(a.is_Matrix for a in x.args),
                     lambda x: MatMul(*x.args).doit())


class MatrixDerivative (Derivative, MatrixExpr):
    """Derivative of a matrix expression by a scalar."""

    shape = property(lambda self: self.expr.shape)

    def doit (self, **hints):
        return _fix_matrix_mul(super().doit(**hints))

    def _entry (self, i, j, **kwargs):
        return MatrixElement(self, i, j)


def _subs_new (cls, expr, variables, point, **assumptions):
    # Same as Subs.__new__ except that the placeholders for matrix
    # variables are MatrixSymbol's.  Quoted from SymPy 1.14.
    if not is_sequence(variables, Tuple):
        variables = [variables]
    variables = Tuple(*variables)
    if has_dups(variables):
        raise ValueError("Variables appear more than once: %s" % variables)
    point = Tuple(*(point if is_sequence(point, Tuple) else [point]))
    if len(point) != len(variables):
        raise ValueError('Number of point values must be the same as '
                         'the number of variables.')
    if not point:
        return sympify(expr)
    if isinstance(expr, Subs):
        variables = expr.variables + variables
        point = expr.point + point
        expr = expr.expr
    else:
        expr = sympify(expr)

    def placeholder (name, v):
        if v.is_Matrix:
            return MatrixSymbol(name, *v.shape)
        return Symbol(name)

    pre = "_"
    while 1:
        reps = [(v, placeholder(pre + str(p), v))
                for v, p in zip(variables, point)]
        if any(r in expr.free_symbols and r in variables
               and placeholder(pre + str(point[variables.index(r)]), r) != r
               for _, r in reps):
            pre += "_"
            continue
        break

    obj = Expr.__new__(cls, expr, variables, point)
    obj._expr = expr.xreplace(dict(reps))
    return obj


def _array_derivative_doit (e, **hints):
    # ArrayDerivative.doit uses convert_array_to_matrix, which may
    # remove dimensions of size 1 and change the rank.
    r = e.expr.doit(**hints)
    for v, c in e.variable_count:
        for i in range(c):
            r = array_derive(convert_matrix_to_array(r), v)
    return r


def _subs_doit (self, **hints):
    # Unlike Subs.doit, keep Subs while the derivatives by the
    # variables can not be evaluated, e.g. Derivative(F(xi), xi).
    e = self.expr
    if hints.get('deep', True):
        if isinstance(e, ArrayDerivative) and isinstance(self, ArraySubs):
            e = _array_derivative_doit(e, **hints)
        else:
            e = _fix_matrix_mul(e.doit(**hints))
    if any(set(d.variables) & set(self.variables)
           for d in e.atoms(Derivative)):
        return self.func(e, self.variables, self.point)
    return e.subs(list(zip(self.variables, self.point)))


class MatrixSubs (Subs, MatrixExpr):
    """Subs of a matrix expression."""

    shape = property(lambda self: self.expr.shape)
    __new__ = _subs_new
    doit = _subs_doit

    def _entry (self, i, j, **kwargs):
        return MatrixElement(self, i, j)

    def _eval_derivative (self, s):
        # s is a scalar.
        e = self.expr
        r = []
        if e.has(s):
            r.append(MatrixSubs(e.diff(s), self.variables, self.point))
        for v, p in zip(self.variables, self.point):
            dp = p.diff(s)
            if dp == 0:
                continue
            if v.is_Matrix:
                raise NotImplementedError()
            r.append(dp * MatrixSubs(MatrixDerivative(e, v),
                                     self.variables, self.point))
        if not r:
            return ZeroMatrix(*self.shape)
        if len(r) == 1:
            return r[0]
        return MatAdd(*r)


class ArraySubs (Subs, _ArrayExpr):
    """Subs of an array expression."""

    shape = property(lambda self: self.expr.shape)
    __new__ = _subs_new
    doit = _subs_doit


def _matrix_function_fdiff (expr, argindex):
    from .matrix_symbol import make_dummy

    if not (1 <= argindex <= len(expr.args)):
        raise ArgumentIndexError(expr, argindex)
    ix = argindex - 1
    A = expr.args[ix]
    # For a matrix A, ArrayDerivative(F(A), A) can't be converted by
    # convert_array_to_matrix, so ArraySubs is always used.
    if (isinstance(A, Dummy) or (A.is_Symbol and not A.is_Matrix)) \
       and not any(A in a.free_symbols
                               for i, a in enumerate(expr.args) if i != ix):
        D = A
        e = expr
    else:
        D = make_dummy(A, 'xi_%i' % argindex, dummy_index=hash(A))
        e = expr.func(*(expr.args[:ix] + (D,) + expr.args[ix + 1:]))
    if A.is_Matrix:
        d = ArrayDerivative(e, D, evaluate=False)
        if D != A:
            d = ArraySubs(d, D, A)
    else:
        d = MatrixDerivative(e, D, evaluate=False)
        if D != A:
            d = MatrixSubs(d, D, A)
    return d


class _DirectionDummy (Dummy):
    """The scalar s of a directional derivative d/ds F(a + s da)|_{s=0}."""
    pass


def _matrix_function_eval_derivative (expr, x):
    # x is a scalar.
    if isinstance(x, _DirectionDummy):
        # Otherwise doit() makes a directional derivative endlessly.
        return MatrixDerivative(expr, x, evaluate=False)
    l = []
    for i, a in enumerate(expr.args, start=1):
        da = a.diff(x)
        if da == 0 or (da.is_Matrix and da.is_ZeroMatrix):
            continue
        if not a.is_Matrix:
            df = _matrix_function_fdiff(expr, i)
            # Don't use da * df, because it calls doit().
            l.append(df if da == 1 else MatMul(da, df))
        else:
            # Directional derivative: d/dx F(a) = d/ds F(a + s da)|_{s=0}.
            s = _DirectionDummy('s_%i' % i, dummy_index=hash((a, da)))
            # Don't use a + s * da, because it calls doit().
            e = expr.func(*(expr.args[:i - 1] + (MatAdd(a, MatMul(s, da)),)
                            + expr.args[i:]))
            l.append(MatrixSubs(MatrixDerivative(e, s, evaluate=False),
                                s, 0))
    if not l:
        return ZeroMatrix(*expr.shape)
    if len(l) == 1:
        return l[0]
    return MatAdd(*l)


@array_derive.register(AppliedMatrixUndef)
def _ (expr: AppliedMatrixUndef, x: Expr):
    # x is a matrix.
    rx = get_rank(x)
    rf = get_rank(expr)
    l = []
    for i, a in enumerate(expr.args, start=1):
        if not a.has(x):
            continue
        df = _matrix_function_fdiff(expr, i)
        if a == x:
            l.append(df)
            continue
        if a.is_Matrix:
            da = array_derive(convert_matrix_to_array(a), x)
            ra = get_rank(a)
            tp = _array_tensor_product(da, df)
            tp = _array_contraction(tp, *[(rx + j, rx + ra + j)
                                          for j in range(ra)])
        else:
            da = convert_matrix_to_array(a.diff(x))
            tp = _array_tensor_product(da, convert_matrix_to_array(df))
        l.append(tp)
    if not l:
        return ZeroArray(*(get_shape(x) + expr.shape))
    return _array_add(*l)


def mat_diff (expr, x):
    """Differentiate expr by x.

    It is the same as expr.diff(x), except that it returns the array
    expression as is when it can't be converted to a matrix.

    convert_array_to_matrix of SymPy 1.14, which expr.diff(x) uses,
    has bugs about arrays of rank 4.  It fails for
    (-W * phi(u)).diff(u), where phi is a MatrixFunction, and it
    returns a wrong array for F(X * Y).diff(X).
    """
    if not x.is_Matrix:
        return expr.diff(x)
    from sympy.tensor.array.expressions.from_array_to_matrix \
        import convert_array_to_matrix
    r = array_derive(convert_matrix_to_array(expr), x)
    try:
        c = convert_array_to_matrix(r)
    except Exception:
        return r
    if isinstance(c, MatrixExpr) or not hasattr(c, 'shape'):
        return c
    return r


def _Transpose_fixed_eval_derivative (self, x):
    # x is a scalar.
    return transpose(self.arg._eval_derivative(x))


def fix_Transpose_eval_derivative ():
    """Fix Transpose._eval_derivative of SymPy 1.14.

    The original returns the derivative of the argument without
    transposing it.  This is called when this module is imported.
    """
    Transpose._eval_derivative = _Transpose_fixed_eval_derivative


fix_Transpose_eval_derivative()


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
