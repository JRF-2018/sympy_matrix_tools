# -*- coding: utf-8 -*-
"""Higher-order pattern unification for SymPy terms.

Wild's are metavariables.  Bound variables of Lambda are compared up
to alpha-equivalence.  A metavariable applied to distinct bound
variables (or other symbols), e.g. PredApply(F_, x, y), is a pattern
in the sense of Miller, and it is solved by abstraction:
F_ := Lambda((x, y), t).

Arguments of commutative operators (And, Or, Add, Mul, MatAdd,
Equivalent) are tried in all permutations, so ho_unify() may yield
more than one unifier.  For associative operators (And, Or, Add, Mul, MatAdd and
MatMul), a metavariable may stand for several arguments, e.g. ?M + ?N
unifies with A + B + C.
"""
__version__ = '0.9.0' # Time-stamp: <2026-11-03T01:08:23Z>

from itertools import permutations, combinations
from sympy import Wild, Lambda, Dummy, And, Or, Add, Mul, Equivalent, \
    MatrixSymbol, MatAdd, MatMul, Transpose
from sympy.logic.boolalg import Boolean
from .apply import Apply
from .matrix_symbol import make_dummy
from .unif_tools import construct_term


_COMMUTATIVE = (And, Or, Equivalent)
MAX_PERMUTATION_ARGS = 5
MAX_AC_ARGS = 8


def eval_matrix_ops (z):
    """Evaluate MatMul, MatAdd and Transpose in z as their operators
    (*, +, .T) do.  Substitutions (xreplace) leave them unevaluated,
    e.g. MatMul(x.T, MatMul(P.T, M, P), x)."""
    if z.is_Atom:
        return z
    args = [eval_matrix_ops(a) for a in z.args]
    if isinstance(z, (MatMul, MatAdd, Transpose)):
        try:
            return z.func(*args).doit(deep=False)
        except Exception:
            pass
    if all(a is b for a, b in zip(args, z.args)):
        return z
    return construct_term(z.func, *args)


def _is_boolean_wild (w):
    # Same as logic_tools.is_symbol_boolean.
    for f in (lambda: w.is_boolean, lambda: w.is_bool,
              lambda: w._assumptions["boolean"],
              lambda: w._assumptions["bool"]):
        try:
            return f()
        except (AttributeError, KeyError):
            pass
    return None


def _is_commutative_op (t):
    if isinstance(t, _COMMUTATIVE):
        return True
    if isinstance(t, (Add, Mul)) and t.is_commutative:
        return True
    if isinstance(t, MatAdd):
        return True
    return False


def _is_ac_op (t):
    # Associative and commutative.
    if isinstance(t, (And, Or, MatAdd)):
        return True
    return isinstance(t, (Add, Mul)) and t.is_commutative


def _subst (t, sigma):
    if not sigma:
        return t
    return eval_matrix_ops(t.xreplace(sigma))


def _is_symbol (a):
    # MatrixSymbol.is_Symbol is False.
    return a.is_Symbol or isinstance(a, MatrixSymbol)


def _is_pattern (t):
    if not isinstance(t, Apply) or not isinstance(t.function, Wild):
        return False
    args = t.arguments
    return len(set(args)) == len(args) \
        and all(_is_symbol(a) and not isinstance(a, Wild) for a in args)


def _ok_value (w, v, local):
    if v.has(w):
        return False                    # occurs check
    if v.free_symbols & local:
        return False                    # a bound variable escapes
    if any(v.has(e) for e in w.exclude):
        return False
    if not all(f(v) for f in w.properties):
        return False
    if _is_boolean_wild(w) and not isinstance(v, Boolean):
        return False
    if w.is_Matrix != getattr(v, 'is_Matrix', False) \
       and not isinstance(v, Lambda):
        return False
    return True


def _bind (w, v, sigma, local):
    v = _subst(v, sigma)
    if v == w:
        return sigma
    if not _ok_value(w, v, local):
        return None
    r = {k: _subst(x, {w: v}) for k, x in sigma.items()}
    r[w] = v
    return r


def _solve_pattern (p, t, sigma, local):
    # p = F(a1, ..., an) with distinct symbols ai.
    xs = tuple(make_dummy(a, a.name if hasattr(a, 'name') else 'x')
               for a in p.arguments)
    body = t.xreplace(dict(zip(p.arguments, xs)))
    # Symbols in local but not in the arguments must not appear.
    return _bind(p.function, Lambda(xs, body), sigma,
                 local - set(p.arguments))


def _unify (s, t, sigma, local):
    s = _subst(s, sigma)
    t = _subst(t, sigma)
    if s == t:
        yield sigma
        return

    if isinstance(s, Wild) and not isinstance(t, Wild):
        r = _bind(s, t, sigma, local)
        if r is not None:
            yield r
        return
    if isinstance(t, Wild):
        r = _bind(t, s, sigma, local)
        if r is not None:
            yield r
        return

    if _is_pattern(s) and not (isinstance(t, Apply)
                               and t.function == s.function):
        r = _solve_pattern(s, t, sigma, local)
        if r is not None:
            yield r
            return
    if _is_pattern(t) and not (isinstance(s, Apply)
                               and s.function == t.function):
        r = _solve_pattern(t, s, sigma, local)
        if r is not None:
            yield r
            return

    if isinstance(s, Lambda) or isinstance(t, Lambda):
        if not (isinstance(s, Lambda) and isinstance(t, Lambda)):
            return
        # Not args, because Lambda(x, x) is IdentityFunction with no args.
        sv, tv = s.signature, t.signature
        if len(sv) != len(tv):
            return
        ds = tuple(make_dummy(v, getattr(v, 'name', 'b')) for v in sv)
        sb = s.expr.xreplace(dict(zip(sv, ds)))
        tb = t.expr.xreplace(dict(zip(tv, ds)))
        yield from _unify(sb, tb, sigma, local | set(ds))
        return

    if s.is_Atom or t.is_Atom:
        return
    if type(s) != type(t):
        return
    if len(s.args) != len(t.args):
        # A metavariable may stand for several arguments, e.g.
        # ?M + ?N and A + B + C, or P.T*?M*P and P.T*A*B*P.
        if len(s.args) > len(t.args):
            s, t = t, s
        if len(t.args) > MAX_AC_ARGS:
            return
        if _is_ac_op(s):
            seen = []
            for r in _unify_ac(s, t, sigma, local):
                if r not in seen:
                    seen.append(r)
                    yield r
        elif isinstance(s, MatMul):
            yield from _unify_assoc(s, t, sigma, local)
        return

    if _is_commutative_op(s) and len(t.args) <= MAX_PERMUTATION_ARGS:
        seen = []
        for targs in permutations(t.args):
            for r in _unify_args(s.args, targs, sigma, local):
                if r not in seen:
                    seen.append(r)
                    yield r
        return
    yield from _unify_args(s.args, t.args, sigma, local)


def _is_var (a):
    return isinstance(a, Wild) or _is_pattern(a)


def _group (t, args):
    if len(args) == 1:
        return args[0]
    return eval_matrix_ops(construct_term(t.func, *args))


def _unify_ac (s, t, sigma, local):
    # s has fewer arguments than t.  Each argument of s that is not a
    # metavariable takes one argument of t, and the metavariables take
    # the rest, at least one each.
    vs = [a for a in s.args if _is_var(a)]
    fixed = [a for a in s.args if not _is_var(a)]
    if not vs:
        return

    def match_fixed (fixed, rest, sigma):
        if not fixed:
            yield rest, sigma
            return
        for i, b in enumerate(rest):
            for r in _unify(fixed[0], b, sigma, local):
                yield from match_fixed(fixed[1:], rest[:i] + rest[i + 1:], r)

    def assign (vs, rest, sigma):
        if len(vs) == 1:
            yield from _unify(vs[0], _group(t, rest), sigma, local)
            return
        n = len(rest)
        for k in range(1, n - len(vs) + 2):
            for idx in combinations(range(n), k):
                g = [rest[i] for i in idx]
                others = [rest[i] for i in range(n) if i not in idx]
                for r in _unify(vs[0], _group(t, g), sigma, local):
                    yield from assign(vs[1:], others, r)

    for rest, r in match_fixed(fixed, list(t.args), sigma):
        if len(rest) >= len(vs):
            yield from assign(vs, rest, r)


def _unify_assoc (s, t, sigma, local):
    # s has fewer factors than t.  A metavariable of s may take a
    # contiguous segment of the factors of t.
    sargs, targs = s.args, t.args

    def go (i, j, sigma):
        if i == len(sargs):
            if j == len(targs):
                yield sigma
            return
        a = sargs[i]
        if _is_var(a):
            for k in range(j + 1, len(targs) - (len(sargs) - i - 1) + 1):
                for r in _unify(a, _group(t, targs[j:k]), sigma, local):
                    yield from go(i + 1, k, r)
        elif j < len(targs):
            for r in _unify(a, targs[j], sigma, local):
                yield from go(i + 1, j + 1, r)

    yield from go(0, 0, sigma)


def _unify_args (sargs, targs, sigma, local):
    if not sargs:
        yield sigma
        return
    for r in _unify(sargs[0], targs[0], sigma, local):
        yield from _unify_args(sargs[1:], targs[1:], r, local)


def ho_unify (s, t, sigma=None):
    """Yield unifiers of s and t as dicts from Wild's to terms.

    >>> from sympy import Symbol, Wild, Lambda
    >>> from sympy_matrix_tools import Apply
    >>> x = Symbol("x")
    >>> F = Wild("F")
    >>> f = Symbol("f")
    >>> list(ho_unify(Lambda(x, Apply(F, x)), Lambda(x, Apply(f, x + 1))))
    [{F_: Lambda(_x, Apply(f, _x + 1))}]
    """
    if sigma is None:
        sigma = {}
    yield from _unify(s, t, sigma, frozenset())
