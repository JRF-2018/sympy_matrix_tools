# -*- coding: utf-8 -*-
"""Higher-order pattern unification for SymPy terms.

Wild's are metavariables.  Bound variables of Lambda are compared up
to alpha-equivalence.  A metavariable applied to distinct bound
variables (or other symbols), e.g. PredApply(F_, x, y), is a pattern
in the sense of Miller, and it is solved by abstraction:
F_ := Lambda((x, y), t).

Arguments of commutative operators (And, Or, Add, Mul, Equivalent)
are tried in all permutations, so ho_unify() may yield more than one
unifier.
"""
__version__ = '0.7.0' # Time-stamp: <2026-10-22T17:31:49Z>

from itertools import permutations
from sympy import Wild, Lambda, Dummy, And, Or, Add, Mul, Equivalent
from sympy.logic.boolalg import Boolean
from .apply import Apply


_COMMUTATIVE = (And, Or, Equivalent)
MAX_PERMUTATION_ARGS = 5


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
    return False


def _subst (t, sigma):
    if not sigma:
        return t
    return t.xreplace(sigma)


def _is_pattern (t):
    if not isinstance(t, Apply) or not isinstance(t.function, Wild):
        return False
    args = t.arguments
    return len(set(args)) == len(args) \
        and all(a.is_Symbol and not isinstance(a, Wild) for a in args)


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
    xs = tuple(Dummy(a.name if hasattr(a, 'name') else 'x')
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
        ds = tuple(Dummy(getattr(v, 'name', 'b')) for v in sv)
        sb = s.expr.xreplace(dict(zip(sv, ds)))
        tb = t.expr.xreplace(dict(zip(tv, ds)))
        yield from _unify(sb, tb, sigma, local | set(ds))
        return

    if s.is_Atom or t.is_Atom:
        return
    if type(s) != type(t) or len(s.args) != len(t.args):
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
