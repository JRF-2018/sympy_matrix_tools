# -*- coding: utf-8 -*-
__version__ = '0.3.3' # Time-stamp: <2022-07-14T13:05:17Z>

import re
from sympy import And, Implies, Predicate, Lambda, AppliedPredicate, Wild,\
    sympify, Basic, Symbol, MatrixSymbol
from sympy.logic.boolalg import BooleanFunction
from sympy.unify import unify
from .matrix_symbol import make_dummy, make_wild
from .unif_tools import rename_bound_symbols_uniquely, construct_term
from .apply import Apply, MatApply, PredApply


ForAll = Predicate("ForAll")


def is_ForAll (z):
    return isinstance(z, AppliedPredicate) and z.function == ForAll


class FrozenGoal (BooleanFunction):
    """Semantically be an identity (lambda x: x)"""
    @classmethod
    def eval (cls, arg):
        pass

    def to_nnf (self, simplify=True):
        return self

    def to_anf (self, deep=True):
        return self

    def melt (self):
        return self.args[0]


def is_symbol_boolean (z):
    try:
        return z.is_boolean
    except:
        pass
    try:
        return z.is_bool
    except:
        pass
    try:
        return z._assumptions["boolean"]
    except:
        pass
    try:
        return z._assumptions["bool"]
    except:
        pass
    return None


def and_implication_normal_form (il, free_wild=True, free_dummy=True):
    if not isinstance(il, Implies) and not is_ForAll(il):
        return rename_bound_symbols_uniquely(il, free_wild=free_wild,
                                             free_dummy=free_dummy)
    prems = []
    vs = []
    while True:
        if isinstance(il, Implies):
            prems.extend(il.args[0:len(il.args) - 1])
            il = il.args[len(il.args) - 1]
        elif is_ForAll(il) \
             and len(il.arguments) == 1 \
             and isinstance(il.arguments[0], Lambda):
            lam = il.arguments[0]
            lvs = lam.args[0]
            larg = lam.args[1]
            rd = {}
            for x in lvs:
                d = make_dummy(x, x.name)
                vs.append(d)
                rd[x] = d
            il = larg.subs(rd)
        else:
            break
    prems = sum([(list(a.args) if isinstance(a, And) else [a])
                 for a in prems], [])
    prems = [and_implication_normal_form(x, free_wild=False,
                                         free_dummy=False) for x in prems]

    r = Implies(And(*prems), il)
    vs = [x for x in vs if x in r.free_symbols]
    if vs:
        vs = sorted(vs, key=lambda x: x.name)
        r = ForAll(Lambda(tuple(vs), r))

    
    return rename_bound_symbols_uniquely(r, free_wild=free_wild,
                                         free_dummy=free_dummy)


def make_outer_forall_wild (z):
    if not is_ForAll(z):
        return z
    if not len(z.arguments) == 1 or not isinstance(z.arguments[0], Lambda):
        return z
    lvs = z.arguments[0].args[0]
    larg = z.arguments[0].args[1]
    d = {}
    for x in lvs:
        d[x] = make_wild(x, x.name)
    return larg.subs(d)


def make_bound_symbols_wild (z, unif=None):
    if unif is None:
        unif = {}
    if z in unif:
        return unif[z]
    if z.is_Atom:
        return z
    inv = {}
    if hasattr(z, 'bound_symbols'):
        b = set(z.bound_symbols)
        f = z.free_symbols & b
        orig_unif = unif
        unif = unif.copy()
        for x in b:
            y = make_wild(x, x.name)
            unif[x] = y
            if x in f:
                inv[y] = x if x not in orig_unif else orig_unif[x]
    return \
        construct_term(z.func,
                       *[make_bound_symbols_wild(a, unif)
                         for a in z.args]).subs(inv)


def get_forall_symbols (z):
    vs = []
    while is_ForAll(z) and isinstance(z.arguments[0], Lambda):
        vs.extend(list(z.arguments[0].args[0]))
        z = z.arguments[0].args[1]
    return tuple(vs), z


def get_syms_prems_concl (z):
    """Get ForAll symbols and premises and conclusion."""
    vs = []
    while is_ForAll(z) and isinstance(z.arguments[0], Lambda):
        vs.extend(list(z.arguments[0].args[0]))
        z = z.arguments[0].args[1]
    if not isinstance(z, Implies):
        return tuple(vs), [], z
    il = z
    prems = []
    while isinstance(il, Implies):
        prems.extend(il.args[0:len(il.args) - 1])
        il = il.args[len(il.args) - 1]
    prems = sum([(list(a.args) if isinstance(a, And) else [a])
                 for a in prems], [])
    return tuple(vs), prems, il


def get_numbered_premise (z, num):
    syms, prems, il = get_syms_prems_concl(z)
    if not prems:
        return None
    if num < len(prems):
        return prems[num]
    return None


def bound_symbols_for_free (z):
    b = set()
    if hasattr(z, 'bound_symbols'):
        b = set(z.bound_symbols)
    d = {}
    if b:
        for v in z.free_symbols:
            d[v] = b
    for x in z.args:
        d1 = bound_symbols_for_free(x)
        for k, v in d1.items():
            if k not in d:
                d[k] = set()
            d[k] = d[k] | v
    return d


def simple_match (m, z):
    if m == z:
        return True
    if isinstance(m, Predicate):
        return m in z.atoms(Predicate)
    if isinstance(m, type):
        return not not z.atoms(m)
    if isinstance(m, Basic):
        return m.match(z) is not None
    if isinstance(m, str):
        vs = z.atoms(Symbol) | z.atoms(MatrixSymbol)
        names = {v.name for v in vs}
        if m in names:
            return True
        names = {re.sub(r'(.*[^01-9])([01-9]+)$', r'\1', s)
                 for s in names}
        return m in names
    return False


def subterm_match (m, z):
    if formula_match(m, z):
        return True
    if z.is_Atom:
        return False
    return any([formula_match(m, a) for a in z.args])


def prem_match (m, z):
    vs, prems, concl = get_syms_prems_concl(z)
    cond, num = separate_cond_num(m)
    if num is None:
        return any([formula_match(cond, a) for a in prems])
    if num < len(prems):
        return formula_match(cond, prems[num])
    return False


def concl_match (m, z):
    vs, prems, concl = get_syms_prems_concl(z)
    return formula_match(m, concl)


def premnum_match (m, z):
    vs, prems, concl = get_syms_prems_concl(z)
    n = len(prems)
    if isinstance(m, int):
        return n == m
    rel, m = m
    if rel == "==":
        return n == m
    if rel == ">=":
        return n >= m
    if rel == ">":
        return n > m
    if rel == "<=":
        return n <= m
    if rel == "<":
        return n < m
    if rel == "!=":
        return n != m
    raise ValueError("%s is not a valid relation" % str(rel))


def not_match (m, z):
    return not formula_match(m, z)


def or_match (m, z):
    return any([formula_match(a, z) for a in m])


_formula_match_table = {
    "prem": prem_match,
    "premnum": premnum_match,
    "sub": subterm_match,
    "concl": concl_match,
    "not": not_match,
    "or": or_match,
    "simple": simple_match
}

def formula_match (m, z):
    if isinstance(m, tuple):
        if len(m) == 2 and m[0] in _formula_match_table:
            m = [m]
        else:
            m = list(m)
    elif not isinstance(m, list):
        m = [m]
    for m1 in m:
        if isinstance(m1, tuple):
            r = False
            if m1[0] in _formula_match_table:
                r = _formula_match_table[m1[0]](m1[1], z)
            if not r:
                return False
        elif not simple_match(m1, z):
            return False
    return True


def separate_cond_num (num):
    if num is None:
        return None, None
    if isinstance(num, int):
        return None, num
    m = num
    if isinstance(m, tuple):
        if len(m) == 2 and m[0] in _formula_match_table:
            m = [m]
        else:
            m = list(m)
    elif not isinstance(m, list):
        m = [m]
    ints = [i for i in m if isinstance(i, int)]
    m = [a for a in m if not isinstance(a, int)]
    i = None
    if ints:
        i = ints[0]
    return m, i


def get_prems_index (prems, index):
    m, i = separate_cond_num(index)
    if m is None and i is None:
        return None
    if m is None and isinstance(i, int):
        if i < len(prems):
            return i
        else:
            return None
    if i is None:
        i = 0
    l = []
    for j, x in enumerate(prems):
        if formula_match(m, x):
            l.append(j)
    if i < len(l):
        return l[i]
    return None


def _contract_double_lift (z, exclude=()):
    exclude = set(exclude)
    
    if isinstance(z, Apply) and isinstance(z.function, Apply) \
       and z.function.function.is_Symbol and z.function.function.is_Wild \
       and z not in exclude:
        symbols = z.function.arguments
        a = set(z.arguments)
        return z.func(z.function.function,
                      *(tuple(z.arguments) + tuple(symbols)))
    if z.is_Atom:
        return z
    return \
        construct_term(z.func,
                       *[_contract_double_lift(a, exclude=exclude)
                         for a in z.args])


def _lift_wild (z, symbols, exclude=()):
    exclude = set(exclude)
    
    if isinstance(z, Apply) and z.function.is_Symbol and z.function.is_Wild \
       and z.function not in exclude:
        a = set(z.arguments)
        return z.func(z.function,
                      *(tuple([_lift_wild(x, symbols, exclude=exclude)
                               for x in z.arguments]) + tuple(symbols)))
    if z.is_Symbol and z.is_Wild and z not in exclude:
        if z.is_Matrix:
            return MatApply(Wild(z.name, exclude=z.exclude), *symbols)
        if is_symbol_boolean(z):
            return PredApply(Wild(z.name, exclude=z.exclude), *symbols)
        return Apply(Wild(z.name, exclude=z.exclude), *symbols)
    if z.is_Atom:
        return z
    inv = {}
    if hasattr(z, 'bound_symbols'):
        b = set(z.bound_symbols)
        f = z.free_symbols & b
        for x in f:
            inv[x] = _lift_wild(x, symbols, exclude)
        exclude = exclude | b
    return \
        construct_term(z.func,
                       *[_lift_wild(a, symbols, exclude)
                         for a in z.args]).subs(inv)


def lift_wild (z, symbols, exclude=()):
    return _contract_double_lift(_lift_wild(z, symbols, exclude=exclude),
                                 exclude=exclude)



def make_proofstate (z, split_prems=False):
    if split_prems:
        z = rename_bound_symbols_uniquely(z)
        f = z.free_symbols
        gvs, prems, il = get_syms_prems_concl(z)
        f = f | set(gvs)
        d = {}
        for v in f:
            d[v] = make_dummy(v, v.name)
        z = z.xreplace(d)
        prems = [y.xreplace(d) for y in prems]
        il = il.xreplace(d)
        return prems, \
            and_implication_normal_form(Implies(il, FrozenGoal(z)),
                                        free_dummy=False)
    return and_implication_normal_form(Implies(z, FrozenGoal(z)))


def _unif_replace (d, inv):
    r = {}
    for k, v in d.items():
        r[k] = v.xreplace(inv)
    return r


def _match_or_unify (z1, z2, exclude=()):
    d = {}
    inv = {}
    for v in set(exclude):
        d[v] = make_dummy(v, v.name)
        inv[d[v]] = v
    z1 = z1.xreplace(d)
    z2 = z2.xreplace(d)
    r = z1.match(z2)
    if r is not None:
        yield _unif_replace(r, inv)
    r = z2.match(z1)
    if r is not None:
        yield _unif_replace(r, inv)
    w = z1.atoms(Wild) | z2.atoms(Wild)
    for r in unify(z1, z2, variables=w):
        yield _unif_replace(r, inv)


def _resolve_implications (x, z, bounded=None, fw=(), gvs=()):
    if bounded is None:
        bounded = {}
    lvs, xprems, xil = get_syms_prems_concl(x)
    z = lift_wild(z, gvs + lvs)
    zvs, zprems, zil = get_syms_prems_concl(z)
    assert not zvs

    for d in _match_or_unify(xil, zil, exclude=fw + gvs + lvs):
        done = True
        for k, v in d.items():
            if k in bounded and bounded[k] & v.free_symbols:
                done = False
                break
        if done:
            nxprems = [y.xreplace(d) for y in xprems]
            nzprems = [y.xreplace(d) for y in zprems]
            if nxprems:
                nzprems = [Implies(And(*nxprems), y) for y in nzprems]
            if lvs:
                nzprems = [ForAll(Lambda(lvs, y)) for y in nzprems]
            yield nzprems, d


def iter_resolve_implications (proofstate, z, index=None, num=0,
                               resolver=_resolve_implications, **kwargs):

    """
    Return the next proofstate which is resolved by proofstate and z.

    I referred to the generic proof assistant Isabelle.

    Parameters
    ==========

    proofstate : an implication whose premises are subgoals and
    conclusion is the goal. Each proofstate is an ordinary theorem
    itself. You can resolve a usual theorem that does not use
    FrozenGoal.

    z : a theorem to be applied to the proofstate.

    index : specifies the subgoal (premise) to have z applied. It may
    be a number or condition to match.

    num: specifies a number to be selected or condition to match.

    Examples
    ========

    >>> from sympy import Symbol, Function, Predicate, Lambda
    >>> from sympy_matrix_tools import ForAll, resolve_implications,\\
    ... print_proofstate

    >>> x = Symbol("x")
    >>> f = Function("f")
    >>> P = Predicate("P")
    >>> Q = Predicate("Q")
    >>> p1 = ForAll(Lambda(x, P(x)))
    >>> p2 = ForAll(Lambda(x, Implies(P(x), Q(f(x)))))
    >>> g3 = ForAll(Lambda(x, Q(f(x))))  # goal

    >>> # At first, make proofstate with specifying a goal.
    >>> z = make_proofstate(g3)
    >>> z
    Implies(Q.ForAll(Lambda(x, Q.Q(f(x)))), FrozenGoal(Q.ForAll(Lambda(x1, Q.Q(f(x1))))))
    >>> # FrozenGoal is semantically an identity (lambda x: x).
    >>> # When printing, it is marked by "*".
    >>> print_proofstate(z)  # Each prem is a subgoal.
    prem 0: Q.ForAll(Lambda(x, Q.Q(f(x))))
    gl*: Q.ForAll(Lambda(x1, Q.Q(f(x1))))
    >>> z = resolve_implications(z, p2)
    >>> print_proofstate(z)
    prem 0: Q.ForAll(Lambda(x, Q.P(x)))
    gl*: Q.ForAll(Lambda(x1, Q.Q(f(x1))))
    >>> z = resolve_implications(z, p1)
    >>> print_proofstate(z)
    gl*: Q.ForAll(Lambda(x, Q.Q(f(x))))
    >>> # To end the proof, melt the FrozenGoal.
    >>> z = melt_theorem(z)
    >>> print_proofstate(z)
    ForAll: x
    gl: Q.Q(f(x))
    >>> z
    Q.ForAll(Lambda(x, Q.Q(f(x))))

    In this framework, Wild symbols behave like meta variables in
    Isabelle. You need to specify boolean=True when creating Wild or
    Symbol as boolean.

    >>> U = Wild("U", boolean=True)
    >>> V = Wild("V", boolean=True)
    >>> W = Wild("W", boolean=True)
    >>> impE = Implies(And(Implies(U, V), U, Implies(V, W)), W)
    >>> impE
    Implies(U_ & (Implies(U_, V_)) & (Implies(V_, W_)), W_)

    """

    gvs = ()
    proofstate, z = rename_bound_symbols_uniquely([proofstate, z],
                                                  free_dummy=False)
    z = make_outer_forall_wild(z)
    z = make_bound_symbols_wild(z)

    fw = set()
    for g in proofstate.atoms(FrozenGoal):
        fw = fw | g.atoms(Wild)
    fw = tuple(fw)
    
    bounded = bound_symbols_for_free(z)
    b1 = bound_symbols_for_free(proofstate)
    for k, v in b1.items():
        if k not in bounded:
            bounded[k] = set()
        bounded[k] = bounded[k] | v

    gvs, prems, il = get_syms_prems_concl(proofstate)

    icond, index = separate_cond_num(index)
    if index is not None and index >= len(prems):
        # raise ValueError("No premise of the index.")
        return
    econd, num = separate_cond_num(num)
    if econd is not None and num is None:
        num = 0

    for i, x in enumerate(prems):
        if icond is None or formula_match(icond, x):
            if isinstance(index, int):
                if index != 0:
                    index -= 1
                    continue
        for r, d in resolver(x, z, bounded=bounded, fw=fw, gvs=gvs, **kwargs):
            nprems = [y.xreplace(d) for y in prems]
            nil = il.xreplace(d)
            nprems = nprems[0:i] + r + nprems[i+1:]
            if nprems:
                r = Implies(And(*nprems), nil)
            else:
                r = nil
            if gvs:
                r = ForAll(Lambda(gvs, r))
            r = and_implication_normal_form(r, free_dummy=False)
            if econd is None or formula_match(econd, r):
                if isinstance(num, int):
                    if num != 0:
                        num -= 1
                        continue
                yield r


def resolve_implications (proofstate, z, index=None, num=0,
                          resolver=_resolve_implications, **kwargs):
    for z1 in iter_resolve_implications(proofstate, z, index=index, num=num,
                                        resolver=resolver,
                                        **kwargs):
        return z1
    return None


def _eresolve_implications (x, z, bounded=None, fw=(), gvs=(), elim_index=0,
                            elim=True):
    if bounded is None:
        bounded = {}
    lvs, xprems, xil = get_syms_prems_concl(x)
    z = lift_wild(z, gvs + lvs)
    zvs, zprems, zil = get_syms_prems_concl(z)
    assert not zvs
    elim_index = get_prems_index(zprems, elim_index)
    if elim_index is None:
        return

    for i, u in enumerate(xprems):
        for d1 in _match_or_unify(u, zprems[elim_index], exclude=fw + gvs + lvs):
            nxprems = [y.xreplace(d1) for y in xprems]
            nzprems = [y.xreplace(d1) for y in zprems]
            nxil = xil.xreplace(d1)
            nzil = zil.xreplace(d1)
            for d2 in _match_or_unify(nxil, nzil, exclude=fw + gvs + lvs):
                mxprems = [y.xreplace(d2) for y in nxprems]
                mzprems = [y.xreplace(d2) for y in nzprems]
                mzprems.pop(elim_index)
                if elim:
                    mxprems = mxprems[0:i] + mxprems[i+1:]
                if mxprems:
                    mzprems = [Implies(And(*mxprems), y) for y in mzprems]
                if lvs:
                    mzprems = [ForAll(Lambda(lvs, y)) for y in mzprems]
                d = d2.copy()
                d.update(d1)
                done = True
                for k, v in d.items():
                    if k in bounded and bounded[k] & v.free_symbols:
                        done = False
                        break
                if done:
                    yield mzprems, d


def iter_eresolve_implications (proofstate, z, index=None, num=0,
                                elim_index=0, elim=True):
    for z1 in iter_resolve_implications(proofstate, z, index=index, num=num,
                                        resolver=_eresolve_implications,
                                        elim_index=elim_index,
                                        elim=elim):
        yield z1


def eresolve_implications (proofstate, z, index=None, num=0,
                           elim_index=0, elim=True):
    return resolve_implications(proofstate, z, index=index, num=num,
                                resolver=_eresolve_implications,
                                elim_index=elim_index,
                                elim=elim)


def _forall_eresolve_implications (x, z, bounded=None, fw=(), gvs=(),
                                   forall_index=None, elim=True):
    if bounded is None:
        bounded = {}
    lvs, xprems, xil = get_syms_prems_concl(x)

    if forall_index is not None:
        m = forall_index
        if isinstance(m, tuple):
            if len(m) == 2 and m[0] in _formula_match_table:
                m = [m]
            else:
                m = list(m)
        elif not isinstance(m, list):
            m = [m]
        m = [ForAll] + m
        forall_index = get_prems_index(xprems, m)
        if forall_index is None:
            return

    for i, u in enumerate(xprems):
        if forall_index is not None and i != forall_index:
            continue
        if is_ForAll(u):
            z = u
            f = z.free_symbols
            z = make_outer_forall_wild(z)
            z = make_bound_symbols_wild(z)
            z = lift_wild(z, gvs + lvs, exclude=f)
            nxprems = xprems.copy()
            nx = xil
            if elim:
                nxprems = nxprems[0:i] + nxprems[i+1:]
            if nxprems:
                nx = Implies(And(*(nxprems + [z])), nx)
            else:
                nx = Implies(z, nx)
            if lvs:
                nx = ForAll(Lambda(lvs, nx))
            yield [nx], {}


def iter_forall_eresolve_implications (proofstate, index=None, num=0,
                                  forall_index=None, elim=True):
    for z1 in iter_resolve_implications(proofstate, sympify(True),
                                        index=index, num=num,
                                        resolver=_forall_eresolve_implications,
                                        forall_index=forall_index,
                                        elim=elim):
        yield z1


def forall_eresolve_implications (proofstate, index=None, num=0,
                                  forall_index=None, elim=True):
    return resolve_implications(proofstate, sympify(True),
                                index=index, num=num,
                                resolver=_forall_eresolve_implications,
                                forall_index=forall_index,
                                elim=elim)


def _dresolve_implications (x, z, bounded=None, fw=(), gvs=(), elim_index=0,
                            elim=True):
    if bounded is None:
        bounded = {}
    lvs, xprems, xil = get_syms_prems_concl(x)
    z = lift_wild(z, gvs + lvs)
    zvs, zprems, zil = get_syms_prems_concl(z)
    assert not zvs
    elim_index = get_prems_index(zprems, elim_index)
    if elim_index is None:
        return

    for i, u in enumerate(xprems):
        for d1 in _match_or_unify(u, zprems[elim_index],
                                  exclude=fw + gvs + lvs):
            nxprems = [y.xreplace(d1) for y in xprems]
            nzprems = [y.xreplace(d1) for y in zprems]
            nxil = xil.xreplace(d1)
            nzil = zil.xreplace(d1)
            nzprems.pop(elim_index)
            nx = nxil
            if elim:
                nxprems = nxprems[0:i] + nxprems[i+1:]
            if nxprems:
                nzprems = [Implies(And(*nxprems), y) for y in nzprems]
                nx = Implies(And(*(nxprems + [nzil])), nx)
            else:
                nx = Implies(nzil, nx)
            if lvs:
                nzprems = [ForAll(Lambda(lvs, y)) for y in nzprems]
                nx = ForAll(Lambda(lvs, nx))
            nzprems = nzprems + [nx]
            done = True
            for k, v in d1.items():
                if k in bounded and bounded[k] & v.free_symbols:
                    done = False
                    break
            if done:
                yield nzprems, d1


def iter_dresolve_implications (proofstate, z, index=None, num=0,
                           elim_index=0, elim=True):
    for z1 in iter_resolve_implications(proofstate, z, index=index, num=num,
                                        resolver=_dresolve_implications,
                                        elim_index=elim_index,
                                        elim=elim):
        yield z1


def dresolve_implications (proofstate, z, index=None, num=0,
                           elim_index=0, elim=True):
    return resolve_implications(proofstate, z, index=index, num=num,
                                resolver=_dresolve_implications,
                                elim_index=elim_index,
                                elim=elim)


def _sresolve_implications (x, z, bounded=None, fw=(), gvs=()):
    if bounded is None:
        bounded = {}
    lvs, xprems, xil = get_syms_prems_concl(x)
    z = lift_wild(z, gvs + lvs)
    zvs, zprems, zil = get_syms_prems_concl(z)
    assert not zvs

    nxprems = xprems
    nzprems = zprems
    nxil = xil
    nzil = zil
    nx = nxil
    if nxprems:
        nzprems = [Implies(And(*nxprems), y) for y in nzprems]
        nx = Implies(And(*(nxprems + [nzil])), nx)
    else:
        nx = Implies(nzil, nx)
    if lvs:
        nzprems = [ForAll(Lambda(lvs, y)) for y in nzprems]
        nx = ForAll(Lambda(lvs, nx))
    nzprems = nzprems + [nx]
    yield nzprems, {}


def iter_sresolve_implications (proofstate, z, index=None, num=0):
    """subgoal_tac"""
    for z1 in iter_resolve_implications(proofstate, z, index=index, num=num,
                                        resolver=_sresolve_implications):
        yield z1


def sresolve_implications (proofstate, z, index=None, num=0):
    """subgoal_tac"""
    return resolve_implications(proofstate, z, index=index, num=num,
                                resolver=_sresolve_implications)


def check_most_trivial (z):
    if z == True:
        return True
    if z == False:
        return False
    lvs, zprems, zil = get_syms_prems_concl(z)
    return any([y for y in zprems if y == zil])


def _remove_trivial_assumptions (x, bounded=None, fw=(), gvs=()):
    if bounded is None:
        bounded = {}
    lvs, xprems, xil = get_syms_prems_concl(x)

    for i, z in enumerate(xprems):
        f = z.free_symbols
        z = make_outer_forall_wild(z)
        z = make_bound_symbols_wild(z)
        z = lift_wild(z, gvs + lvs, exclude=f)
        for d in _match_or_unify(xil, z, exclude=fw + gvs + lvs):
            done = True
            for k, v in d.items():
                if k in bounded and bounded[k] & v.free_symbols:
                    done = False
                    break
            if done:
                yield [], d
        zvs, zprems, zil = get_syms_prems_concl(z)
        assert not zvs
        if zprems:
            for d in _match_or_unify(xil, zil, exclude=fw + gvs + lvs):
                done = True
                for k, v in d.items():
                    if k in bounded and bounded[k] & v.free_symbols:
                        done = False
                        break
                if not done:
                    continue
                nzprems = [y.xreplace(d) for y in zprems]
                nxprems = [y.xreplace(d) for y in xprems]
                if nxprems:
                    nzprems = [Implies(And(*nxprems), y) for y in nzprems]
                # if lvs:
                #     nzprems = [ForAll(Lambda(lvs, y)) for y in nzprems]
                yield nzprems, d


def iter_remove_trivial_assumptions (proofstate, index=None, num=None):

    gvs, prems, il = get_syms_prems_concl(proofstate)

    fw = set()
    for g in proofstate.atoms(FrozenGoal):
        fw = fw | g.atoms(Wild)
    fw = tuple(fw)

    bounded = bound_symbols_for_free(proofstate)
    if not prems:
        yield il
        return
    icond, index = separate_cond_num(index)
    if index is not None and index >= len(prems):
        # raise ValueError("No premise of the index.")
        return
    econd, num = separate_cond_num(num)
    for i, x in enumerate(prems):
        if icond is None or formula_match(icond, x):
            if isinstance(index, int):
                if index != 0:
                    index -= 1
                    continue
        for r, d in _remove_trivial_assumptions(x, bounded=bounded,
                                                fw=fw, gvs=gvs):
            nprems = [y.xreplace(d) for y in prems]
            il = il.xreplace(d)
            nprems = nprems[0:i] + r + nprems[i+1:]
            nprems = [y for y in nprems if not check_most_trivial(nprems)]
            if nprems:
                r = Implies(And(*nprems), il)
            else:
                r = il
            if gvs:
                r = ForAll(Lambda(gvs, r))
            r = and_implication_normal_form(r, free_dummy=False)
            if econd is None or formula_match(econd, r):
                if isinstance(num, int):
                    if num != 0:
                        num -= 1
                        continue
                yield r


def remove_trivial_assumptions (proofstate, index=None, num=None):
    l = []
    for r in iter_remove_trivial_assumptions(proofstate,
                                             index=index, num=num):
        l.append(r)

    if not l:
        return None

    def prems_num (x):
        _, prems, _ = get_syms_prems_concl(x)
        if prems is None:
            return 0
        return len(prems)

    return sorted(l, key=prems_num)[0]

    

def try_remove_trivial_assumptions (proofstate, index=None, num=None):
    z = proofstate
    while True:
        z = remove_trivial_assumptions(proofstate, index, num)
        if z is None or z == proofstate:
            return proofstate
        proofstate = z


def melt_theorem (z):
    if isinstance(z, FrozenGoal):
        z = z.melt()
        z = rename_bound_symbols_uniquely(z, free_dummy=False)
        d = {}
        for v in z.free_symbols:
            if v.is_Dummy:
                d[v] = make_wild(v, v.name)
        z = z.xreplace(d)
        z = rename_bound_symbols_uniquely(z)
        return z
    if z.atoms(FrozenGoal):
        raise ValueError("There remain subgoals.")
    raise ValueError("Already melted.")


def print_proofstate (z):
    gvs = ()
    gvs, prems, il = get_syms_prems_concl(z)
    if gvs:
        print("ForAll:", *gvs)
    if prems:
        for i, x in enumerate(prems):
            print("prem %d:" % i, x)
    if isinstance(il, FrozenGoal):
        print("gl*:", il.args[0])
    else:
        print("gl:", il)
