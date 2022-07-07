# -*- coding: utf-8 -*-
__version__ = '0.2.5' # Time-stamp: <2022-07-07T14:12:08Z>

from sympy import And, Implies, Predicate, Lambda, AppliedPredicate, Wild,\
    sympify

from sympy.unify import unify
from .matrix_symbol import make_dummy, make_wild
from .unif_tools import rename_bound_symbols_uniquely, construct_term
from .apply import Apply, MatApply, PredApply

ForAll = Predicate("ForAll")


def is_ForAll (z):
    return isinstance(z, AppliedPredicate) and z.function == ForAll


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
        return il
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


def _resolve_implications (x, z, bounded=None, gvs=()):
    if bounded is None:
        bounded = {}
    lvs, xprems, xil = get_syms_prems_concl(x)
    z = lift_wild(z, gvs + lvs)
    zvs, zprems, zil = get_syms_prems_concl(z)
    assert not zvs

    for d in _match_or_unify(xil, zil, exclude=gvs + lvs):
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


def resolve_implications (proofstate, z, index=None, goal=False,
                          resolver=_resolve_implications, **kwargs):

    """Return the next proofstate which is resolved by proofstate and z.

    I referred to the generic proof assistant Isabelle.

    Parameters
    ==========

    proofstate : a implication whose premises are subgoals and
    conclusion is the goal.

    z : a theorem to be applied to the proofstate.

    index : specifies the subgoal (premise) to have z applied.

    goal : to make the first goal, specify goal=True.

    Examples
    ========

    >>> from sympy import Symbol, Function, Predicate, Lambda
    >>> from sympy_matrix_tools import ForAll

    >>> x = Symbol("x")
    >>> f = Function("f")
    >>> P = Predicate("P")
    >>> Q = Predicate("Q")
    >>> p1 = ForAll(Lambda(x, P(x)))
    >>> p2 = ForAll(Lambda(x, Implies(P(x), Q(f(x)))))
    >>> g3 = ForAll(Lambda(x, Q(f(x))))  # goal

    >>> z = resolve_implications(g3, p2, goal=True)
    >>> z
    Q.ForAll(Lambda(x, Implies(Q.P(x), Q.Q(f(x)))))
    >>> z = resolve_implications(z, p1)
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
    proofstate, z = rename_bound_symbols_uniquely([proofstate, z])
    z = make_outer_forall_wild(z)
    z = make_bound_symbols_wild(z)

    bounded = bound_symbols_for_free(z)
    b1 = bound_symbols_for_free(proofstate)
    for k, v in b1.items():
        if k not in bounded:
            bounded[k] = set()
        bounded[k] = bounded[k] | v

    gvs, proofstate = get_forall_symbols(proofstate)
    if goal:
        prems = [proofstate]
        il = proofstate
    else:
        vs, prems, il = get_syms_prems_concl(proofstate)
        assert not vs
    for i, x in enumerate(prems):
        if index is not None and i != index:
            continue
        for r, d in resolver(x, z, bounded=bounded, gvs=gvs, **kwargs):
            prems = [y.xreplace(d) for y in prems]
            il = il.xreplace(d)
            pre = prems[0:i]
            post = prems[i+1:]
            prems = pre + r + post
            if prems:
                r = Implies(And(*prems), il)
            else:
                r = il
            if gvs:
                r = ForAll(Lambda(gvs, r))
            return and_implication_normal_form(r)
    return None


def _eresolve_implications (x, z, bounded=None, gvs=(), elim_index=0,
                            elim=True):
    if bounded is None:
        bounded = {}
    lvs, xprems, xil = get_syms_prems_concl(x)
    z = lift_wild(z, gvs + lvs)
    zvs, zprems, zil = get_syms_prems_concl(z)
    assert not zvs
    if len(zprems) < elim_index + 1:
        return

    for i, u in enumerate(xprems):
        for d1 in _match_or_unify(u, zprems[elim_index], exclude=gvs + lvs):
            nxprems = [y.xreplace(d1) for y in xprems]
            nzprems = [y.xreplace(d1) for y in zprems]
            nxil = xil.xreplace(d1)
            nzil = zil.xreplace(d1)
            for d2 in _match_or_unify(nxil, nzil, exclude=gvs + lvs):
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
                yield mzprems, d


def eresolve_implications (proofstate, z, index=None, goal=False,
                           elim_index=0, elim=True):
    return resolve_implications(proofstate, z, index=index, goal=goal,
                                resolver=_eresolve_implications,
                                elim_index=elim_index,
                                elim=elim)


def _forall_eresolve_implications (x, z, bounded=None, gvs=(),
                                   forall_index=None, elim=True):
    if bounded is None:
        bounded = {}
    lvs, xprems, xil = get_syms_prems_concl(x)

    for i, u in enumerate(xprems):
        if is_ForAll(u):
            if forall_index is not None and i != forall_index:
                continue
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
            print("nx", nx)
            yield [nx], {}


def forall_eresolve_implications (proofstate, index=None, goal=False,
                                  forall_index=None, elim=True):
    return resolve_implications(proofstate, sympify(True),
                                index=index, goal=goal,
                                resolver=_forall_eresolve_implications,
                                forall_index=forall_index,
                                elim=elim)


def _dresolve_implications (x, z, bounded=None, gvs=(), elim_index=0,
                            elim=True):
    if bounded is None:
        bounded = {}
    lvs, xprems, xil = get_syms_prems_concl(x)
    z = lift_wild(z, gvs + lvs)
    zvs, zprems, zil = get_syms_prems_concl(z)
    assert not zvs
    if len(zprems) < elim_index + 1:
        return

    for i, u in enumerate(xprems):
        for d1 in _match_or_unify(u, zprems[elim_index], exclude=gvs + lvs):
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
            yield nzprems, d1


def dresolve_implications (proofstate, z, index=None, goal=False,
                           elim_index=0, elim=True):
    return resolve_implications(proofstate, z, index=index, goal=goal,
                                resolver=_dresolve_implications,
                                elim_index=elim_index,
                                elim=elim)


def _remove_trivial_assumptions (x, bounded=None, gvs=()):
    if bounded is None:
        bounded = {}
    lvs, xprems, xil = get_syms_prems_concl(x)

    for i, z in enumerate(xprems):
        f = z.free_symbols
        z = make_outer_forall_wild(z)
        z = make_bound_symbols_wild(z)
        z = lift_wild(z, gvs + lvs, exclude=f)
        for d in _match_or_unify(xil, z, exclude=gvs + lvs):
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
            for d in _match_or_unify(xil, zil, exclude=gvs + lvs):
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


def remove_trivial_assumptions (proofstate, index=None, num=None):

    gvs, prems, il = get_syms_prems_concl(proofstate)
    l = []
    bounded = bound_symbols_for_free(proofstate)
    if not prems:
        return il
    for i, x in enumerate(prems):
        if index is not None and i != index:
            continue
        for r, d in _remove_trivial_assumptions(x, bounded=bounded, gvs=gvs):
            if num is not None:
                if num != 0:
                    num -= 1
                    continue
                num -= 1
            nprems = [y.xreplace(d) for y in prems]
            il = il.xreplace(d)
            pre = nprems[0:i]
            post = nprems[i+1:]
            nprems = pre + r + post
            if nprems:
                r = Implies(And(*nprems), il)
            else:
                r = il
            if gvs:
                r = ForAll(Lambda(gvs, r))
            l.append(and_implication_normal_form(r))
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


def print_proofstate (z):
    gvs = ()
    gvs, prems, il = get_syms_prems_concl(z)
    if gvs:
        print("ForAll:", *gvs)
    if prems:
        for i, x in enumerate(prems):
            print("prem %d:" % i, x)
    print("gl:", il)
