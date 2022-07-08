# -*- coding: utf-8 -*-
__version__ = '0.3.0' # Time-stamp: <2022-07-08T17:06:49Z>

import re
from sympy import preorder_traversal, Eq, Equivalent, And, Implies
from sympy import Sum, Product, ITE, Piecewise, Not, Lambda, Basic
from sympy.unify import unify, rewriterule
from sympy import Symbol, MatrixSymbol, Wild, Dummy
from sympy.unify.usympy import basic_new_legal, eval_false_legal, illegal


def apply_rewriterules (rules, z):
    l = []
    for x in preorder_traversal(z):
        for r in rules:
            for y in r(x):
                l.append((x, y))
                break
    return z.xreplace(dict(l))


def rewrite_subterm (dest, from_term, to_term, **kwargs):
    l = []
    rl = rewriterule(from_term, to_term, **kwargs)
    for x in preorder_traversal(dest):
        for y in rl(x):
            l.append((x, y))
            break
    return dest.xreplace(dict(l))


def rewriterules_from_eqs (eqs, variables=(), **kwargs):
    l = []
    vs = []
    vd = []
    for x in eqs:
        if isinstance(x, Eq) or isinstance(x, Equivalent):
            p = x.args[len(x.args) - 1]
            for q in x.args[0:len(x.args) - 1]:
                l.append(rewriterule(q, p, variables=variables, **kwargs))
        else:
            raise ValueError("%s is not Eq or Equivalent" % str(x))
    return l


def try_rewriterule (source, target, variables=(),
                     condition=None, assume=None):
    rl = rewriterule(source, target, variables=variables,
                     condition=condition, assume=assume)
    def apply_rl (z):
        for x in rl(z):
            return x
        return z
    return apply_rl


## resolve modus ponens.
def resolve_mp (impl, conds, variables=(), complete=False, simplify=False):
    if simplify is None:
        simplify = lambda x: x
    elif not simplify:
        simplify = lambda x: x.doit(deep=True)
    elif simplify == True:
        simplify = lambda x: x.simplify()
    
    cs = sum([(list(a.args) if isinstance(a, And) else [a]) for a in conds], [])
    il = impl
    prems = []
    while isinstance(il, Implies):
        prems.extend(il.args[0:len(il.args) - 1])
        il = il.args[len(il.args) - 1]
    prems = sum([(list(a.args) if isinstance(a, And) else [a]) for a in prems], [])
    while prems and cs:
        c = cs.pop(0)
        u = None
        for p in prems:
            for x in unify(c, p, variables=variables):
                u = x
                break
        if u is not None:
            new_prems = []
            c = simplify(c.subs(u))
            il = simplify(il.subs(u))
            cs = [simplify(x.subs(u)) for x in cs]
            cs = [x for x in cs if x != True]
            for p in prems:
                q = p.subs(u)
                if q != c and q != True:
                    new_prems.append(q)
            prems = new_prems
    if complete and prems:
        raise ValueError("Unresolved Premises: %s" % str(prems))
    if prems:
        return Implies(And(*prems), il)
    return il


def renamed_symbol (z, new_name):
    if z.is_Dummy:
        if hasattr(z, 'as_symbol'):
            return z.as_symbol(new_name)
        elif isinstance(z, MatrixSymbol):
            return MatrixSymbol(new_name, *z.shape)
        else:
            return Symbol(new_name, **z.assumptions0)
    elif hasattr(z, 'renamed'):
        return z.renamed(new_name)
    elif isinstance(z, MatrixSymbol):
        return MatrixSymbol(new_name, *z.shape)
    elif isinstance(z, Wild):
        return z.func(new_name, z.exclude, z.properties, **z.assumptions0)
    elif isinstance(z, Symbol):
        return z.func(new_name, **z.assumptions0)
    else:
        raise TypeError("%s is not a Symbol." % str(z))


def _numbered_name (prefix='x', start=0, exclude=()):
    exclude = set(exclude or [])
    while True:
        if start == 0 and prefix != "":
            name = prefix
        else:
            name = '%s%s' % (prefix, start)
        if name not in exclude:
            yield name
        start += 1


def clean_addressed_symbols (z, exclude=()):
    if isinstance(z, list) or isinstance(z, tuple):
        zl = z
        ret = 1
    else:
        zl = [z]
        ret = 0
    a = set()
    for z in zl:
        a = a | z.atoms(Symbol) | z.atoms(MatrixSymbol)

    names = {x.name for x in a} | {x.name for x in exclude}
    q = re.compile(r'(.*)@\[([01-9,]*)\]@([01-9]+)$')
    d = {}
    dl = []
    for x in a:
        m = q.match(x.name)
        if m:
            n = m.group(1)
            l = [int(y) for y in m.group(2).split(',')]
            r = int(m.group(3))
            d[x] = (n, '[%s]@%s' % (','.join([str(i) for i in l]),
                                    str(r)))
    l = sorted(sorted(d.keys(), key=lambda x: d[x][1]),
               key=lambda x: len(d[x][1]))
    for x in l:
        n = d[x][0]
        nn = _numbered_name(n)
        y = next(nn)
        while y in names:
            y = next(nn)
        names.add(y)
        y = renamed_symbol(x, y)
        dl.append((x, y))

    r = [z.xreplace(dict(dl)) for z in zl]
    if ret == 1:
        return r
    else:
        return r[0]


def construct_term (func, *args):
    if any(issubclass(func, cls) for cls in eval_false_legal):
        return func(*args, evaluate=False)
    elif any(issubclass(func, cls) for cls in basic_new_legal):
        return Basic.__new__(func, *args)
    else:
        return func(*args)


def _rename_bound_symbols_uniquely (term, addr, unif, exclude=()):
    if term in unif:
        return unif[term]
    if term.is_Atom:
        return term
    inv = {}
    if hasattr(term, 'bound_symbols'):
        b = set(term.bound_symbols)
        f = term.free_symbols & b
        if unif is None:
            orig_unif = {}
            unif = {}
        else:
            orig_unif = unif
            unif = unif.copy()
        for i, x in enumerate(b):
            if x in exclude:
                continue
            m = re.match(r'(.*[^01-9])?([01-9]+)$', x.name)
            if m:
                xr = m.group(1) or ""
            else:
                xr = x.name
            y = renamed_symbol(x, "%s@[%s]@%s" %
                               (xr,
                                ",".join([str(j) for j in addr]),
                                str(i)))
            unif[x] = y
            if x in f:
                inv[y] = x if x not in orig_unif else orig_unif[x]
    return \
        construct_term(term.func,
                       *[_rename_bound_symbols_uniquely(
                           a, addr + [i], unif, exclude=exclude)
                         for i, a in enumerate(term.args)]).subs(inv)


def rename_bound_symbols_uniquely (term_or_list_of_list,
                                   free_wild=True, free_dummy=True,
                                   exclude=()):
    """
    Return a term, a list of terms or a list of lists of terms whose
    bound symbols are renamed uniquely.

    Elements of the outer list don't share free symbols.  Elements of
    inner lists share free symbols.

    Parameters
    ==========

    term_or_list_of_list : term or list pf terms or list of lists of terms.

    free_wild: True when handling free wild symbols like bound symbols.

    free_dummy: True when handling free dummy symbols like bound symbols.

    exclude: excluding symbols of renaming.

    Examples
    ========

    >>> from sympy import Symbol, Function, Lambda
    >>> from sympy_matrix_tools import rename_bound_symbols_uniquely

    >>> x = Symbol("x")
    >>> x1 = Symbol("x1")
    >>> f = Function("f")
    >>> rename_bound_symbols_uniquely(x + x1 + f(Lambda(x1, x1 + 1)))
    x + x1 + f(Lambda(x2, x2 + 1))

    """
    l = term_or_list_of_list
    if l == [] or l == [[]]:
        return l
    if not (isinstance(l, list) or isinstance(l, tuple)):
        l = [[l]]
        ret = 0
    elif len(l) < 1 or not (isinstance(l[0], list) or isinstance(l, tuple)):
        l = [[z] for z in l]
        ret = 1
    else:
        ret = 2

    l = [[clean_addressed_symbols(z, exclude=exclude) for z in l1] for l1 in l]
    r = []
    for i1, l1 in enumerate(l):
        a = set()
        for z in l1:
            a = a | z.free_symbols
        a1 = set()
        if free_wild:
            a1 = a1 | {x for x in a if x.is_Wild}
        if free_dummy:
            a1 = a1 | {x for x in a if x.is_Dummy}
        unif = {}
        for i, x in enumerate(a1):
            if x in exclude:
                continue
            m = re.match(r'(.*[^01-9])?([01-9]+)$', x.name)
            if m:
                xr = m.group(1) or ""
            else:
                xr = x.name
            y = renamed_symbol(x, "%s@[%s]@%s" % (xr, str(i1), str(i)))
            unif[x] = y
        r1 = []
        for i2, z in enumerate(l1):
            r2 = _rename_bound_symbols_uniquely(
                z, [i1, i2], unif, exclude=exclude)
            r1.append(r2)
        r.append(r1)

    l = r
    r = []
    exclude = set(exclude)
    for l1 in l:
        r1 = clean_addressed_symbols(l1, exclude)
        for r2 in r1:
            exclude = exclude | r2.atoms(Symbol) | r2.atoms(MatrixSymbol)
        r.append(r1)

    if ret == 2:
        return r
    elif ret == 1:
        return [x[0] for x in r]
    else:
        return r[0][0]


def _remove_conds (conds, vs):
    if isinstance(conds, set):
        return conds - vs
    else:
        conds = conds.copy()
        while conds and conds[0] in vs:
            conds.pop(0)
        return conds


def conditional_apply (z, func, conds, match=None):
    if not conds:
        if match is None:
            return func(z)
        if callable(match) and not isinstance(match, Basic):
            if match(z):
                return func(z)
        if z.match(match) is not None:
            return func(z)
        if z.is_Atom:
            return z
        return construct_term(z.func,
                              *[conditional_apply(x, func, conds,
                                                  match=match)
                                for x in z.args])
    if z.is_Atom:
        return z
    
    if isinstance(z, Sum) or isinstance(z, Product):
        b = set(z.bound_symbols)
        c = _remove_conds(conds, b)
        return construct_term(z.func,
                              conditional_apply(z.function, func, c,
                                                match=match),
                              *[conditional_apply(x, func, conds,
                                                  match=match)
                                for x in z.limits])
    if isinstance(z, Lambda):
        b = set(z.bound_symbols)
        c = _remove_conds(conds, b)
        return construct_term(z.func,
                              conditional_apply(z.args[0], func, conds,
                                                match=match),
                              conditional_apply(z.args[1], func, c,
                                                match=match))

    if isinstance(z, ITE):
        a, b, c = z.args
        c1 = _remove_conds(conds, {a})
        c2 = _remove_conds(conds, {Not(a)})
        return construct_term(z.func,
                              conditional_apply(a, func, conds,
                                                match=match),
                              conditional_apply(b, func, c1,
                                                match=match),
                              conditional_apply(c, func, c2,
                                                match=match))
    if isinstance(z, Piecewise):
        l = []
        orig_conds = conds
        for a, b in z.args:
            c = _remove_conds(conds, {b})
            l.append((conditional_apply(a, func, c, match=match),
                      conditional_apply(b, func, orig_conds, match=match)))
            conds = _remove_conds(conds, {Not(b)})
        return construct_term(z.func, *l)

    if isinstance(z, Implies):
        prems = z.args[0:len(z.args) - 1]
        x = z.args[len(z.args) - 1]
        prems = sum([(list(a.args) if isinstance(a, And) else [a]) \
                     for a in prems], [])

        c = _remove_conds(conds, set(prems))
        return construct_term(z.func,
                              *[conditional_apply(y, func, conds,
                                                  match=match)
                                for y in z.args[0:len(z.args) - 1]],
                              conditional_apply(x, func, c, match=match))

    if hasattr(z, 'bound_symbols'):
        b = set(z.bound_symbols)
        c = _remove_conds(conds, b)
        return construct_term(z.func,
                              *[conditional_apply(x, func, c,
                                                  match=match)
                                for x in z.args])

    return construct_term(z.func,
                          *[conditional_apply(x, func, conds,
                                              match=match)
                            for x in z.args])
