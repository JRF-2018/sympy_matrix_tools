# -*- coding: utf-8 -*-
__version__ = '0.2.0' # Time-stamp: <2022-06-19T01:27:52Z>

from sympy import preorder_traversal, Eq, Equivalent, And, Implies
from sympy.unify import unify, rewriterule


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
