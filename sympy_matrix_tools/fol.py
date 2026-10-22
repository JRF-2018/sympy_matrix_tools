# -*- coding: utf-8 -*-
__version__ = '0.5.0' # Time-stamp: <2026-10-22T16:34:30Z>

"""First-order logic on the minimal core of logic_tools.

The core of logic_tools has only ForAll, Implies and And (as the
meta-level of Isabelle).  Here the other connectives are defined
impredicatively as in higher-order logic:

  Exists(Lambda(x, P(x)))  ==  ForAll R. (ForAll x. P(x) -> R) -> R
  A | B                    ==  ForAll R. (A -> R) & (B -> R) -> R
  Bot                      ==  ForAll R. R
  neg(A)                   ==  A -> Bot

Each definition is given as a pair of axioms (fold and unfold).  The
classical logic is given by the axiom ccontr: ((A -> Bot) -> Bot) ->
A.  The other rules are theorems proved from them.

SymPy's own false and Not are not used for Bot and neg, because SymPy
evaluates Implies(A, false) to ~A at once.  SymPy's Or is used for
disjunction; its automatic simplification (commutativity etc.) is
classically sound.

Predicates are Symbol's applied by PredApply; use PredApply(P, x, y)
for a binary predicate.
"""

from sympy import Symbol, Lambda, Implies, And, Or, Predicate
from .apply import PredApply
from .logic_tools import ForAll, make_proofstate, melt_theorem, auto, \
    and_implication_normal_form, get_syms_prems_concl, formula_str


Exists = Predicate("Exists")
Bot = Symbol("Bot", boolean=True)


def neg (a):
    """Negation as an abbreviation: neg(A) = Implies(A, Bot)."""
    return Implies(a, Bot)


_x = Symbol("x")
_y = Symbol("y")
_P = Symbol("P")
_A = Symbol("A", boolean=True)
_B = Symbol("B", boolean=True)
_C = Symbol("C", boolean=True)
_R = Symbol("R", boolean=True)


def _ex (P, v=_y):
    return Exists(Lambda(v, PredApply(P, v)))


class Theorem:
    """A named statement with its proof.

    proof is a function that takes the initial proofstate and returns
    a proofstate with no subgoals.  An axiom has no proof.
    """

    def __init__ (self, name, statement, proof=None, axiom=False):
        self.name = name
        self.statement = and_implication_normal_form(statement)
        self.proof = proof
        self.axiom = axiom
        self.verified = axiom

    def verify (self):
        if self.axiom:
            return True
        z = self.proof(make_proofstate(self.statement))
        if z is None:
            raise ValueError("%s: the proof failed." % self.name)
        z = melt_theorem(z)
        if and_implication_normal_form(z) != self.statement:
            raise ValueError("%s: the proof proved %s." % (self.name, z))
        self.verified = True
        return True

    def __repr__ (self):
        kind = "axiom" if self.axiom else "theorem"
        return "%s %s: %s" % (kind, self.name,
                              formula_str(self.statement))


def _st (r):
    return getattr(r, 'statement', r)


def by_auto (*rules, depth=6, max_nodes=2000, safe=()):
    """A proof by the auto tactic with rules (Theorem's or formulas).

    safe is a list of (rule, elim_index) given to auto.  Use SAFE for
    the eliminations of Exists and Or in assumptions.  A proof must
    use only theorems proved before (no circularity).
    """
    def proof (z):
        return auto(z, [_st(r) for r in rules], depth=depth,
                    max_nodes=max_nodes,
                    safe=[(_st(r), i) for r, i in safe])
    return proof


SAFE = []


# Definitions.

ex_fold = Theorem("ex_fold", ForAll(Lambda(_P, Implies(
    ForAll(Lambda(_R, Implies(
        ForAll(Lambda(_x, Implies(PredApply(_P, _x), _R))), _R))),
    _ex(_P)))), axiom=True)

ex_unfold = Theorem("ex_unfold", ForAll(Lambda(_P, Implies(
    _ex(_P),
    ForAll(Lambda(_R, Implies(
        ForAll(Lambda(_x, Implies(PredApply(_P, _x), _R))), _R)))))),
                    axiom=True)

or_fold = Theorem("or_fold", ForAll(Lambda((_A, _B), Implies(
    ForAll(Lambda(_R, Implies(And(Implies(_A, _R), Implies(_B, _R)), _R))),
    Or(_A, _B)))), axiom=True)

or_unfold = Theorem("or_unfold", ForAll(Lambda((_A, _B), Implies(
    Or(_A, _B),
    ForAll(Lambda(_R, Implies(And(Implies(_A, _R), Implies(_B, _R)),
                              _R)))))), axiom=True)

bot_unfold = Theorem("bot_unfold", ForAll(Lambda(_R, Implies(Bot, _R))),
                     axiom=True)

ccontr = Theorem("ccontr", ForAll(Lambda(_A, Implies(neg(neg(_A)), _A))),
                 axiom=True)


# Derived rules.

exI = Theorem("exI", ForAll(Lambda((_P, _x), Implies(PredApply(_P, _x),
                                                     _ex(_P)))),
              by_auto(ex_fold))

exE = Theorem("exE", ForAll(Lambda((_P, _C), Implies(
    And(_ex(_P), ForAll(Lambda(_x, Implies(PredApply(_P, _x), _C)))),
    _C))), by_auto(ex_unfold))

disjI1 = Theorem("disjI1", ForAll(Lambda((_A, _B), Implies(_A, Or(_A, _B)))),
                 by_auto(or_fold))

disjI2 = Theorem("disjI2", ForAll(Lambda((_A, _B), Implies(_B, Or(_A, _B)))),
                 by_auto(or_fold))

disjE = Theorem("disjE", ForAll(Lambda((_A, _B, _C), Implies(
    And(Or(_A, _B), Implies(_A, _C), Implies(_B, _C)), _C))),
                by_auto(or_unfold))

botE = Theorem("botE", ForAll(Lambda(_R, Implies(Bot, _R))),
               by_auto(bot_unfold))

excluded_middle = Theorem("excluded_middle", ForAll(Lambda(_A, Or(
    _A, neg(_A)))), by_auto(ccontr, disjI1, disjI2))


# Not the indices, because And sorts the premises.
SAFE.extend([(exE, Exists), (disjE, Or)])

AXIOMS = [ex_fold, ex_unfold, or_fold, or_unfold, bot_unfold, ccontr]
THEOREMS = [exI, exE, disjI1, disjI2, disjE, botE, excluded_middle]
