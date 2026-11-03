# -*- coding: utf-8 -*-
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

Theorem's are proved by tactics: by_auto (proof search), by_rule,
by_rewrite (by an equation), by_cas (a computation of SymPy) and
by_steps, or forward by a structured Proof (have ... show).  The
statements may be about matrices, and Derivation.theorem turns the
facts of a derivation ledger into Theorem's.  explain() shows what a
theorem rests on.  See docs/fol.md.
"""
__version__ = '0.10.0' # Time-stamp: <2026-11-03T03:41:31Z>

from sympy import Symbol, Lambda, Implies, And, Or, Predicate, \
    MatrixSymbol, Wild, Dummy, sympify
from sympy import preorder_traversal
from .apply import PredApply, Apply
from .matrix_symbol import fix_Basic_MatrixSymbol_about_dummy
from .unif_tools import rename_bound_symbols_uniquely
from .logic_tools import ForAll, is_ForAll, make_proofstate, \
    melt_theorem, auto, and_implication_normal_form, formula_str, \
    get_forall_symbols, get_syms_prems_concl, resolve_implications, \
    rewrite_goal, rewrite_goal_by_assumption, cas_goal, \
    _split_concl, _map_terms, is_symbol_boolean
from sympy.logic.boolalg import Boolean


# Formulas that bind matrices (e.g. ForAll M. psd(M)) need MatrixDummy
# in the canonical forms of SymPy.
fix_Basic_MatrixSymbol_about_dummy()


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


def _shape_symbols (exprs):
    # The symbols in the dimensions of matrices, and the predicates and
    # functions applied by PredApply and Apply.
    s = set()
    for e in exprs:
        for m in preorder_traversal(e):
            if getattr(m, 'is_Matrix', False) and hasattr(m, 'shape'):
                for d in m.shape:
                    s |= sympify(d).free_symbols
            if isinstance(m, Apply):
                s.add(m.function)
    return s


def _variables (exprs, exclude=()):
    # Free Symbol's and MatrixSymbol's, without dimensions, predicates,
    # Wild's and Dummy's, sorted by name.
    exclude = set(exclude) | _shape_symbols(exprs)
    vs = set()
    for e in exprs:
        vs |= {v for v in e.free_symbols | e.atoms(MatrixSymbol)
               if isinstance(v, (Symbol, MatrixSymbol))
               and not isinstance(v, (Wild, Dummy))
               and v not in exclude and v in e.free_symbols}
    return sorted(vs, key=lambda v: str(v.name))


class Theorem:
    """A named statement with its proof.

    statement is a formula, or a string in the logical notation (see
    logic_parser.parse_formula, to which symbols is given).  proof is
    a function (a tactic such as by_auto) that takes the initial
    proofstate and returns a proofstate with no subgoals.  An axiom
    has no proof.

    kind is "axiom", "theorem", "instance" (made by instance()) or
    "ledger" (a fact of a Derivation, made by Derivation.theorem()).
    uses are the theorems that the proof uses; they are taken from the
    tactic.  dependencies() and explain() show what a theorem rests
    on.
    """

    def __init__ (self, name, statement, proof=None, axiom=False,
                  symbols=None, kind=None, uses=()):
        if isinstance(statement, str):
            # e.g. "∀x. P(x) → (∃y. P(y))"; see parse_formula.
            from .logic_parser import parse_formula
            statement = parse_formula(statement, symbols=symbols)
        self.name = name
        self.statement = and_implication_normal_form(statement)
        self.proof = proof
        self.axiom = axiom or kind == "axiom"
        self.kind = kind or ("axiom" if axiom else "theorem")
        self.uses = list(uses) + [u for u in getattr(proof, 'uses', ())
                                  if isinstance(u, Theorem)]
        self.verified = self.axiom
        self.ledger = None      # (Derivation, fact name) for "ledger"

    def verify (self):
        """Run the proof (once) and check that it proves the
        statement.  The theorems used are verified first.  Raise
        ValueError if the proof fails."""
        if self.verified:
            return True
        for u in self.uses:
            u.verify()
        if self.kind == "instance":
            self.verified = True
            return True
        if self.kind == "ledger":
            d, name = self.ledger
            ok, msg = d._fact_status(name)
            if not ok:
                raise ValueError("%s: %s" % (self.name, msg))
            self.verified = True
            return True
        z = self.proof(make_proofstate(self.statement))
        if z is None:
            raise ValueError("%s: the proof failed." % self.name)
        z = melt_theorem(z)
        if and_implication_normal_form(z) != self.statement:
            raise ValueError("%s: the proof proved %s." % (self.name, z))
        self.verified = True
        return True

    def variables (self):
        """The outer bound variables of the statement."""
        return get_forall_symbols(self.statement)[0]

    def instance (self, name=None, over=None, **subs):
        """The instance where the bound variables named by the keys
        of subs are replaced by the values (forall elimination).

        The new free symbols of the values are bound again by ForAll
        (forall introduction), unless over gives the symbols to bind.
        An instance is useful when unification cannot find the terms,
        e.g. ForAll x. phi(x) instantiated with x := P*x.

        >>> x, y = Symbol("x"), Symbol("y")
        >>> P = Symbol("P")
        >>> th = Theorem("all_P", ForAll(Lambda(x, PredApply(P, x))),
        ...              axiom=True)
        >>> th.instance(x=2*y)
        instance all_P[x := 2*y]: ∀y. P(2*y)
        """
        vs, body = get_forall_symbols(self.statement)
        byname = {str(v.name): v for v in vs}
        d = {}
        for k, t in subs.items():
            if k not in byname:
                raise ValueError("%s has no variable %s; the variables "
                                 "are %s." % (self.name, k,
                                              ", ".join(byname)))
            v = byname[k]
            t = sympify(t)
            if getattr(v, 'is_Matrix', False) \
               != bool(getattr(t, 'is_Matrix', False)) \
               or (getattr(v, 'is_Matrix', False)
                   and tuple(v.shape) != tuple(t.shape)):
                raise ValueError("%s := %s: the shapes differ."
                                 % (k, t))
            # A Symbol is a Boolean in SymPy; see its assumption.
            if is_symbol_boolean(v) and not (
                    isinstance(t, Boolean)
                    and (not t.is_Symbol or is_symbol_boolean(t))):
                raise ValueError("%s := %s: %s is a formula." % (k, t, k))
            d[v] = t
        rest = tuple(v for v in vs if v not in d)
        # Rename the inner bound symbols away from the new terms.
        body = rename_bound_symbols_uniquely([[body] + list(d.values())],
                                             free_wild=False,
                                             free_dummy=False)[0][0]
        body = body.xreplace(d)
        if over is None:
            over = _variables(list(d.values()),
                              exclude=self.statement.free_symbols)
        over = [v for v in over if v not in rest]
        bound = rest + tuple(over)
        st = ForAll(Lambda(bound, body)) if bound else body
        if name is None:
            name = "%s[%s]" % (self.name, ", ".join(
                "%s := %s" % (str(v.name), t) for v, t in d.items()))
        return Theorem(name, st, kind="instance", uses=[self])

    def cas (self, func="doit", name=None):
        """The theorem whose terms are computed by SymPy (a function or
        a method name such as "doit" or "expand"), like simplified of
        Isabelle.  The trusted base includes SymPy.

        >>> x, y = Symbol("x"), Symbol("y")
        >>> P = Symbol("P")
        >>> th = Theorem("all_P", ForAll(Lambda(x, PredApply(P, x))),
        ...              axiom=True)
        >>> th.instance(x=(y + 1)**2 - y**2).cas("expand")
        instance all_P[x := -y**2 + (y + 1)**2].expand: ∀y. P(2*y + 1)
        """
        f = (lambda e: getattr(e, func)()) if isinstance(func, str) \
            else func
        fname = func if isinstance(func, str) \
            else getattr(func, '__name__', 'cas')
        st = _map_terms(self.statement, f)
        return Theorem(name or "%s.%s" % (self.name, fname), st,
                       kind="instance", uses=[self])

    def dependencies (self):
        """All the theorems that this theorem rests on, the used ones
        first."""
        seen = []
        def visit (t):
            for u in t.uses:
                if u not in seen:
                    visit(u)
                    seen.append(u)
        visit(self)
        return seen

    def explain (self):
        """A text that shows the statement, how it is proved, and the
        axioms and ledger facts that it rests on."""
        lines = [repr(self)]
        if self.kind == "theorem":
            how = getattr(self.proof, 'desc', None)
            if how:
                lines.append("  proof: %s" % how)
        if self.kind == "ledger":
            d, name = self.ledger
            lines.append("  by %s" % d._fact_status(name)[1])
        deps = self.dependencies()
        for kind, label in (("axiom", "axioms"),
                            ("ledger", "ledger facts")):
            l = [u for u in deps if u.kind == kind]
            if l:
                lines.append("  %s: %s" % (label,
                                           ", ".join(u.name for u in l)))
        l = [u for u in deps if u.kind == "theorem"
             and getattr(u, 'lemma_of', None) is None]
        if l:
            lines.append("  theorems: %s" % ", ".join(u.name for u in l))
        lines.append("  verified: %s" % ("yes" if self.verified else "no"))
        return "\n".join(lines)

    def __repr__ (self):
        return "%s %s: %s" % (self.kind, self.name,
                              formula_str(self.statement))


def _st (r):
    return getattr(r, 'statement', r)


def _name (r):
    return r.name if isinstance(r, Theorem) else formula_str(r)


def _rules (rules):
    # A theorem whose conclusion is an And or an Equivalent is used as
    # its parts.
    l = []
    for r in rules:
        s = _st(r)
        sp = _split_concl(s)
        if sp is None:
            l.append(s)
        else:
            l.extend(and_implication_normal_form(x) for x in sp)
    return l


def _tactic (f, desc, uses):
    f.desc = desc
    f.uses = list(uses)
    return f


def by_auto (*rules, depth=6, max_nodes=2000, safe=()):
    """A proof by the auto tactic with rules (Theorem's or formulas).

    A rule whose conclusion is an And or an Equivalent (a definition
    A <-> B) is used as its parts.  safe is a list of (rule,
    elim_index) given to auto.  Use SAFE for the eliminations of
    Exists and Or in assumptions.  A proof must use only theorems
    proved before (no circularity).
    """
    def proof (z):
        return auto(z, _rules(rules), depth=depth,
                    max_nodes=max_nodes,
                    safe=[(_st(r), i) for r, i in safe])
    return _tactic(proof, "auto(%s)" % ", ".join(_name(r) for r in rules),
                   list(rules) + [r for r, i in safe])


def by_rule (th, index=0, num=0):
    """A step that resolves the subgoal of the index with the theorem
    th (backward reasoning, rule of Isabelle).  It may leave subgoals,
    so use it in by_steps."""
    def proof (z):
        k = num
        for r in _rules([th]):
            y = resolve_implications(z, r, index=index, num=k)
            if y is not None:
                return y
        return None
    return _tactic(proof, "rule(%s)" % _name(th), [th])


def by_rewrite (eq, index=0, reverse=False, num=0, occurrence=None):
    """A step that rewrites the conclusion of the subgoal of the index
    by the equation theorem eq (see logic_tools.rewrite_goal).  All
    occurrences of the matched subterm are rewritten, or only the one
    of the number occurrence (from 0).  It may leave subgoals, so use
    it in by_steps."""
    def proof (z):
        return rewrite_goal(z, _st(eq), index=index, reverse=reverse,
                            num=num, occurrence=occurrence)
    opt = "" if occurrence is None else ", occurrence=%d" % occurrence
    return _tactic(proof, "rewrite%s(%s%s)" % ("_rev" if reverse else "",
                                               _name(eq), opt), [eq])


def by_rewrite_asm (num=0, index=0, reverse=False, occurrence=None):
    """A step that rewrites the conclusion of the subgoal of the index
    by its own premise, the num-th equation among its premises, e.g.
    an induction hypothesis (see
    logic_tools.rewrite_goal_by_assumption)."""
    def proof (z):
        return rewrite_goal_by_assumption(z, index=index, num=num,
                                          reverse=reverse,
                                          occurrence=occurrence)
    return _tactic(proof, "rewrite%s_asm(%d)" % ("_rev" if reverse else "",
                                                 num), [])


def by_induct (rule, var, index=0):
    """A step of induction on the bound variable named var of the
    subgoal of the index, like induct of Isabelle.

    rule is an induction theorem (usually an axiom of your own) with
    one predicate variable, e.g.

      nat_induct: ∀Phi k. Phi(0) ∧ (∀j. Phi(j) → Phi(j + 1)) → Phi(k)

    The predicate variable is instantiated with Lambda(var, C), where C
    is the conclusion of the subgoal (with the premises that contain
    var), and the subgoal is resolved with the instance.  The subgoals
    of the base case and of the step remain.
    """
    if not isinstance(rule, Theorem):
        rule = Theorem("induction", rule, axiom=True)
    pvs = [v for v in rule.variables()
           if any(a.function == v for a in rule.statement.atoms(PredApply))]
    if len(pvs) != 1:
        raise ValueError("%s must have one predicate variable." % rule.name)
    def proof (z):
        gvs, prems, il = get_syms_prems_concl(z)
        if index >= len(prems):
            return None
        lvs, xprems, xil = get_syms_prems_concl(prems[index])
        ks = [v for v in lvs if str(v.name) == var]
        if not ks:
            return None
        k = ks[0]
        dep = [p for p in xprems if p.has(k)]
        body = Implies(And(*dep), xil) if dep else xil
        inst = rule.instance(**{str(pvs[0].name): Lambda(k, body)})
        return resolve_implications(z, inst.statement, index=index)
    return _tactic(proof, "induct(%s, %s)" % (rule.name, var), [rule])


def by_cas (func="doit", index=0):
    """A step that applies a computation of SymPy (a function or a
    method name such as "doit" or "expand") to the terms in the
    conclusion of the subgoal of the index (see
    logic_tools.cas_goal).  It fails if nothing changes."""
    def proof (z):
        return cas_goal(z, func, index=index)
    name = func if isinstance(func, str) else getattr(func, '__name__',
                                                       str(func))
    return _tactic(proof, "cas(%s)" % name, [])


def by_steps (*tactics):
    """A proof by the tactics applied in order.  It fails if a tactic
    fails.  The last tactic must solve all the subgoals."""
    def proof (z):
        for t in tactics:
            if z is None:
                return None
            z = t(z)
        return z
    uses = sum([list(getattr(t, 'uses', ())) for t in tactics], [])
    return _tactic(proof, "; ".join(getattr(t, 'desc', '?')
                                    for t in tactics), uses)


def _peel (z):
    # Fixed variables, assumptions and the conclusion, as in the normal
    # form.
    fixed, prems = [], []
    while True:
        if is_ForAll(z) and isinstance(z.arguments[0], Lambda):
            fixed.extend(z.arguments[0].signature)
            z = z.arguments[0].expr
        elif isinstance(z, Implies):
            for a in z.args[:-1]:
                prems.extend(a.args if isinstance(a, And) else [a])
            z = z.args[-1]
        else:
            return fixed, prems, z


class Proof:
    """A structured proof in forward style, like Isar of Isabelle.

    The statement ForAll xs. H1 & ... -> C is proved by fixing xs,
    assuming H1, ... (named h1, h2, ... or by assume) and deriving
    facts by have(name, formula, by, using) until show(by, using)
    derives C.  Each step is a lemma ForAll xs. (used facts -> formula)
    proved by its tactic at once, so a wrong step raises ValueError
    where it is.  show() returns the theorem, whose proof combines the
    lemmas.

    Formulas may be strings; the names of xs in them are the fixed
    variables, and symbols gives the other names (see parse_formula).
    """

    def __init__ (self, name, statement, symbols=None, assume=None):
        from .logic_parser import parse_formula
        self.symbols = dict(symbols or {})
        if isinstance(statement, str):
            statement = parse_formula(statement, symbols=self.symbols)
        self.name = name
        self.statement = statement
        self.fixed, prems, self.concl = _peel(statement)
        for v in self.fixed:
            self.symbols[str(v.name)] = v
        names = list(assume) if assume is not None \
            else ["h%d" % (i + 1) for i in range(len(prems))]
        if len(names) != len(prems):
            raise ValueError("The statement has %d assumptions."
                             % len(prems))
        self.facts = dict(zip(names, prems))
        self.assumptions = list(names)
        self.lemmas = []
        self.lines = []
        self.theorem = None

    def _lemma (self, name, formula, by, using):
        from .logic_parser import parse_formula
        if isinstance(formula, str):
            formula = parse_formula(formula, symbols=self.symbols)
        used = list(self.facts) if using is None else list(using)
        for u in used:
            if u not in self.facts:
                raise ValueError("No fact named %s." % u)
        hyps = [self.facts[u] for u in used]
        body = Implies(And(*hyps), formula) if hyps else formula
        vs = tuple(v for v in self.fixed if v in body.free_symbols)
        st = ForAll(Lambda(vs, body)) if vs else body
        th = Theorem("%s.%s" % (self.name, name), st, by)
        th.lemma_of = self.name
        th.verify()
        self.lemmas.append(th)
        self.lines.append((name, formula, used,
                           getattr(by, 'desc', '?')))
        return formula

    def have (self, name, formula, by, using=None):
        """Derive formula from the facts named in using (all the
        facts if None) by the tactic by, and name it."""
        if name in self.facts:
            raise ValueError("%s is already a fact." % name)
        self.facts[name] = self._lemma(name, formula, by, using)

    def show (self, by, using=None):
        """Derive the conclusion, and return the proved Theorem."""
        self._lemma("show", self.concl, by, using)
        lemmas = list(self.lemmas)
        def proof (z):
            return auto(z, _rules(lemmas),
                        depth=len(lemmas) + 2,
                        max_nodes=200 * (len(lemmas) + 1))
        proof = _tactic(proof, "structured proof (%d steps)"
                        % len(lemmas), lemmas)
        th = Theorem(self.name, self.statement, proof)
        th.verify()
        self.theorem = th
        return th

    def __str__ (self):
        l = ["proof %s: %s" % (self.name, formula_str(self.statement))]
        if self.fixed:
            l.append("  fix " + " ".join(str(v.name) for v in self.fixed))
        for a in self.assumptions:
            l.append("  assume %s: %s" % (a, formula_str(self.facts[a])))
        for name, f, used, desc in self.lines:
            head = "show" if name == "show" else "have %s:" % name
            l.append("  %s %s" % (head, formula_str(f)))
            l.append("      using %s  by %s" % (", ".join(used) or "-",
                                                desc))
        if self.theorem is not None:
            l.append("qed")
        return "\n".join(l)

    __repr__ = __str__


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
