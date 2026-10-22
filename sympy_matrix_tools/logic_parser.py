# -*- coding: utf-8 -*-
"""Parse formulas written in the logical notation.

parse_formula reads the notation that formula_str and
print_proofstate(..., pretty=True) print, e.g.
"∀x. P(x) → (∃y. P(y))", and returns the formula of logic_tools and
fol (ForAll, Exists, Implies, And, Or, Bot, PredApply, Wild).

See docs/fol.md.
"""
__version__ = '0.8.0' # Time-stamp: <2026-10-22T18:27:21Z>

import re
from sympy import Symbol, Wild, Lambda, Implies, And, Or, Equivalent, \
    sympify
from sympy.parsing.sympy_parser import parse_expr
from .apply import Apply, PredApply
from .logic_tools import ForAll, FrozenGoal
from .fol import Exists, Bot, neg


_TOKENS = [
    ('FORALL', r'∀|\bALL\b|\bforall\b'),
    ('EXISTS', r'∃|\bEX\b|\bexists\b'),
    ('IFF', r'↔|<->'),
    ('IMP', r'→|->'),
    ('AND', r'∧|&'),
    ('OR', r'∨|\|'),
    ('NOT', r'¬|~'),
    ('FROZEN', r'\*'),
    ('BOT', r'⊥|\bBot\b'),
    ('LPAREN', r'\('),
    ('RPAREN', r'\)'),
    ('DOT', r'\.'),
    ('COMMA', r','),
    ('META', r'\?[A-Za-z_][A-Za-z0-9_]*'),
    ('NAME', r'[A-Za-z_][A-Za-z0-9_]*'),
    ('LBRACKET', r'\['),
]
_TOKEN_RE = re.compile('|'.join('(?P<%s>%s)' % t for t in _TOKENS))


class FormulaSyntaxError (ValueError):
    pass


class _Parser:
    """A recursive-descent parser into a small tree of tuples."""

    def __init__ (self, text):
        self.text = text
        self.pos = 0

    def error (self, msg):
        raise FormulaSyntaxError("%s at %d: %r" % (msg, self.pos,
                                                  self.text[self.pos:]))

    def skip (self):
        while self.pos < len(self.text) and self.text[self.pos].isspace():
            self.pos += 1

    def peek (self):
        self.skip()
        m = _TOKEN_RE.match(self.text, self.pos)
        if m is None:
            return None, None
        return m.lastgroup, m.group()

    def take (self, kind=None):
        k, v = self.peek()
        if k is None or (kind is not None and k != kind):
            self.error("Expected %s" % (kind or "a token"))
        self.pos += len(v)
        return v

    def raw_until_close (self, open_c, close_c):
        # The position is just after open_c.  Return the text until the
        # matching close_c, split at the commas of the top level.
        depth = 1
        start = self.pos
        args = []
        while self.pos < len(self.text):
            c = self.text[self.pos]
            if c in "([{":
                depth += 1
            elif c in ")]}":
                depth -= 1
                if depth == 0:
                    args.append(self.text[start:self.pos])
                    self.pos += 1
                    return args
            elif c == "," and depth == 1 and close_c == ")":
                args.append(self.text[start:self.pos])
                start = self.pos + 1
            self.pos += 1
        self.error("Unbalanced %s" % open_c)

    def parse (self):
        r = self.formula()
        self.skip()
        if self.pos != len(self.text):
            self.error("Unexpected text")
        return r

    def formula (self):
        k, v = self.peek()
        if k in ('FORALL', 'EXISTS'):
            self.take()
            vs = []
            while True:
                k2, v2 = self.peek()
                if k2 == 'NAME':
                    vs.append(self.take())
                elif k2 == 'COMMA':
                    self.take()
                elif k2 == 'DOT':
                    self.take()
                    break
                else:
                    self.error("Expected a variable or '.'")
            if not vs:
                self.error("No variable")
            return ('forall' if k == 'FORALL' else 'exists', vs,
                    self.formula())
        return self.iff()

    def iff (self):
        a = self.imp()
        if self.peek()[0] == 'IFF':
            self.take()
            return ('iff', a, self.imp())
        return a

    def imp (self):
        a = self.disj()
        if self.peek()[0] == 'IMP':
            self.take()
            return ('imp', a, self.imp_rhs())
        return a

    def imp_rhs (self):
        # A binder may follow an operator: "A → ∀x. P(x)".
        if self.peek()[0] in ('FORALL', 'EXISTS'):
            return self.formula()
        return self.imp()

    def disj (self):
        l = [self.conj()]
        while self.peek()[0] == 'OR':
            self.take()
            l.append(self.conj())
        return l[0] if len(l) == 1 else ('or', l)

    def conj (self):
        l = [self.unary()]
        while self.peek()[0] == 'AND':
            self.take()
            l.append(self.unary())
        return l[0] if len(l) == 1 else ('and', l)

    def unary (self):
        k, v = self.peek()
        if k == 'NOT':
            self.take()
            return ('not', self.unary())
        if k == 'FROZEN':
            self.take()
            return ('frozen', self.unary())
        if k in ('FORALL', 'EXISTS'):
            return self.formula()
        return self.atom()

    def atom (self):
        k, v = self.peek()
        if k == 'LPAREN':
            self.take()
            r = self.formula()
            self.take('RPAREN')
            return r
        if k == 'BOT':
            self.take()
            return ('bot',)
        if k == 'LBRACKET':
            self.take()
            return ('sympy', self.raw_until_close('[', ']')[0])
        if k in ('NAME', 'META'):
            self.take()
            meta = k == 'META'
            name = v[1:] if meta else v
            if self.text.startswith('(', self.pos):
                self.pos += 1
                return ('pred', name, meta, self.raw_until_close('(', ')'))
            return ('prop', name, meta)
        self.error("Expected a formula")


def _props (node):
    """The names used as propositions (not as predicates or terms)."""
    t = node[0]
    if t == 'prop':
        return {node[1]}
    if t in ('forall', 'exists'):
        return _props(node[2]) - set(node[1])
    if t in ('imp', 'iff'):
        return _props(node[1]) | _props(node[2])
    if t in ('and', 'or'):
        return set().union(*[_props(a) for a in node[1]])
    if t in ('not', 'frozen'):
        return _props(node[1])
    return set()


class _Builder:
    def __init__ (self, symbols, predicates):
        self.symbols = dict(symbols or {})
        self.predicates = dict(predicates or {})
        self.wilds = {}

    def wild (self, name, boolean):
        key = (name, boolean)
        if key not in self.wilds:
            self.wilds[key] = Wild(name, boolean=True) if boolean \
                else Wild(name)
        return self.wilds[key]

    def term (self, text, env):
        # ?F(a, b) is Apply(Wild("F"), a, b), and ?x is Wild("x").
        local = {}
        local.update(self.symbols)
        local.update(env)

        def meta (m):
            name = "__meta_" + m.group(1)
            local[name] = self.wild(m.group(1), False)
            if m.group(2):
                return "__apply(%s, " % name
            return name
        text = re.sub(r'\?([A-Za-z_][A-Za-z0-9_]*)(\()?', meta, text)
        local['__apply'] = lambda f, *a: Apply(f, *a)
        return parse_expr(text, local_dict=local)

    def build (self, node, env):
        t = node[0]
        if t in ('forall', 'exists'):
            props = _props(node[2])
            vs = [Symbol(v, boolean=True) if v in props else Symbol(v)
                  for v in node[1]]
            env2 = dict(env)
            env2.update(zip(node[1], vs))
            body = self.build(node[2], env2)
            if t == 'forall':
                return ForAll(Lambda(tuple(vs), body))
            for v in reversed(vs):
                body = Exists(Lambda(v, body))
            return body
        if t == 'imp':
            return Implies(self.build(node[1], env), self.build(node[2], env))
        if t == 'iff':
            return Equivalent(self.build(node[1], env),
                              self.build(node[2], env))
        if t == 'and':
            return And(*[self.build(a, env) for a in node[1]])
        if t == 'or':
            return Or(*[self.build(a, env) for a in node[1]])
        if t == 'not':
            return neg(self.build(node[1], env))
        if t == 'frozen':
            return FrozenGoal(self.build(node[1], env))
        if t == 'bot':
            return Bot
        if t == 'sympy':
            local = dict(self.symbols)
            local.update(env)
            return sympify(parse_expr(node[1], local_dict=local))
        if t == 'prop':
            name, meta = node[1], node[2]
            if meta:
                return self.wild(name, True)
            if name in env:
                return env[name]
            if name in self.symbols:
                return self.symbols[name]
            return Symbol(name, boolean=True)
        if t == 'pred':
            name, meta, args = node[1], node[2], node[3]
            args = [self.term(a, env) for a in args if a.strip()]
            if meta:
                return PredApply(self.wild(name, False), *args)
            if name not in env and name in self.predicates:
                return self.predicates[name](*args)
            head = env.get(name) or self.symbols.get(name) or Symbol(name)
            return PredApply(head, *args)
        raise ValueError(node)


def parse_formula (text, symbols=None, predicates=None):
    """Parse a formula written in the logical notation.

    The notation is that of formula_str, with ASCII alternatives:

      * ∀x y. φ (ALL x y. φ), ∃x. φ (EX x. φ): binders extend as far
        to the right as possible;
      * φ ↔ ψ (<->), φ → ψ (->, right associative), φ ∨ ψ (|),
        φ ∧ ψ (&), ¬φ (~), ⊥ (Bot), in this order of precedence;
      * P(t, ...) is PredApply(P, t, ...), where the arguments are
        SymPy expressions; a bare name A is a proposition
        Symbol("A", boolean=True);
      * ?x is a metavariable (Wild), ?P(x) is PredApply(Wild("P"), x),
        and ?f(x) in a term is Apply(Wild("f"), x);
      * [expr] is any SymPy expression, e.g. [x < 3];
      * *φ is FrozenGoal(φ), the goal of a proofstate.

    ¬φ is neg(φ) = Implies(φ, Bot), and ∃ is fol.Exists.  A bound
    variable is a Symbol, with boolean=True if it is used as a
    proposition.  symbols maps names to the objects to use for the free
    names, and predicates maps names to functions (e.g. SymPy
    Predicate's) to apply instead of PredApply.

    >>> from sympy_matrix_tools import formula_str
    >>> f = parse_formula("ALL x. P(x) -> (EX y. Q(f(x), y))")
    >>> f
    Q.ForAll(Lambda(x, Implies(PredApply(P, x), Q.Exists(Lambda(y, PredApply(Q, f(x), y))))))
    >>> formula_str(f)
    '∀x. P(x) → (∃y. Q(f(x), y))'
    >>> parse_formula(formula_str(f)) == f
    True
    >>> parse_formula("∀A. ¬¬A → A")
    Q.ForAll(Lambda(A, Implies(Implies(Implies(A, Bot), Bot), A)))
    """
    tree = _Parser(text).parse()
    return _Builder(symbols, predicates).build(tree, {})
