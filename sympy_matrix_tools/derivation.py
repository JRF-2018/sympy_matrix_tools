# -*- coding: utf-8 -*-
__version__ = '0.6.0' # Time-stamp: <2026-10-22T16:54:53Z>

"""A ledger of derivations with numerical spot checks.

A Derivation records facts (hypotheses, definitions and derived
equations) and the steps of each derivation.  Each step says why it
is believed:

  * cas: a computation by SymPy (expand, doit, mat_diff, ...),
  * rule: a rewrite by a named fact (a hypothesis, a definition or a
    derived fact),
  * nothing: the step is only asserted.

Every step can also be checked numerically: the symbols are replaced
by random rational numbers and matrices that satisfy the hypotheses
(symmetric, positive definite, invertible, definitions), and both
sides are compared exactly.  The numerical check finds wrong steps
whoever wrote them, a human, an AI or a buggy CAS.

As in the calc of Lean, the intermediate expressions are written in
the script, so that the script can be read without running it.
"""

import random
from sympy import Basic, Eq, Symbol, MatrixSymbol, MatrixExpr, \
    ImmutableMatrix, Matrix, Rational, Integer, ZeroMatrix, Identity, \
    MatAdd, Q, diff, sympify
from sympy import latex as _latex
from sympy.assumptions.assume import AppliedPredicate
from .matrix_tools import mat_expand, mat_coeff
from .matrix_function import mat_diff


class Grad (MatrixExpr):
    """The gradient of a scalar (or 1x1 matrix) f by x, unevaluated.

    Grad(f, x).doit() is mat_diff(f, x).  Its shape is x.shape.  f
    must be written explicitly in x: a symbol defined by x (see
    Derivation.define) is treated as independent of x.
    """

    def __new__ (cls, f, x):
        return Basic.__new__(cls, sympify(f), x)

    f = property(lambda self: self.args[0])
    x = property(lambda self: self.args[1])
    shape = property(lambda self: self.x.shape)

    def doit (self, **hints):
        f = self.f.doit(**hints) if hints.get('deep', True) else self.f
        r = mat_diff(f, self.x)
        if r == 0:
            return ZeroMatrix(*self.shape)
        return r

    def _latex (self, printer):
        return r"\frac{\partial}{\partial %s}\left(%s\right)" \
            % (printer._print(self.x), printer._print(self.f))


# Justifications.

class Justification:
    def apply (self, expr, ledger):
        """Return candidates of the result of applying this to expr."""
        raise NotImplementedError()


class cas (Justification):
    """A computation by SymPy.

    func is a function or a method name such as "doit" or "expand",
    called as func(expr, *args).
    """

    def __init__ (self, func, *args, label=None):
        self.func = func
        self.args = args
        if label is None:
            if isinstance(func, str):
                label = func
            else:
                label = getattr(func, '__name__', str(func))
            if args:
                label += "(%s)" % ", ".join(str(a) for a in args)
        self.label = label

    def apply (self, expr, ledger):
        if isinstance(self.func, str):
            return [getattr(expr, self.func)(*self.args)]
        return [self.func(expr, *self.args)]

    def __str__ (self):
        return "cas " + self.label

    deps = property(lambda self: ())


class rule (Justification):
    """A rewrite by the named fact lhs = rhs (in both directions)."""

    def __init__ (self, name):
        self.name = name

    def apply (self, expr, ledger):
        prop = ledger.facts[self.name].prop
        if not isinstance(prop, Eq):
            raise ValueError("%s is not an equation." % self.name)
        l, r = prop.args
        pairs = [(l, r), (r, l)]
        if getattr(l, 'is_Matrix', False):
            # The transposes are distributed by SymPy, e.g.
            # (A*y0 + B*x).T is x.T*B.T + y0.T*A.T.
            lt, rt = l.T.doit(), r.T.doit()
            pairs += [(lt, rt), (rt, lt)]
        cands = []
        for e in (expr, mat_expand(expr)):
            # Rewrite both lhs and its transpose in one direction.
            for group in (pairs[0::2], pairs[1::2]):
                c = e
                for a, b in group:
                    c = c.subs(a, b)
                if c not in cands:
                    cands.append(c)
            for a, b in pairs:
                for c in (e.subs(a, b), e.xreplace({a: b})):
                    if c not in cands:
                        cands.append(c)
        return cands

    def __str__ (self):
        return "rule " + self.name

    deps = property(lambda self: (self.name,))


def _is_zero (e):
    if e == 0:
        return True
    if isinstance(e, ZeroMatrix):
        return True
    if isinstance(e, (Matrix, ImmutableMatrix)):
        return e.is_zero_matrix
    return False


def symbolic_equal (a, b):
    """Check a == b symbolically (incompletely)."""
    if a == b:
        return True
    try:
        if mat_expand(a) == mat_expand(b):
            return True
        return _is_zero(mat_expand((a - b).doit()))
    except Exception:
        return False


# Facts and steps.

class Step:
    def __init__ (self, expr, justifications):
        self.expr = sympify(expr)
        self.justifications = list(justifications)
        self.symbolic = None      # True, False or None (not checked)
        self.numeric = None       # (passed, trials), "n/a" or None
        self.counterexample = None

    @property
    def deps (self):
        return sum([list(j.deps) for j in self.justifications], [])

    def label (self):
        if not self.justifications:
            return "asserted"
        return "; ".join(str(j) for j in self.justifications)


class Fact:
    KINDS = ("hypothesis", "definition", "derived")

    def __init__ (self, name, prop, kind, steps=(), deps=(),
                  obligations=(), note=None):
        self.name = name
        self.prop = prop
        self.kind = kind
        self.steps = list(steps)
        self.deps = list(deps)
        self.obligations = list(obligations)   # [(prop, discharged_by)]
        self.note = note
        self.numeric = None


def _numeric_equal (a, b):
    d = a - b
    if hasattr(d, 'applyfunc'):
        return _is_zero(d.applyfunc(lambda e: e.simplify()))
    return sympify(d).simplify() == 0


def _unit (e):
    # (1, 1) matrices are compared as scalars.
    if isinstance(e, (Matrix, ImmutableMatrix)) and e.shape == (1, 1):
        return e[0, 0]
    return e


class Derivation:
    """A ledger of a derivation.

    >>> from sympy import symbols, MatrixSymbol, Eq, Q
    >>> n = symbols("n", integer=True, positive=True)
    >>> A = MatrixSymbol("A", n, n)
    >>> X = MatrixSymbol("X", n, n)
    >>> d = Derivation("small example")
    >>> d.hypothesis(Eq(A.T, A), "A_sym")
    >>> d.define(X, A * A.T, "X_def")
    >>> d.calc("X_sq", X,
    ...        (A * A.T, rule("X_def")),
    ...        (A * A, rule("A_sym")))
    >>> d.check(trials=2)
    True
    >>> print(d.report())
    Derivation: small example
    [hypothesis] A_sym: A.T = A
    [definition] X_def: X = A*A.T
    [derived] X_sq: X = A**2
      X
      = A*A.T    by rule X_def    symbolic: ok    numeric: ok 2/2
      = A**2    by rule A_sym    symbolic: ok    numeric: ok 2/2
    Summary: 2 steps, 0 asserted, 0 symbolic failures, 0 numeric failures, 0 open obligations.
    """

    def __init__ (self, title=""):
        self.title = title
        self.facts = {}
        self.order = []
        self.aliases = {}

    def alias (self, expr, name):
        """Show expr as name in reports (only for display).  The
        transpose of expr is shown as name.T."""
        expr = sympify(expr)
        if getattr(expr, 'is_Matrix', False):
            s = MatrixSymbol(name, *expr.shape)
            self.aliases[expr] = s
            self.aliases[expr.T.doit()] = s.T
        else:
            self.aliases[expr] = Symbol(name)

    def _show (self, e, printer=str):
        e = sympify(e)
        if isinstance(e, Eq):
            # Don't rebuild Eq, which evaluates its arguments.
            return "%s = %s" % (self._show(e.lhs, printer),
                                self._show(e.rhs, printer))
        for k in sorted(self.aliases, key=lambda k: -len(str(k))):
            e = e.xreplace({k: self.aliases[k]})
        return printer(e)

    def _add (self, fact):
        if fact.name in self.facts:
            raise ValueError("%s is already used." % fact.name)
        self.facts[fact.name] = fact
        self.order.append(fact.name)
        return fact

    def hypothesis (self, prop, name):
        """Add a hypothesis.  Eq(X.T, X) and Q.symmetric(X),
        Q.positive_definite(X) and Q.invertible(...) are respected in
        the numerical checks.  Other hypotheses are constraints that
        random values do not satisfy."""
        self._add(Fact(name, sympify(prop), "hypothesis"))

    def define (self, symbol, expr, name=None):
        """Define symbol as an abbreviation of expr."""
        if name is None:
            name = str(symbol) + "_def"
        self._add(Fact(name, Eq(symbol, expr, evaluate=False),
                       "definition"))

    def calc (self, name, start, *steps):
        """Add a derived fact Eq(start, last) by a chain of equalities.

        Each step is (expr, justification, ...) or expr (asserted).
        The justifications are applied in order to the previous
        expression, and one of the results must be equal to expr.
        """
        l = []
        for s in steps:
            if isinstance(s, tuple):
                l.append(Step(s[0], s[1:]))
            else:
                l.append(Step(s, ()))
        deps = sorted(set(sum([st.deps for st in l], [])))
        prop = Eq(start, l[-1].expr, evaluate=False)
        fact = self._add(Fact(name, prop, "derived", steps=l, deps=deps))
        fact.start = sympify(start)
        self._check_symbolic(fact)

    def solve_linear (self, name, fact_name, x):
        """Derive x = -K**-1 * c from the fact K*x + c = 0.

        The invertibility of K is an obligation.  It is discharged by
        a hypothesis Q.invertible(K), where K may be written by a
        defined symbol.  If K is a defined expression, the defined
        symbol is used in the solution.
        """
        prop = self.facts[fact_name].prop
        e = prop.lhs if _is_zero(prop.rhs) else (prop.lhs - prop.rhs)
        for e in (e, mat_expand(e.doit())):
            # Try the form written first, and then the expanded form.
            terms = e.args if isinstance(e, MatAdd) else (e,)
            K = []
            c = []
            for t in terms:
                if t.has(x):
                    k = mat_coeff(t, x)
                    if k is None or k.has(x):
                        K = None
                        break
                    K.append(Identity(x.rows) if k == 1 else k)
                else:
                    c.append(t)
            if K:
                break
        if not K:
            raise ValueError("%s is not linear in %s." % (prop, x))
        K = MatAdd(*K).doit() if len(K) > 1 else K[0]
        c = MatAdd(*c).doit() if len(c) > 1 else \
            (c[0] if c else ZeroMatrix(*x.shape))
        Ksym = K
        for f in self.facts.values():
            if f.kind == "definition" and symbolic_equal(f.prop.rhs, K):
                Ksym = f.prop.lhs
        obligation = Q.invertible(Ksym)
        by = None
        for f in self.facts.values():
            if f.kind == "hypothesis" and f.prop == obligation:
                by = f.name
            elif f.kind == "hypothesis" \
                 and isinstance(f.prop, AppliedPredicate) \
                 and f.prop.function == Q.invertible \
                 and symbolic_equal(f.prop.arguments[0], K):
                by = f.name
        sol = -Ksym.inverse() * c
        fact = self._add(Fact(name, Eq(x, sol, evaluate=False), "derived",
                              deps=[fact_name] + ([by] if by else []),
                              obligations=[(obligation, by)],
                              note="solve_linear(%s, %s)" % (fact_name, x)))
        fact.solved = (fact_name, x)

    # Symbolic checks.

    def _check_symbolic (self, fact):
        prev = fact.start
        for st in fact.steps:
            if not st.justifications:
                st.symbolic = None
            else:
                cands = [prev]
                try:
                    for j in st.justifications:
                        cands = sum([j.apply(c, self) for c in cands], [])
                    st.symbolic = any(symbolic_equal(c, st.expr)
                                      for c in cands)
                except Exception:
                    st.symbolic = False
            prev = st.expr

    # Numerical checks.

    def _symbols (self):
        s = set()
        for f in self.facts.values():
            s |= f.prop.atoms(MatrixSymbol, Symbol)
            for st in f.steps:
                s |= st.expr.atoms(MatrixSymbol, Symbol)
            if hasattr(f, 'start'):
                s |= f.start.atoms(MatrixSymbol, Symbol)
        return s

    def _generatable (self, name, memo=None):
        if memo is None:
            memo = {}
        if name in memo:
            return memo[name]
        f = self.facts[name]
        if f.kind == "definition":
            r = True
        elif f.kind == "hypothesis":
            r = self._is_respected(f.prop)
        else:
            # A fact by solve_linear constrains the variable, so it is
            # not generatable either (it depends on the constraint).
            r = all(self._generatable(d, memo) for d in f.deps)
        memo[name] = r
        return r

    def _is_respected (self, prop):
        if isinstance(prop, Eq):
            # X.T = X or X = X.T
            return any(isinstance(a, MatrixSymbol) and b == a.T
                       for a, b in (prop.args, prop.args[::-1]))
        if isinstance(prop, AppliedPredicate):
            return prop.function in (Q.symmetric, Q.positive_definite,
                                     Q.invertible)
        return False

    def _constraints (self):
        sym, pd, inv = set(), set(), []
        for f in self.facts.values():
            if f.kind != "hypothesis":
                continue
            p = f.prop
            if isinstance(p, Eq):
                for a, b in (p.args, p.args[::-1]):
                    if isinstance(a, MatrixSymbol) and b == a.T:
                        sym.add(a)
            elif isinstance(p, AppliedPredicate):
                if p.function == Q.symmetric:
                    sym.add(p.arguments[0])
                elif p.function == Q.positive_definite:
                    pd.add(p.arguments[0])
                elif p.function == Q.invertible:
                    inv.append(p.arguments[0])
        return sym, pd, inv

    def _random_rational (self, rng):
        return Rational(rng.randint(-9, 9), rng.randint(1, 3))

    def instance (self, rng):
        """Return random values of the symbols that satisfy the
        hypotheses respected and the definitions."""
        syms = self._symbols()
        defs = {f.prop.lhs: f.prop.rhs for f in self.facts.values()
                if f.kind == "definition"}
        sym, pd, inv = self._constraints()
        dims = set()
        for s in syms:
            if isinstance(s, MatrixSymbol):
                for d in s.shape:
                    dims |= d.free_symbols
        values = {}
        for d in sorted(dims, key=str):
            values[d] = Integer(rng.randint(2, 4))
        for s in sorted(syms, key=str):
            if s in dims or s in defs or s in values:
                continue
            if isinstance(s, MatrixSymbol):
                r, c = [int(sympify(k).xreplace(values)) for k in s.shape]
                M = Matrix(r, c, lambda i, j: self._random_rational(rng))
                if s in pd:
                    M = M * M.T + Matrix.eye(r)
                elif s in sym:
                    M = (M + M.T) / 2
                values[s] = ImmutableMatrix(M)
            elif isinstance(s, Symbol) and not s.is_integer:
                values[s] = self._random_rational(rng)
        pending = dict(defs)
        while pending:
            done = False
            for k, v in list(pending.items()):
                if not (v.free_symbols & set(pending)) \
                   and not any(v.has(p) for p in pending):
                    values[k] = self.evaluate(v, values)
                    del pending[k]
                    done = True
            if not done:
                raise ValueError("Circular definitions.")
        for k in inv:
            if self.evaluate(k, values).det() == 0:
                return None
        return values

    def evaluate (self, expr, values):
        """Evaluate expr exactly with values (see instance)."""
        expr = sympify(expr)
        for g in sorted(expr.atoms(Grad), key=lambda g: -len(str(g))):
            expr = expr.xreplace({g: self._evaluate_grad(g, values)})
        r = expr.xreplace(values).doit()
        if isinstance(r, MatrixExpr):
            r = r.as_explicit()
        if isinstance(r, (Matrix, ImmutableMatrix)):
            r = ImmutableMatrix(r).applyfunc(lambda a: a.doit())
        return _unit(r)

    def _evaluate_grad (self, g, values):
        x = g.x
        v = dict(values)
        xv = v.pop(x)
        r, c = xv.shape
        X = ImmutableMatrix(r, c, lambda i, j: Symbol("_g_%d_%d" % (i, j)))
        v[x] = X
        f = self.evaluate(g.f, v)
        G = ImmutableMatrix(r, c, lambda i, j: diff(f, X[i, j]))
        return G.xreplace(dict(zip(X, xv)))

    def check (self, trials=3, seed=0):
        """Check every step and fact numerically.  Return True if no
        failure is found."""
        rng = random.Random(seed)
        instances = []
        tries = 0
        while len(instances) < trials and tries < trials * 20:
            tries += 1
            v = self.instance(rng)
            if v is not None:
                instances.append(v)
        ok = True
        for name in self.order:
            f = self.facts[name]
            if f.kind != "derived":
                continue
            if hasattr(f, 'solved'):
                f.numeric = self._check_solved(f, instances)
                ok = ok and f.numeric[0] == f.numeric[1]
                continue
            prev = f.start
            for st in f.steps:
                if not all(self._generatable(d) for d in st.deps):
                    st.numeric = "n/a"
                else:
                    n = 0
                    for v in instances:
                        try:
                            same = _numeric_equal(self.evaluate(prev, v),
                                                  self.evaluate(st.expr, v))
                        except Exception:
                            same = False
                        if same:
                            n += 1
                        elif st.counterexample is None:
                            st.counterexample = v
                    st.numeric = (n, len(instances))
                    ok = ok and n == len(instances)
                prev = st.expr
        return ok

    def _check_solved (self, f, instances):
        fact_name, x = f.solved
        prop = self.facts[fact_name].prop
        n = 0
        for v in instances:
            w = dict(v)
            w[x] = ImmutableMatrix(self.evaluate(f.prop.rhs, v))
            try:
                d = self.evaluate(prop.lhs - prop.rhs, w)
                same = _is_zero(d)
            except Exception:
                same = False
            if same:
                n += 1
        return (n, len(instances))

    # Reports.

    def _numeric_str (self, num):
        if num is None:
            return "not checked"
        if num == "n/a":
            return "n/a"
        n, t = num
        return ("ok %d/%d" if n == t else "FAILED %d/%d") % (n, t)

    def _symbolic_str (self, s):
        return {True: "ok", False: "FAILED", None: "none"}[s]

    def summary (self):
        steps = [st for f in self.facts.values() for st in f.steps]
        asserted = sum(1 for st in steps if not st.justifications)
        sfail = sum(1 for st in steps if st.symbolic is False)
        nfail = sum(1 for st in steps if isinstance(st.numeric, tuple)
                    and st.numeric[0] != st.numeric[1])
        nfail += sum(1 for f in self.facts.values()
                     if isinstance(f.numeric, tuple)
                     and f.numeric[0] != f.numeric[1])
        opened = sum(1 for f in self.facts.values()
                     for p, by in f.obligations if by is None)
        return dict(steps=len(steps), asserted=asserted,
                    symbolic_failures=sfail, numeric_failures=nfail,
                    open_obligations=opened)

    def report (self, printer=str):
        """Return the report as a plain text."""
        _p = printer
        printer = lambda e: self._show(e, _p)
        lines = ["Derivation: %s" % self.title]
        for name in self.order:
            f = self.facts[name]
            lines.append("[%s] %s: %s" % (f.kind, name, printer(f.prop)))
            if f.steps:
                lines.append("  %s" % printer(f.start))
                for st in f.steps:
                    lines.append("  = %s    by %s    symbolic: %s"
                                 "    numeric: %s"
                                 % (printer(st.expr), st.label(),
                                    self._symbolic_str(st.symbolic),
                                    self._numeric_str(st.numeric)))
            if f.note:
                lines.append("  by %s    numeric: %s"
                             % (f.note, self._numeric_str(f.numeric)))
            for p, by in f.obligations:
                lines.append("  obligation %s: %s"
                             % (printer(p), ("discharged by " + by)
                                if by else "OPEN"))
        s = self.summary()
        lines.append("Summary: %(steps)d steps, %(asserted)d asserted, "
                     "%(symbolic_failures)d symbolic failures, "
                     "%(numeric_failures)d numeric failures, "
                     "%(open_obligations)d open obligations." % s)
        return "\n".join(lines)

    def _repr_latex_ (self):
        latex = lambda e: self._show(e, _latex)
        rows = [r"\textbf{%s}" % self.title.replace("_", r"\_")]
        for name in self.order:
            f = self.facts[name]
            nm = name.replace("_", r"\_")
            rows.append(r"\text{[%s] %s:}\ %s" % (f.kind, nm, latex(f.prop)))
            if f.steps:
                rows.append(r"\quad %s" % latex(f.start))
                for st in f.steps:
                    lab = st.label().replace("_", r"\_")
                    rows.append(r"\quad = %s & \text{%s; %s, %s}"
                                % (latex(st.expr), lab,
                                   self._symbolic_str(st.symbolic),
                                   self._numeric_str(st.numeric)))
            for p, by in f.obligations:
                rows.append(r"\quad\text{obligation } %s: \text{%s}"
                            % (latex(p), ("discharged by " + by)
                               .replace("_", r"\_") if by else "OPEN"))
        return r"$\begin{array}{ll}" + r"\\".join(rows) \
            + r"\end{array}$"
