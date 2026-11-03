# -*- coding: utf-8 -*-
__version__ = '0.10.0' # Time-stamp: <2026-11-03T03:41:31Z>

"""Logic with matrices: parse_formula with matrix binders and
equations, rewrite_goal, cas_goal, Theorem.instance, Theorem.cas,
ledger theorems and structured proofs (Proof)."""

import pytest
from sympy import Symbol, MatrixSymbol, Lambda, Eq, Implies, And, \
    ZeroMatrix, Function
from sympy_matrix_tools import PredApply, ForAll, make_proofstate, \
    rewrite_goal, cas_goal, formula_str, parse_formula, \
    get_syms_prems_concl
from sympy_matrix_tools.fol import Theorem, Proof, by_auto, by_rule, \
    by_rewrite, by_steps, by_cas
from sympy_matrix_tools.derivation import Derivation, cas, rule
from sympy_matrix_tools.matrix_tools import mat_expand

n = Symbol("n", integer=True, positive=True)
S = {"n": n}
S.update({k: MatrixSymbol(k, n, n) for k in "A B C M N P R W".split()})
S.update({k: MatrixSymbol(k, n, 1) for k in ["x", "y"]})
S.update({k: MatrixSymbol(k, 1, 1) for k in ["a", "b"]})
A, B, C, M, N, P, R, W, x = [S[k] for k in "A B C M N P R W x".split()]


def axioms ():
    return dict(
        nonneg_add=Theorem("nonneg_add",
                           "∀a b. nonneg(a) ∧ nonneg(b) → nonneg(a + b)",
                           axiom=True, symbols=S),
        psd_def=Theorem("psd_def", "∀M. psd(M) ↔ (∀x. nonneg(x.T*M*x))",
                        axiom=True, symbols=S),
        le_def=Theorem("le_def", "∀A B. mat_le(A, B) ↔ psd(B - A)",
                       axiom=True, symbols=S))


def test_parse_matrix_binders ():
    f = parse_formula("∀M x. nonneg(x.T*M*x)", symbols=S)
    lam = f.arguments[0]
    assert lam.signature == (M, x)
    assert parse_formula(formula_str(f), symbols=S) == f
    # Without symbols, the bound variables are scalar symbols.
    g = parse_formula("∀M. P(M)")
    assert g.arguments[0].signature == (Symbol("M"),)


def test_parse_equation ():
    x1, y1 = Symbol("x"), Symbol("y")
    f = Function("f")
    e = parse_formula("∀x. f(x) = 2*x ∧ P(x)")
    body = e.arguments[0].expr
    assert Eq(f(x1), 2 * x1) in body.args
    assert formula_str(e) == "∀x. f(x) = 2*x ∧ P(x)" \
        or formula_str(e) == "∀x. P(x) ∧ f(x) = 2*x"
    assert parse_formula(formula_str(e)) == e
    # Comparisons are not equations.
    assert parse_formula("[x <= 3] → P(x)") \
        == Implies(x1 <= 3, PredApply(Symbol("P"), x1))
    assert parse_formula("x = y → (y = x)") \
        == Implies(Eq(x1, y1), Eq(y1, x1))
    q = parse_formula("∀M N x. x.T*(M + N)*x = x.T*M*x + x.T*N*x",
                      symbols=S)
    assert get_syms_prems_concl(q)[2] \
        == Eq(x.T*(M + N)*x, x.T*M*x + x.T*N*x)


def test_ac_unification ():
    from sympy import Wild
    from sympy_matrix_tools import ho_unify, MatrixWild
    Mw, Nw = MatrixWild("M", n, n), MatrixWild("N", n, n)
    rs = list(ho_unify(Mw + Nw, A + B + C))
    assert len(rs) == 6
    assert all((r[Mw] + r[Nw]).doit() == A + B + C for r in rs)
    rs = list(ho_unify(P.T*Mw*P, P.T*A*B*P))
    assert rs == [{Mw: A*B}]
    assert list(ho_unify(P.T*Mw*P, P.T*A*B)) == []
    u, v, w = Symbol("u"), Symbol("v"), Symbol("w")
    U = Wild("U")
    assert len(list(ho_unify(U + w, u + v + w))) == 1


def test_rewrite_and_cas_goal ():
    x1, y1 = Symbol("x"), Symbol("y")
    f = Function("f")
    Ps = Symbol("P")
    double = ForAll(Lambda(x1, Eq(f(x1), 2 * x1)))
    z = make_proofstate(ForAll(Lambda(y1, PredApply(Ps, f(y1) + 1))))
    r = rewrite_goal(z, double)
    assert formula_str(get_syms_prems_concl(r)[1][0]) == "∀y. P(2*y + 1)"
    r2 = rewrite_goal(r, double, reverse=True)
    assert formula_str(get_syms_prems_concl(r2)[1][0]) == "∀y. P(f(y) + 1)"
    assert rewrite_goal(z, ForAll(Lambda(x1, Eq(f(x1 + 1), x1)))) is None
    z = make_proofstate(PredApply(Ps, (x1 + 1)**2 - x1**2))
    r = cas_goal(z, "expand")
    assert get_syms_prems_concl(r)[1][0] == PredApply(Ps, 2 * x1 + 1)
    assert cas_goal(r, "expand") is None


def test_instance ():
    ax = axioms()
    psd_unfold = Theorem("psd_unfold", "∀M x. psd(M) → nonneg(x.T*M*x)",
                         by_auto(ax["psd_def"]), symbols=S)
    i = psd_unfold.instance(x=P*x)
    assert i.kind == "instance"
    assert formula_str(i.statement) \
        == "∀M P x. psd(M) → nonneg(x.T*P.T*M*P*x)"
    with pytest.raises(ValueError):
        psd_unfold.instance(x=P)            # the shapes differ
    with pytest.raises(ValueError):
        psd_unfold.instance(z=x)            # no such variable
    # The inner bound variables are renamed away from the new terms.
    u, v = Symbol("u"), Symbol("v")
    Qs, Rs = Symbol("Q"), Symbol("R")
    th = Theorem("t", ForAll(Lambda(u, Implies(
        ForAll(Lambda(v, PredApply(Qs, u, v))), PredApply(Rs, u)))),
                 axiom=True)
    s = formula_str(th.instance(u=v).statement)
    assert s.startswith("∀v. (∀") and "Q(v, v)" not in s
    assert th.instance(u=v).verify()
    from sympy_matrix_tools.fol import ccontr
    with pytest.raises(ValueError):
        ccontr.instance(A=u)                # A is a formula


def test_cas_theorem ():
    x1, y1 = Symbol("x"), Symbol("y")
    Ps = Symbol("P")
    th = Theorem("all_P", ForAll(Lambda(x1, PredApply(Ps, x1))), axiom=True)
    t = th.instance(x=(y1 + 1)**2 - y1**2).cas("expand")
    assert formula_str(t.statement) == "∀y. P(2*y + 1)"
    assert th in t.dependencies()


def test_ledger_theorem ():
    d = Derivation("ledger")
    d.calc("quad_add", x.T*(M + N)*x, (x.T*M*x + x.T*N*x, cas(mat_expand)))
    d.hypothesis(Eq(W.T, W), "W_sym")
    d.calc("sym", x.T*W.T*x, (x.T*W*x, rule("W_sym")))
    d.calc("bad", P.T*(B - A)*P, (P.T*B*P - P*A*P.T, cas(mat_expand)))
    d.calc("asserted", P.T*(B - A)*P, P.T*B*P - P.T*A*P)
    t = d.theorem("quad_add")
    assert t.kind == "ledger"
    assert formula_str(t.statement) \
        == "∀M N x. x.T*(M + N)*x = x.T*M*x + x.T*N*x"
    assert t.verify()
    assert "checked symbolically and numerically" in t.explain()
    # A hypothesis becomes a premise.
    t2 = d.theorem("sym")
    assert formula_str(t2.statement) == "∀W x. W.T = W → x.T*W.T*x = x.T*W*x"
    assert t2.verify()
    with pytest.raises(ValueError, match="numerical"):
        d.theorem("bad").verify()
    with pytest.raises(ValueError, match="asserted"):
        d.theorem("asserted").verify()
    with pytest.raises(ValueError):
        d.theorem("W_sym")


def test_backward_and_forward ():
    ax = axioms()
    d = Derivation("quadratic forms")
    d.calc("quad_add", x.T*(M + N)*x, (x.T*M*x + x.T*N*x, cas(mat_expand)))
    d.calc("congr_diff", P.T*(B - A)*P,
           (P.T*B*P - P.T*A*P, cas(mat_expand)))
    quad_add = d.theorem("quad_add")
    congr_diff = d.theorem("congr_diff")
    psd_unfold = Theorem("psd_unfold", "∀M x. psd(M) → nonneg(x.T*M*x)",
                         by_auto(ax["psd_def"]), symbols=S)
    psd_add = Theorem("psd_add", "∀M N. psd(M) ∧ psd(N) → psd(M + N)",
                      by_steps(by_rule(ax["psd_def"]), by_rewrite(quad_add),
                               by_auto(ax["nonneg_add"], psd_unfold)),
                      symbols=S)
    assert psd_add.verify()
    psd_congr = Theorem("psd_congr", "∀M P. psd(M) → psd(P.T*M*P)",
                        by_auto(ax["psd_def"], psd_unfold.instance(x=P*x)),
                        symbols=S)
    le_trans = Theorem(
        "le_trans", "∀A B C. mat_le(A, B) ∧ mat_le(B, C) → mat_le(A, C)",
        by_auto(ax["le_def"], psd_add.instance(M=C - B, N=B - A)),
        symbols=S)
    assert le_trans.verify()
    deps = le_trans.dependencies()
    assert {t.name for t in deps if t.kind == "axiom"} \
        == {"le_def", "psd_def", "nonneg_add"}
    assert [t.name for t in deps if t.kind == "ledger"] == ["quad_add"]

    pf = Proof("lyap_mono",
               "∀A B P R. mat_le(A, B) → mat_le(P.T*A*P + R, P.T*B*P + R)",
               symbols=S)
    assert pf.assumptions == ["h1"]
    pf.have("d", "psd(B - A)", by_auto(ax["le_def"]), using=["h1"])
    pf.have("e", "psd(P.T*(B - A)*P)", by_auto(psd_congr), using=["d"])
    with pytest.raises(ValueError):
        pf.have("e", "psd(A)", by_auto())                  # used name
    with pytest.raises(ValueError):
        pf.have("g", "psd(P.T*A*P)", by_auto(ax["le_def"], max_nodes=50),
                using=["d"])                               # wrong step
    pf.have("f", "psd(P.T*B*P - P.T*A*P)",
            by_steps(by_rewrite(congr_diff, reverse=True), by_auto()),
            using=["e"])
    th = pf.show(by_auto(ax["le_def"]), using=["f"])
    assert th.verified and pf.theorem is th
    assert th.statement == and_implication(
        "∀A B P R. mat_le(A, B) → mat_le(P.T*A*P + R, P.T*B*P + R)")
    s = str(pf)
    assert "assume h1: mat_le(A, B)" in s and s.endswith("qed")
    assert "ledger facts: congr_diff" in th.explain()

    # A non-theorem is rejected.
    with pytest.raises(ValueError):
        Theorem("le_sym", "∀A B. mat_le(A, B) → mat_le(B, A)",
                by_auto(ax["le_def"], psd_add, psd_congr, max_nodes=200),
                symbols=S).verify()


def and_implication (s):
    return Theorem("tmp", s, symbols=S).statement


def test_by_cas ():
    x1 = Symbol("x")
    Ps = Symbol("P")
    th = Theorem("p", PredApply(Ps, 2 * x1 + 1), axiom=True)
    t = Theorem("q", PredApply(Ps, (x1 + 1)**2 - x1**2),
                by_steps(by_cas("expand"), by_auto(th)))
    assert t.verify()
    assert t.explain().splitlines()[1] == "  proof: cas(expand); auto(p)"


def test_induction ():
    from sympy_matrix_tools.fol import by_induct, by_rewrite_asm
    k = Symbol("k", integer=True, nonnegative=True)
    s = Function("s")
    T = {"k": k, "s": s}
    nat_induct = Theorem(
        "nat_induct", "∀Phi k. Phi(0) ∧ (∀k. Phi(k) → Phi(k + 1)) → Phi(k)",
        axiom=True, symbols=T)
    s_0 = Theorem("s_0", "s(0) = 0", axiom=True, symbols=T)
    s_step = Theorem("s_step", "∀k. s(k + 1) = s(k) + 2*k + 1",
                     axiom=True, symbols=T)
    # The instance does not bind the predicate s or Eq.
    i = nat_induct.instance(Phi=Lambda(k, Eq(s(k), k**2)))
    assert [str(v) for v in i.variables()] == ["k"]
    tac = by_steps(by_induct(nat_induct, "k"), by_rewrite(s_0),
                   by_rewrite(s_step), by_rewrite_asm(), by_cas("expand"))
    sq = Theorem("sq", "∀k. s(k) = k**2", tac, symbols=T)
    assert sq.verify()
    assert {t.name for t in sq.dependencies()} \
        == {"nat_induct", "s_0", "s_step"}
    with pytest.raises(ValueError):
        Theorem("cube", "∀k. s(k) = k**3", tac, symbols=T).verify()
    with pytest.raises(ValueError):
        by_induct(s_step, "k")          # no predicate variable


def test_rewrite_occurrence ():
    x1 = Symbol("x")
    f = Function("f")
    Ps = Symbol("P")
    double = ForAll(Lambda(x1, Eq(f(x1), 2 * x1)))
    z = make_proofstate(PredApply(Ps, f(x1), f(x1) + 1))
    r = rewrite_goal(z, double, occurrence=1)
    assert get_syms_prems_concl(r)[1][0] == PredApply(Ps, f(x1), 2*x1 + 1)
    r = rewrite_goal(z, double)
    assert get_syms_prems_concl(r)[1][0] == PredApply(Ps, 2*x1, 2*x1 + 1)
    assert rewrite_goal(z, double, occurrence=2) is None
