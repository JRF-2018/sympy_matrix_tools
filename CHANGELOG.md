# Changelog

All notable changes are recorded here.  The format is based on
[Keep a Changelog](https://keepachangelog.com/).  Each version is
tagged as `vX.Y.Z` on GitHub.


## [0.10.0] - 2026-11-03

### Added
- `by_induct(rule, var)`: induction on a bound variable by an
  induction axiom with a predicate variable, like `induct` of
  Isabelle.
- `by_rewrite_asm` (`logic_tools.rewrite_goal_by_assumption`):
  rewriting by an equation among the premises of a subgoal, e.g. an
  induction hypothesis.
- `by_rewrite(..., occurrence=i)` and `rewrite_goal(..., occurrence=i)`
  rewrite only one occurrence.
- `examples/loewner_order.ipynb`: Vk ≼ Vk+1 and 0 ≼ Vk for all k, with
  induction introduced as an axiom.
- `docs/ai_guide.md`: "Axioms are fine", on what the framework is for
  (plausible arguments with every step checked, not proofs from first
  principles), and Recipe 6 (induction).

### Fixed
- `Theorem.instance` bound the predicates and the dimensions in the
  new terms (e.g. `mat_le` and `n` in a `Lambda`) as variables.


## [0.9.0] - 2026-11-03

### Added
- Logic with matrices: `parse_formula` and `Theorem` take `symbols`
  whose `MatrixSymbol`s may be bound (`∀M. psd(M) → ...`), and read
  equations `s = t`; `formula_str` prints `Eq` as `s = t`.
- Tactics in `sympy_matrix_tools.fol`: `by_rule` (one backward step),
  `by_rewrite` (rewriting by an equation theorem, as the substitution
  of equals), `by_cas` (a computation of SymPy on the terms of a
  subgoal) and `by_steps` (a sequence).  `by_auto` uses a theorem whose
  conclusion is `↔` or `∧` as its parts.  In `logic_tools`:
  `rewrite_goal`, `iter_rewrite_goal` and `cas_goal`.
- `Theorem.instance(x=t)` (forall elimination with explicit terms) and
  `Theorem.cas(func)` (terms computed by SymPy).
- `Proof`: structured proofs in forward style like Isar (`have` ...
  `show`); each step is proved at once.
- `Derivation.theorem(name)`: a derived fact of a ledger as a
  `Theorem` of kind "ledger", whose `verify()` requires symbolic and
  numerical checks of all its steps.
- `Theorem.dependencies()` and `Theorem.explain()`: the axioms, ledger
  facts and theorems that a theorem rests on.
- `ho_unify`: associative-commutative unification (`?M + ?N` with
  `A + B + C`) and associative unification of `MatMul` (`P.T*?M*P`
  with `P.T*A*B*P`).
- `examples/loewner_order.ipynb`: the Loewner order defined in the
  notebook, with two ledgers and backward and forward proofs.
- `docs/ai_guide.md`: a recipe for statements about matrices and a
  workflow for AI agents asked to prove or check an argument; the
  README has a prompt template.
- The weekly workflow has a manual input `force_failure`, which makes
  the tests fail on purpose to try the path that opens an issue.

### Changed
- Importing `sympy_matrix_tools.fol` applies
  `fix_Basic_MatrixSymbol_about_dummy()`.
- The normal form of formulas evaluates `MatMul`, `MatAdd` and
  `Transpose` as their operators do (e.g. `(P*x).T` is `x.T*P.T`).
- `Theorem.verify()` verifies the theorems used first, and runs a
  proof only once.

### Fixed
- Metavariables of matrices were lifted to a `MatApply` without its
  shape, and patterns with `MatrixSymbol` arguments were not solved,
  so the proofstates could not be used with matrices.
- The fallback matchers of SymPy could raise errors on matrix
  expressions during unification.


## [0.8.0] - 2026-10-23

### Added
- `parse_formula` (`sympy_matrix_tools.logic_parser`): reads formulas
  in the logical notation that `formula_str` prints (with ASCII
  alternatives), including metavariables and proofstates.
- `Theorem` accepts a string as its statement.
- `.github/workflows/sympy-latest.yml`: weekly tests with the latest
  SymPy, which open an issue on failure (only one at a time) and close
  it when the tests pass again.


## [0.7.1] - 2026-10-23

### Changed
- The documents in `docs/` explain in words what each function and
  class is for (MatrixWild, MatrixDummy, Apply, the steps of the
  proofstates, the parameters of `Sum_step_forward`, ...), with tables
  of the main names.

### Fixed
- `make_wild` raised `TypeError` for a `MatrixSymbol`.


## [0.7.0] - 2026-10-23

### Added
- Documents split into `docs/*.md` by topic, with a short README.
- `docs/ai_guide.md` and `llms.txt` for AI agents.
- This changelog.
- Example notebooks in `examples/` (`lq_control.ipynb`,
  `chain_rule.ipynb`, `first_order_logic.ipynb`), executed by
  `tests/test_examples.py`.

- `Derivation` checks expressions with undefined functions
  (`MatrixFunction` and `Function`) numerically, by substituting random
  polynomials for them.
- `derivative_str`: shows derivatives of undefined functions as
  `x'(t)`, `F'(a)`, `D_2F(a, b)` and `DF(a)[d]`; used in the reports
  of `Derivation`.
- Module docstrings and docstrings of the main functions.

### Changed
- `rule(name)` also rewrites a part of a sum and a contiguous part of
  a product, and the symbolic check identifies a (1, 1) matrix with
  its transpose using the symmetry facts cited in the step.
- The example of the LQ control derives also the optimal value.
- `pytest.ini` runs all the tests including doctests of the documents.

### Removed
- `examples/lq_control.py`, replaced by `examples/lq_control.ipynb`.
- pytest runs the doctests in `README.md` and `docs/*.md`.


## [0.6.0] - 2026-10-23

### Added
- `sympy_matrix_tools.derivation`: `Derivation`, a ledger of
  derivations.  Steps are justified by `cas(...)`, `rule(name)` or
  nothing (asserted), and checked symbolically and numerically with
  random rational matrices.  `solve_linear` with obligations, text and
  LaTeX reports, display aliases.
- `Grad`, an unevaluated gradient by a matrix.
- `examples/lq_control.py`.


## [0.5.0] - 2026-10-23

### Added
- `sympy_matrix_tools.hounif`: `ho_unify`, higher-order pattern
  unification with alpha-equivalence of `Lambda`.
- `auto` and `iter_auto`: proof search by backtracking with iterative
  deepening, with safe eliminations.
- `split_goals`, `formula_str` and `print_proofstate(..., pretty=True)`.
- `sympy_matrix_tools.fol`: first-order logic on the minimal core.
  `Exists`, `Or` and `Bot` are defined impredicatively; the rules and
  basic theorems are proved by `auto`.

### Fixed
- An unsound behavior of `remove_trivial_assumptions`: bound variables
  of a premise were dropped and could be unified with outer
  metavariables.
- The unifiers of SymPy's `match` and `unify` are validated before
  use.
- `check_most_trivial` was called with a wrong argument.


## [0.4.0] - 2026-10-22

### Added
- The chain rule for `MatrixFunction` by a scalar and by a matrix
  (`MatrixDerivative`, `MatrixSubs`, `ArraySubs`).
- `mat_diff`, which avoids bugs of `convert_array_to_matrix`.
- `mat_expand`, since `MatMul.expand()` of SymPy 1.14 does not expand
  scalar coefficients.
- `fix_Transpose_eval_derivative`, applied on import.

### Changed
- Requires SymPy >= 1.14; branches for older SymPy are removed.
- `fix_MatMul_args_cnc` does nothing (fixed in SymPy 1.11).

### Fixed
- Import errors and assumptions of `MatrixDummy` and `MatrixWild`
  with SymPy 1.14.
- A typo in `make_wild`.


## [0.3.6] and earlier (2022)

The history is logged in Japanese at
http://jrf.cocolog-nifty.com/statuses/2022/06/post-2c0e50.html .

- 0.3.x (2022-07 to 2022-10): Isabelle-style proofstates
  (`make_proofstate`, `melt_theorem`, `iter_resolve_implications`,
  `eresolve_implications`, `dresolve_implications`,
  `remove_trivial_assumptions`, ...); support for SymPy 1.11 (0.3.6).
- 0.2.x (2022-06 to 2022-07): `MatrixWild`, `MatrixDummy`, `Apply`,
  `MatApply`, `PredApply`, unification and rewriting
  (`rewriterules_from_eqs`, `resolve_mp`), renaming bound symbols,
  `conditional_apply`, and the first predicate logic with `ForAll`
  (0.2.7).
- 0.1.x (2022-05): `MatrixFunction`, functions for sequences
  (`Sum_step_forward`, ...), `MatSum`, `MatProduct`, the fix of
  `ExpectationMatrix.expand`.
- 0.0.x (2022-04): the first versions: `mat_coeff`, `mat_collect`,
  `mat_divide` and the fix of `MatMul.args_cnc`.

[0.8.0]: https://github.com/JRF-2018/sympy_matrix_tools/compare/v0.7.1...v0.8.0
[0.7.1]: https://github.com/JRF-2018/sympy_matrix_tools/compare/v0.7.0...v0.7.1
[0.7.0]: https://github.com/JRF-2018/sympy_matrix_tools/compare/v0.6.0...v0.7.0
[0.6.0]: https://github.com/JRF-2018/sympy_matrix_tools/compare/v0.5.0...v0.6.0
[0.5.0]: https://github.com/JRF-2018/sympy_matrix_tools/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/JRF-2018/sympy_matrix_tools/compare/v0.3.6...v0.4.0
[0.3.6]: https://github.com/JRF-2018/sympy_matrix_tools/releases/tag/v0.3.6
