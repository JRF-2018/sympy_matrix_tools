# -*- coding: utf-8 -*-
"""Functions for non-commutative matrix expressions.

mat_coeff, mat_collect, mat_divide and mat_trivial_divide find
coefficients, collect terms and divide in matrix expressions, where
the order of factors matters.  mat_expand expands scalar coefficients
too, and partial_apply applies a function to a subexpression.  It also
has fix_ExpectationMatrix_expand for sympy.stats.

See docs/matrix_tools.md.
"""
__version__ = '0.7.0' # Time-stamp: <2026-10-22T17:31:49Z>

from sympy import Basic, Mul, Add, MatMul, MatAdd, Integer, \
  Identity, MatPow, Pow, Inverse, Transpose, powsimp


def fix_MatMul_args_cnc ():
  # MatMul.args_cnc was fixed in SymPy 1.11.  Left for compatibility.
  pass


def _ExpectationMatrix_fixed_expand (self, **hints):
  from sympy.stats import Expectation
  from sympy.stats.rv import is_random
  from sympy.core.function import expand as _expand
  expr = self.args[0]
  condition = self._condition
  if not is_random(expr):
      return expr

  if isinstance(expr, Add):
      return Add.fromiter(Expectation(a, condition=condition).expand()
              for a in expr.args)

  elif isinstance(expr, Transpose):
      return expr.func(Expectation(expr.args[0],
                                   condition=condition).expand())

  expand_expr = _expand(expr)
  if isinstance(expand_expr, Add):
      return Add.fromiter(Expectation(a, condition=condition).expand()
              for a in expand_expr.args)

  elif isinstance(expand_expr, Transpose):
      return expand_expr.func(Expectation(expand_expr.args[0],
                                          condition=condition).expand())

  elif isinstance(expr, (Mul, MatMul)):
      rv = []
      nonrv = []
      postnon = []

      for a in expr.args:
          if is_random(a):
              if rv:
                  rv.extend(postnon)
              else:
                  nonrv.extend(postnon)
              postnon = []
              rv.append(a)
          elif a.is_Matrix:
              postnon.append(a)
          else:
              nonrv.append(a)

      # In order to avoid infinite-looping (MatMul may call .doit() again),
      # do not rebuild
      if len(nonrv) == 0 and len(postnon) == 0:
          return self
      return Mul.fromiter(nonrv)*Expectation(Mul.fromiter(rv),
              condition=condition).expand()*Mul.fromiter(postnon)

  return self


def fix_ExpectationMatrix_expand ():
  from sympy.stats import ExpectationMatrix
  ExpectationMatrix.expand = _ExpectationMatrix_fixed_expand


def mat_expand (X, **hints):
  """expand() that also expands and powsimp's scalar coefficients of MatMul.

  MatMul.expand() of SymPy >= 1.12 leaves scalar coefficients
  unexpanded, and x1*x1**(N - 2) is no longer combined.

  >>> from sympy import Symbol, MatrixSymbol, Mul
  >>> x1 = Symbol("x1")
  >>> N = Symbol("N", integer=True)
  >>> M1 = MatrixSymbol("M1", N, N)
  >>> mat_expand(x1 * x1 ** (N - 2) * (M1 + M1 ** 2))
  x1**(N - 1)*M1**2 + x1**(N - 1)*M1
  """
  X = X.expand(**hints).doit()
  return X.replace(
    lambda e: isinstance(e, MatMul),
    lambda e: MatMul(*[powsimp(x.expand(**hints)) if x.is_commutative else x
                       for x in e.args]).doit())


def partial_apply (X, arg, applyf):
  """Replace the subexpression arg of X by applyf(arg).

  >>> from sympy import MatrixSymbol, Symbol
  >>> N = Symbol("N", integer=True)
  >>> M1 = MatrixSymbol("M1", N, N)
  >>> M2 = MatrixSymbol("M2", N, N)
  >>> z = (M1 + M2) * M1 + M2
  >>> partial_apply(z, (M1 + M2) * M1, lambda y: y.expand())
  (M1**2 + M2*M1) + M2
  """
  repl = applyf(arg)
  return X.subs(arg, repl)

def _expand_pow (X, expand_pow=True, right=True):
  assert X.func == MatPow or X.func == Pow
  l = []
  pow = X.args[1]
  q = X.args[0]
  if not (isinstance(expand_pow, dict) and X in expand_pow
          and expand_pow[X] is not True) \
     and pow.is_integer and pow.is_number:
    if pow == 0:
      if q.is_commutative:
        q = Integer(1)
      else:
        q = Identity(q.rows)
      pow = 1
    elif pow < 0:
      q = q ** -1
      pow = -pow
    return [q] * pow
  else:
    if isinstance(expand_pow, dict):
      if X not in expand_pow:
        if 'default' in expand_pow:
          expand_pow = expand_pow['default']
        if True in expand_pow:
          expand_pow = expand_pow[True]
        else:
          return [X]
      else:
        expand_pow = expand_pow[X]
    if expand_pow is True:
      return [X]
    elif isinstance(expand_pow, Basic) and not expand_pow.is_number:
      if expand_pow == pow:
        return [X]
      if right:
        return [q ** (pow - expand_pow),  q ** expand_pow]
      else:
        return [q ** expand_pow, q ** (pow - expand_pow)]
    elif expand_pow is False or expand_pow == 0:
      return [X]
    elif expand_pow > 0:
      if expand_pow == pow:
        return [q] * expand_pow
      if right:
        return [q ** (pow - expand_pow)] + [q] * expand_pow
      else:
        return [q] * expand_pow + [q ** (pow - expand_pow)]
    else:
      if expand_pow == pow:
        return [q ** -1] * (- expand_pow)
      if right: 
        return [q ** (pow + expand_pow)] + [q ** -1] * (-expand_pow)
      else:
        return [q ** -1] * (-expand_pow) + [q ** (pow + expand_pow)]


def mat_separate_cnc (X, expand_pow=False, cset=False, right=True):
  if X.func == MatMul:
    args = X.args
  else:
    args = [X]
  Xc = [x for x in args if x.is_commutative]
  Xnc = [x for x in args if not x.is_commutative]
  l = []
  for x in Xc:
    if x.func == Mul:
      l += x.args
    else:
      l.append(x)
  Xc = l
  if expand_pow:
    l = []
    for x in Xnc:
      if x.func == MatPow and x.args[1].is_integer:
        l += _expand_pow(x, expand_pow=expand_pow, right=right)
        if cset:
          raise ValueError("cset must be False.")
      else:
        l.append(x)
    Xnc = l
    l = []
    for x in Xc:
      if x.func == Pow and x.args[1].is_integer:
        l += _expand_pow(x, expand_pow=expand_pow, right=right)
        if cset:
          raise ValueError("cset must be False.")
      else:
        l.append(x)
    Xc = l
  if cset:
    Xc = set(Xc)
  return Xc, Xnc
    

def mat_coeff (X, Y, right=False, expand_pow=False, nth=0):
  """Return the coefficient C of Y in a term C*Y of X.

  If right is True, look for a term Y*C instead.  The coefficient of
  the nth matching term (from 0) is returned, not the sum; use
  mat_collect to collect all terms.  If expand_pow is True, powers
  such as M**2 are treated as products M*M.  Return None if no term
  matches.

  >>> from sympy import MatrixSymbol, Symbol
  >>> N = Symbol("N", integer=True)
  >>> M1 = MatrixSymbol("M1", N, N)
  >>> M2 = MatrixSymbol("M2", N, N)
  >>> z = M1 * M2 + 3 * M2 * M1 + M1 ** 2
  >>> mat_coeff(z, M1)
  3*M2
  >>> mat_coeff(z, M1, right=True)
  M2
  >>> mat_coeff(z, M1, nth=1, expand_pow=True)
  M1
  """
  Yc, Ync = mat_separate_cnc(Y, expand_pow=expand_pow, cset=False,
                             right=not right)
  if X.func == MatAdd:
    args = X.args
  else:
    args = [X]
  for Z in args:
    Zc, Znc = mat_separate_cnc(Z, expand_pow=expand_pow, cset=False,
                               right=not right)
    R = None
    if len(Znc) >= 1:
      if not right and Y == Znc[-1]:
        R = Zc + Znc[0:-1]
      if right and Y == Znc[0]:
        R = Zc + Znc[1:]
    if R is None:
      done = True
      for x in Yc:
        try:
          Zc.remove(x)
        except:
          done = False
          break
      if not done:
        continue
    if R is None and len(Znc) >= len(Ync):
      if not right and Ync == Znc[len(Znc)-len(Ync):]:
        R = Zc + Znc[0:len(Znc)-len(Ync)]
      if right and Ync == Znc[0:len(Ync)]:
        R = Zc + Znc[len(Ync):]
    if R is not None:
      if nth != 0:
        nth -= 1
        continue
      if len(R) == 0:
        return Integer(1)
      elif len(R) == 1:
        return R[0]
      else:
        return MatMul(*R)
  return None


def mat_collect (X, Y, right=False, expand_pow=False, expand_inner=False):
  """Collect the terms of the form C*Y in X into (sum of C)*Y.

  If right is True, collect Y*C into Y*(sum of C).  Scalar terms are
  multiplied by the identity matrix.  expand_pow treats powers as
  products (see also the dict form in docs/matrix_tools.md), and
  expand_inner expands the collected coefficient by mat_expand.

  >>> from sympy import MatrixSymbol, Symbol
  >>> N = Symbol("N", integer=True)
  >>> M1 = MatrixSymbol("M1", N, N)
  >>> M2 = MatrixSymbol("M2", N, N)
  >>> mat_collect(M1 * M2 + M2 * M2 + 3 * M2, M2, expand_pow=True)
  (3*I + M1 + M2)*M2
  """
  Yc, Ync = mat_separate_cnc(Y, expand_pow=expand_pow, cset=False,
                             right=not right)
  if X.func != MatAdd:
    return X
  lco = []
  lnco = []
  for Z in X.args:
    Zc, Znc = mat_separate_cnc(Z, expand_pow=expand_pow, cset=False,
                               right=not right)
    R = None
    if len(Znc) >= 1:
      if not right and Y == Znc[-1]:
        R = Zc + Znc[0:-1]
      if right and Y == Znc[0]:
        R = Zc + Znc[1:]
    if R is None:
      done = True
      for x in Yc:
        try:
          Zc.remove(x)
        except:
          done = False
          break
      if not done:
        lnco.append(Z)
        continue
    if R is None and len(Znc) >= len(Ync):
      if not right and Ync == Znc[len(Znc)-len(Ync):]:
        R = Zc + Znc[0:len(Znc)-len(Ync)]
      if right and Ync == Znc[0:len(Ync)]:
        R = Zc + Znc[len(Ync):]
    if R is not None:
      if len(R) == 0:
        lco.append(Integer(1))
      elif len(R) == 1:
        lco.append(R[0])
      else:
        lco.append(Mul(*R))
    else:
        lnco.append(Z)
  lc = [x for x in lco if x.is_commutative]
  lnc = [x for x in lco if not x.is_commutative]
  if lnc:
    if len(lnc) == 1:
      A = lnc[0]
    else:
      A = Add(*lnc)
    if lc:
      if not A.is_square:
        raise ValueError("When adding, scalar is added by a non squared matrix.")
      A = A + Add(*list(lc)) * Identity(A.rows)
    if expand_inner:
      A = mat_expand(A)
    if right:
      return Add(Mul(Y, A), *lnco)
    else:
      return Add(Mul(A, Y), *lnco)
  elif lc:
    A = Add(*list(lc))
    if expand_inner:
      A = mat_expand(A)
    return Add(Mul(A, Y), *lnco)
  else:
    return X


def mat_divide (X, Y, right=False):
  """Multiply X by Y**-1 from the left (from the right if right is
  True), cancelling Y in each term where it is possible.

  >>> from sympy import MatrixSymbol, Symbol
  >>> N = Symbol("N", integer=True)
  >>> M1 = MatrixSymbol("M1", N, N)
  >>> M2 = MatrixSymbol("M2", N, N)
  >>> mat_divide(M1 * M2, M1)
  M2
  >>> mat_divide(M1 * M2, M2)
  M2**(-1)*M1*M2
  """
  assert Y.is_commutative or Y.is_square
  if not Y.is_commutative:
    if right:
      assert X.cols == Y.rows
    if not right:
      assert X.rows == Y.cols
  Yc, Ync = mat_separate_cnc(Y, cset=False, expand_pow=True,
                             right=right)
  if X.func == MatAdd:
    args = X.args
  else:
    args = [X]
  l = []
  for Z in args:
    Zc, Znc = mat_separate_cnc(Z, cset=False, expand_pow=True,
                               right=right)
    R = None
    if len(Znc) >= 1:
      Q = Mul(*Zc)
      if right and Y == Znc[-1]:
        R = ([Q] if Q != 1 else []) + Znc[0:-1]
      if not right and Y == Znc[0]:
        R = ([Q] if Q != 1 else []) + Znc[1:]
    if R is None and len(Znc) >= len(Ync):
      Q = Mul(*Zc) / Mul(*Yc)
      if right and Ync == Znc[len(Znc)-len(Ync):]:
        R = ([Q] if Q != 1 else []) + Znc[0:len(Znc)-len(Ync)]
      if not right and Ync == Znc[0:len(Ync)]:
        R = ([Q] if Q != 1 else []) + Znc[len(Ync):]
    if R is not None:
      if len(R) == 0:
        R2 = Integer(1)
      elif len(R) == 1:
        R2 = R[0]
      else:
        R2 = Mul(*R)
      if R2.is_commutative:
        l.append(R2 * Identity(Y.rows))
      else:
        l.append(R2)
    else:
      if right:
        l.append(Mul(Z, Y ** -1))
      else:
        l.append(Mul(Y ** -1, Z))
  return Add(*l)


def mat_trivial_divide (X):
    """Cancel trivial products of inverses such as (A*B)**-1*A*B in X.

    >>> from sympy import MatrixSymbol, Symbol
    >>> N = Symbol("N", integer=True)
    >>> N2 = Symbol("N2", integer=True)
    >>> M1 = MatrixSymbol("M1", N, N)
    >>> M3 = MatrixSymbol("M3", N, N2)
    >>> M4 = MatrixSymbol("M4", N2, N)
    >>> mat_trivial_divide(M1 + ((M3 * M4) ** -2) * M3 * M4)
    (M3*M4)**(-1) + M1
    """
    d = {}
    for Z in X.atoms(MatMul):
        for Y in Z.atoms(Inverse):
            if issubclass(Y.args[0].func, MatMul):
                Q = Y.args[0].args
            else:
                Q = [Y.args[0]]
            d[MatMul(*Q, Y)] = Identity(Y.rows)
            d[MatMul(Y, *Q)] = Identity(Y.rows)
        for Y in Z.atoms(MatPow):
            if Y.args[1].is_number and Y.args[1] < 0:
                if issubclass(Y.args[0].func, MatMul):
                    Q = Y.args[0].args
                else:
                    Q = [Y.args[0]]
                d[MatMul(*Q, Y)] = Y.args[0] ** (Y.args[1] + 1)
                d[MatMul(Y, *Q)] = Y.args[0] ** (Y.args[1] + 1)
    X = X.subs(d)
    return X
