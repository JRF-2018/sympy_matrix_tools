# -*- coding: utf-8 -*-
__version__ = '0.0.8' # Time-stamp: <2022-04-29T11:54:11Z>
## Language: Japanese/UTF-8

import sympy
from sympy import Mul, Add, MatMul, MatAdd, Integer, Identity, MatPow, Pow


def _MatMul_fixed_args_cnc (self, cset=False, warn=True, split_1=True):
  c = [x for x in self.args if x.is_commutative]
  nc = [x for x in self.args if not x.is_commutative]
  if cset:
    clen = len(c)
    c = set(c)
    if clen and warn and len(c) != clen:
      raise ValueError('repeated commutative arguments: %s' %
                       [ci for ci in c if list(self.args).count(ci) > 1])
  return [c, nc]


def fix_MatMul_args_cnc ():
  MatMul.args_cnc = _MatMul_fixed_args_cnc

  
def partial_apply (X, arg, applyf):
  repl = applyf(arg)
  return X.subs(arg, repl)

def _expand_pow (X, expand_pow=True, right=False):
  assert X.func == MatPow or X.func == Pow
  l = []
  pow = X.args[1]
  q = X.args[0]
  if pow.is_integer and pow.is_number:
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
    elif isinstance(expand_pow, sympy.Basic) and not expand_pow.is_number:
      if expand_pow == pow:
        return [X]
      if right:
        return [q ** expand_pow, q ** (pow - expand_pow)]
      else:
        return [q ** (pow - expand_pow),  q ** expand_pow]
    elif expand_pow is False or expand_pow == 0:
      return [X]
    elif expand_pow > 0:
      if right:
        return [q] * expand_pow + [q ** (pow - expand_pow)]
      else:
        return [q ** (pow - expand_pow)] + [q] * expand_pow
    else:
      if right: 
        return [q ** -1] * (-expand_pow) + [q ** (pow + expand_pow)]
      else:
        return [q ** (pow + expand_pow)] + [q ** -1] * (-expand_pow)


def mat_separate_cnc (X, expand_pow=False, cset=False, right=False):
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
          ValueError("cset must be False.")
      else:
        l.append(x)
    Xnc = l
    l = []
    for x in Xc:
      if x.func == Pow and x.args[1].is_integer:
        l += _expand_pow(x, expand_pow=expand_pow, right=right)
        if cset:
          ValueError("cset must be False.")
      else:
        l.append(x)
    Xc = l
  if cset:
    Xc = set(Xc)
  return Xc, Xnc
    

def mat_coeff (X, Y, right=False, expand_pow=False, nth=0):
  Yc, Ync = mat_separate_cnc(Y, expand_pow=expand_pow, cset=False, right=right)
  if X.func == MatAdd:
    args = X.args
  else:
    args = [X]
  for Z in args:
    Zc, Znc = mat_separate_cnc(Z, expand_pow=expand_pow, cset=False, right=right)
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
  Yc, Ync = mat_separate_cnc(Y, expand_pow=expand_pow, cset=False, right=right)
  if X.func != MatAdd:
    return X
  lco = []
  lnco = []
  for Z in X.args:
    Zc, Znc = mat_separate_cnc(Z, expand_pow=expand_pow, cset=False, right=right)
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
        ValueError("When adding, scalar is added by a non squared matrix.")
      A = A + Add(*list(lc)) * Identity(A.rows)
    if expand_inner:
      A = A.expand().doit()
    return Add(Mul(A, Y), *lnco)
  elif lc:
    A = Add(*list(lc))
    if expand_inner:
      A = A.expand().doit()
    return Add(Mul(A, Y), *lnco)
  else:
    return X


def mat_divide (X, Y, right=False):
  assert Y.is_commutative or Y.is_square
  if not Y.is_commutative:
    if right:
      assert X.rows == Y.cols
    if not right:
      assert X.cols == Y.rows
  Yc, Ync = mat_separate_cnc(Y, cset=False, expand_pow=True, right=not right)
  if X.func == MatAdd:
    args = X.args
  else:
    args = [X]
  l = []
  for Z in args:
    Zc, Znc = mat_separate_cnc(Z, cset=False, expand_pow=True, right=not right)
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
