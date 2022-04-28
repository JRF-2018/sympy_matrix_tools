# -*- coding: utf-8 -*-
__version__ = '0.0.1' # Time-stamp: <2022-04-28T15:54:36Z>
## Language: Japanese/UTF-8

from sympy import Mul, Add, MatMul, MatAdd, Integer, Identity, MatPow, Pow


def partial_apply(X, arg, applyf):
  repl = applyf(arg)
  return X.subs(arg, repl)


def mat_separate_cnc(X, expand_pow=False, cset=True):
  assert X.func == MatMul
  Xc = [x for x in X.args if x.is_commutative]
  Xnc = [x for x in X.args if not x.is_commutative]
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
      if x.func == MatPow and x.args[1].is_integer and x.args[1].is_positive:
        l += [x.args[0]] * x.args[1]
        if cset:
          ValueError("cset must be False.")
      else:
        l.append(x)
    Xnc = l
    l = []
    for x in Xc:
      if x.func == Pow and x.args[1].is_integer and x.args[1].is_positive:
        l += [x.args[0]] * x.args[1]
        if cset:
          ValueError("cset must be False.")
      else:
        l.append(x)
    Xc = l
  if cset:
    Xc = set(Xc)
  return Xc, Xnc
    

def mat_coeff(X, Y, right=False):
  if Y.func == MatMul:
    Yc, Ync = mat_separate_cnc(Y)
  else:
    Yc = set()
    Ync = [Y]
  if X.func == MatAdd:
    args = X.args
  else:
    args = [X]
  for Z in args:
    if Z.func == MatMul:
      Zc, Znc = mat_separate_cnc(Z)
    else:
      Zc = set([])
      Znc = [Z]
    R = None
    if Yc <= Zc and len(Znc) >= len(Ync):
      if not right and Ync == Znc[-len(Ync):]:
        R = list(Zc - Yc) + Znc[0:-len(Ync)]
      if right and Ync == Znc[0:len(Ync)]:
        R = list(Zc - Yc) + Znc[len(Ync):]
    if R is not None:
      if len(R) == 0:
        return 1
      elif len(R) == 1:
        return R[0]
      else:
        return MatMul(*R)
  return None


def mat_collect(X, Y, right=False, expand_pow=False):
  if Y.func == MatMul:
    Yc, Ync = mat_separate_cnc(Y, expand_pow=expand_pow, cset=False)
  else:
    Yc = []
    Ync = [Y]
  if X.func != MatAdd:
    return X
  lco = []
  lnco = []
  for Z in X.args:
    if Z.func == MatMul:
      Zc, Znc = mat_separate_cnc(Z, expand_pow=expand_pow, cset=False)
    else:
      Zc = []
      Znc = [Z]
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
    R = None
    if len(Znc) >= len(Ync):
      if not right and Ync == Znc[-len(Ync):]:
        R = Zc + Znc[0:-len(Ync)]
      if right and Ync == Znc[0:len(Ync)]:
        R = Zc + Znc[len(Ync):]
    if R is not None:
      if len(R) == 0:
        lco.append(Integer(1))
      elif len(R) == 1:
        lco.append(R[0])
      else:
        lco.append(MatMul(*R))
    else:
        lnco.append(Z)
  lc = [x for x in lco if x.is_commutative]
  lnc = [x for x in lco if not x.is_commutative]
  if lnc:
    A = MatAdd(*lnc)
    if lc:
      if not A.is_square:
        ValueError("When adding, scalar is added by a non squared matrix.")
      A = MatAdd(A, MatMul(Add(*list(lc)), Identity(A.rows)))
    return MatAdd(MatMul(A, Y), *lnco)
  elif lc:
    return MatAdd(MatMul(Add(*list(lc)), Y), *lnco)
  else:
    return X


def mat_divide(X, Y, right=False):
  assert Y.is_square
  if right:
    assert X.rows == Y.cols
  if not right:
    assert X.cols == Y.rows
  if Y.func == MatMul:
    Yc, Ync = mat_separate_cnc(Y)
  else:
    Yc = set()
    Ync = [Y]
  if X.func == MatAdd:
    args = X.args
  else:
    args = [X]
  l = []
  for Z in args:
    if Z.func == MatMul:
      Zc, Znc = mat_separate_cnc(Z)
    else:
      Zc = set([])
      Znc = [Z]
    R = None
    if len(Znc) >= len(Ync):
      Q = Mul(*Zc) / Mul(*Yc)
      if right and Ync == Znc[-len(Ync):]:
        R = ([Q] if Q != 1 else []) + Znc[0:-len(Ync)]
      if not right and Ync == Znc[0:len(Ync)]:
        R = ([Q] if Q != 1 else []) + Znc[len(Ync):]
    if R is not None:
      if len(R) == 0:
        l.append(Identity(Y.rows))
      elif len(R) == 1:
        l.append(R[0])
      else:
        l.append(MatMul(*R))
    else:
      if right:
        l.append(MatMul(Z, Y ** -1))
      else:
        l.append(MatMul(Y ** -1, Z))
  return MatAdd(*l)
