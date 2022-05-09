# -*- coding: utf-8 -*-
__version__ = '0.1.4' # Time-stamp: <2022-05-09T10:59:22Z>

import sympy
from .matrix_tools import *
from sympy import Mul, Add, MatMul, preorder_traversal, Sum, Product,\
  Dummy, sympify


def atoms_list (z, *types):
  if types:
    types = tuple(
      [t if isinstance(t, type) else type(t) for t in types])
  nodes = preorder_traversal(z)
  if types:
    result = [node for node in nodes if isinstance(node, types)]
  else:
    result = [node for node in nodes if not node.args]
  return result


def nth_Sum (z, idx):
  return atoms_list(z, Sum)[idx]


def nth_Product (z, idx):
  return atoms_list(z, Product)[idx]


def Sum_Product_expand (z, expand_inner=True):
  assert issubclass(z.func, Sum) or issubclass(z.func, Product)
  if issubclass(z.func, Sum):
    op1 = Add
  else:
    op1 = Mul
  if not issubclass(z.args[0].func, op1):
    return z
  if issubclass(z.args[0].func, MatMul):
    c, nc = mat_separate_cnc(z.args[0])
  else:
    c = z.args[0].args
    nc = []
  l = [z.func(i, z.args[1:]) for i in c]
  if nc:
    l += [z.func(z.args[0].func(*nc), z.args[1:])]
  if expand_inner:
    l = [i.expand().doit() for i in l]
  return op1(*l)


def Sum_expand (z, expand_inner=True):
  for x in atoms_list(z, Sum):
    z = z.subs(x, Sum_Product_expand(x, expand_inner=expand_inner))
  return z


def Product_expand (z, expand_inner=True):
  for x in atoms_list(z, Product):
    z = z.subs(x, Sum_Product_expand(x, expand_inner=expand_inner))
  return z


def Sum_collect (z, level=1):
  if not issubclass(z.func, Add) or level < 1:
    return z
  d = {}
  lo = []
  for x in z.args:
    if issubclass(x.func, Sum) and len(x.args) >= level + 1:
      xargs = x.args
      v, beg, ed = xargs[-level]
      t = (beg, ed) + xargs[len(xargs) - level + 1:]
      if t not in d:
        d[t] = []
      d[t] += [(xargs[0:len(xargs) - level], v)]
    else:
      lo.append(x)
  for t in d:
    l = d[t]
    beg, ed, rargs = t[0], t[1], t[2:]
    if len(l) == 1:
      x, largs, v = l[0][0][0], l[0][0][1:], l[0][1]
      args = list(largs) + [(v, beg, ed)] + list(rargs)
      lo.append(Sum(x, *args))
    else:
      v = l[0][1]
      if all([v1 == v for xl, v1 in l]):
        xi = v
      else:
        xi = Dummy('xi')
      l2 = []
      for xl, v in l :
        x, largs = xl[0], xl[1:]
        if largs:
          x = Sum(x, largs)
        l2.append(x.subs(v, xi))
      lo.append(Sum(Add(*l2), (xi, beg, ed), *rargs))
  return Add(*lo)


def Sum_swap (z, nth=None, first=1):
  y = None
  assert first >= 1
  if nth is not None:
    y = nth_Sum(z, nth)
    if y.func != Sum or len(y.args) < first + 2:
      raise ValueError("No swappables.")
  else:
    for x in atoms_list(z, Sum):
      if x.func == Sum and len(x.args) >= first + 2:
        y = x
        break
  if y is None:
    raise ValueError("No swappables.")
  if first == 1:
    x = y.args[0]
  else:
    x = Sum(*y.args[0:first])
  fv, fbg, fed = y.args[first]
  sv, sbg, sed = y.args[first + 1]
  rargs = y.args[first + 2:]
  if sv in fbg.free_symbols \
     or sv in fed.free_symbols:
    raise ValueError("The variable is not free.")
  if fv in fbg.free_symbols \
     or fv in fed.free_symbols \
     or fv in sbg.free_symbols \
     or fv in sed.free_symbols:
    xi = Dummy("xi")
    y2 = Sum(x.subs(fv, xi), (sv, sbg, sed), (xi, fbg, fed ), *rargs)
  else:
    y2 = Sum(x, (sv, sbg, sed), (fv, fbg, fed), *rargs)
  return z.subs(y, y2)


def Sum_Product_step (z, step=1, begin=None, end=None,
                      minus=False, where=-1, forward=True):
  assert issubclass(z.func, Sum) or issubclass(z.func, Product)
  if where < 0:
    where = len(z.args) + where
  assert 1 <= where < len(z.args)
  step = sympify(step)
  if begin is not None:
    begin = sympify(begin)
  if end is not None:
    end = sympify(end)

  op0 = z.func
  if issubclass(z.func, Sum):
    op1 = Add
    inv = lambda x: -x
  else:
    op1 = Mul
    inv = lambda x: x ** (-1)
  if where > 1:
    y = op0(*z.args[0:where])
  else:
    y = z.args[0]
  v, beg, ed = z.args[where]
  rargs = z.args[where + 1:]
  l = []
  r = []
  if step.is_integer and step.is_number:
    if step >= 0:
      for i in range(step):
        if forward:
          l.append(y.subs(v, beg))
          beg += 1
        else:
          r.append(y.subs(v, ed))
          ed -= 1
    else:
      for i in range(-step):
        if forward:
          l.append(inv(y.subs(v, beg - 1)))
          beg -= 1
        else:
          r.append(inv(y.subs(v, ed + 1)))
          ed += 1
    if beg == ed:
      q = op1(*(l + [y.subs(v, beg)] +  r))
      if rargs:
        return op0(q, *rargs)
      else:
        return q
    if begin is not None and beg != begin:
      y = y.subs(v, v - begin  + beg)
      ed = ed + begin - beg
      beg = begin
    if end is not None and ed != end:
      y = y.subs(v, v - end + ed)
      beg = beg + end - ed
      ed = end
    q = op1(*(l + [op0(y, (v, beg, ed))] + r))
    if rargs:
      return op0(q, *rargs)
    else:
      return q
  else: # step may be symbol.
    if (minus is not None and minus) \
       or (minus is None and step < 0):
      if forward:
        yc = y
        yl = inv(y)
        begl, edl = beg + step, beg - 1
        begc, edc = beg + step, ed
        if begin is not None:
          if begl != begin:
            yl = yl.subs(v, v - begin  + begl)
            edl = edl + begin - begl
            begl = begin
          if begc != begin:
            yc = yc.subs(v, v - begin  + begc)
            edc = edc + begin - begc
            begc = begin
        if end is not None:
          if edl != end:
            yl = yl.subs(v, v - end + edl)
            begl = begl + end - edl
            edl = end
          if edc != end:
            yc = yc.subs(v, v - end + edc)
            begc = begc + end - edc
            edc = end
        q = op1(*[op0(yl, (v, begl, edl)), op0(yc, (v, begc, edc))])
        if rargs:
          return op0(q, *rargs)
        else:
          return q
      else:
        yc = y
        yr = inv(y)
        begc, edc = beg, ed - step
        begr, edr = ed + 1, ed - step
        if begin is not None:
          if begc != begin:
            yc = yc.subs(v, v - begin  + begc)
            edc = edc + begin - begc
            begc = begin
          if begr != begin:
            yr = yr.subs(v, v - begin  + begr)
            edr = edr + begin - begr
            begr = begin
        if end is not None:
          if edc != end:
            yc = yc.subs(v, v - end + edc)
            begc = begc + end - edc
            edc = end
          if edr != end:
            yr = yr.subs(v, v - end + edr)
            begr = begr + end - edr
            edr = end
        q = op1(*[op0(yc, (v, begc, edc)), op0(yr, (v, begr, edr))])
        if rargs:
          return op0(q, *rargs)
        else:
          return q
    else:
      if forward:
        yl, yc = y, y
        begl, edl = beg, beg + step
        begc, edc = beg + step + 1, ed
        if begin is not None:
          if begl != begin:
            yl = yl.subs(v, v - begin  + begl)
            edl = edl + begin - begl
            begl = begin
          if begc != begin:
            yc = yc.subs(v, v - begin  + begc)
            edc = edc + begin - begc
            begc = begin
        if end is not None:
          if edl != end:
            yl = yl.subs(v, v - end + edl)
            begl = begl + end - edl
            edl = end
          if edc != end:
            yc = yc.subs(v, v - end + edc)
            begc = begc + end - edc
            edc = end
        q = op1(*[op0(yl, (v, begl, edl)), op0(yc, (v, begc, edc))])
        if rargs:
          return op0(q, *rargs)
        else:
          return q
      else:
        yr, yc = y, y
        begc, edc = beg, ed - step - 1
        begr, edr = ed - step, ed
        if begin is not None:
          if begc != begin:
            yc = yc.subs(v, v - begin  + begc)
            edc = edc + begin - begc
            begc = begin
          if begr != begin:
            yr = yr.subs(v, v - begin  + begr)
            edr = edr + begin - begr
            begr = begin
        if end is not None:
          if edc != end:
            yc = yc.subs(v, v - end + edc)
            begc = begc + end - edc
            edc = end
          if edr != end:
            yr = yr.subs(v, v - end + edr)
            begr = begr + end - edr
            edr = end
        q = op1(*[op0(yc, (v, begc, edc)), op0(yr, (v, begr, edr))])
        if rargs:
          return op0(q, *rargs)
        else:
          return q


def Sum_step_forward (z, step=1, begin=None, end=None,
                      minus=False, where=-1, nth=0):
  y = nth_Sum(z, nth)
  if y is None:
    raise ValueError("Illegal index.")
  y2 = Sum_Product_step(y, step=step, begin=begin, end=end,
                        minus=minus, where=where, forward=True)
  return z.subs(y, y2)


def Sum_step_backward (z, step=1, begin=None, end=None,
                       minus=False, where=-1, nth=0):
  y = nth_Sum(z, nth)
  if y is None:
    raise ValueError("Illegal index.")
  y2 = Sum_Product_step(y, step=step, begin=begin, end=end,
                        minus=minus, where=where, forward=False)
  return z.subs(y, y2)


def Product_step_forward (z, step=1, begin=None, end=None,
                          minus=False, where=-1, nth=0):
  y = nth_Product(z, nth)
  if y is None:
    raise ValueError("Illegal index.")
  y2 = Sum_Product_step(y, step=step, begin=begin, end=end,
                        minus=minus, where=where, forward=True)
  return z.subs(y, y2)


def Product_step_backward (z, step=1, begin=None, end=None,
                           minus=False, where=-1, nth=0):
  y = nth_Product(z, nth)
  if y is None:
    raise ValueError("Illegal index.")
  y2 = Sum_Product_step(y, step=step, begin=begin, end=end,
                        minus=minus, where=where, forward=False)
  return z.subs(y, y2)
