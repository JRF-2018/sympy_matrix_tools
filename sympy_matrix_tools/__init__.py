# -*- coding: utf-8 -*-
__version__ = '0.2.0' # Time-stamp: <2022-06-19T07:59:00Z>

import sympy
from packaging.version import parse as pv

from .matrix_tools import *
from .matrix_function import *
from .seq_tools import *
from .matrix_concrete import *
from .apply import *
from .unif_tools import *
if pv(sympy.__version__) >= pv('1.10.1'):
    from .matrix_symbol import *
