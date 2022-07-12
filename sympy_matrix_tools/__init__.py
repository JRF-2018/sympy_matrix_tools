# -*- coding: utf-8 -*-
__version__ = '0.3.2' # Time-stamp: <2022-07-12T22:45:05Z>

import sympy
from packaging.version import parse as parse_version

from .matrix_tools import *
from .matrix_function import *
from .seq_tools import *
from .matrix_concrete import *
from .unif_tools import *
if parse_version(sympy.__version__) >= parse_version('1.10.1'):
    from .matrix_symbol import *
    from .apply import *
    from .logic_tools import *
