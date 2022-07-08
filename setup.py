# -*- coding: utf-8 -*-
__version__ = '0.3.0' # Time-stamp: <2022-07-08T03:18:16Z>
## Language: Japanese/UTF-8

from setuptools import setup, find_packages

setup(name='sympy_matrix_tools',
      version=__version__,
      description='Some tools for sympy matrices.',
      author='JRF',
      url='https://github.com/JRF-2018/sympy_matrix_tools',
      packages=find_packages(exclude=['tests']),)
