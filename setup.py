# -*- coding: utf-8 -*-
__version__ = '0.7.1' # Time-stamp: <2026-10-22T18:04:40Z>
## Language: Japanese/UTF-8

from setuptools import setup, find_packages

setup(name='sympy_matrix_tools',
      version=__version__,
      description='Some tools for sympy matrices.',
      author='JRF',
      url='https://github.com/JRF-2018/sympy_matrix_tools',
      packages=find_packages(exclude=['tests']),
      install_requires=['sympy>=1.14'],)
