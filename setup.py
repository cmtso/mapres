"""
https://github.com/cmtso/map_res

You must cite:
Timothy C. Johnson, Glenn E. Hammond, Xingyuan Chen,
PFLOTRAN-E4D: A parallel open source PFLOTRAN module for simulating time-lapse electrical resistivity data,
Computers & Geosciences,Volume 99,2017,Pages 72-80,https://doi.org/10.1016/j.cageo.2016.09.006
"""

from setuptools import find_packages

from numpy.distutils.core import setup, Extension

setup(name="mapres",
      version="0.0.1",
    sources=['src/mapres/test_interp2.f90','src/mapres/mapit.f90'],
    f2py_options=['--quiet'],
      ext_modules=[])
