from setuptools import find_packages

from numpy.distutils.core import setup, Extension

setup(name="map_res",
      version="0.0.1",
    sources=['src/map_res/test_interp2.f90','src/map_res/mapit.f90'],
    f2py_options=['--quiet'],
      ext_modules=[])