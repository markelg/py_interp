#!/usr/bin/env python

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='py_interp',
    platforms=['GNU/Linux'],
    version='1.0.0',
    author='Marker Garcia, Carlos Blanco',
    keywords=['climate', 'weather'],
    install_requires=[ 'numpy', 'netCDF4'],
    url='https://github.com/markelg/py_interp',
    packages=['py_interp'],
    license='MIT',
    include_package_data=True,
    package_data={'py_interp': [
        'README',
    ]
    },
    scripts=['bin/py_interp.py'], 
)

from numpy.distutils.core import Extension, setup
# Fortran extension
setup(
  ext_modules = [Extension( 'py_interp_fortran', ['src/py_interp_fortran.F90'] )],
)
