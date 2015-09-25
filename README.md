# py_interp

Authors: Markel García-Díez (1), Carlos Blanco (2)

1 Institut Català de Ciències del Clima, Barcelona, Spain
2 Department of Applied Mathematics and CC, University of Cantabria, Santander, Spain

## Description

py_interp is a command line tool for interpolating WRF output files to pressure levels, as well as computing a series of diagnostics such as low, medium and high cloud cover. py_interp is a version of p_interp (see in http://www2.mmm.ucar.edu/wrf/users/) written in python. It uses the f2py interface to call fortran subroutines extracted from the original p_interp. These are located in src/py_interp_fortran.F90 inside a module call "routines".

## Install

py_interp can be installed downloading the source and typing python setup.py install, but the recommended way is to use PyPI.

## Usage

Please type python py_interp.py -h for help.

## Standalone executable

It is possible to make py_interp to be a standalone execulable, though it is heavy (14M), using pyinstaller. Use make_executable.sh to build it.
