# py_interp

Authors: Markel García-Díez (1), Carlos Blanco (2)

1. Institut Català de Ciències del Clima, Barcelona, Spain
2. Department of Applied Mathematics and CC, University of Cantabria, Santander, Spain

Thanks to:  Cindy Bruyere, Lluis Fita and Jesús Fernandez 

## Description

py_interp is a command line tool for interpolating WRF output files to pressure levels, as well as computing a series of diagnostics such as low, medium and high cloud cover. py_interp is a python version of p_interp, originally written in fortran, (see in http://www2.mmm.ucar.edu/wrf/users/). It uses the f2py interface to call fortran subroutines extracted from the original p_interp.F90. 

The reason of writing py_interp is that python enables to implement the same functionality in much less lines than fortran, with less redundancies. Thanks to this, is way easier to add and new diagnostics and to debug them. Thanks to f2py, there is not a noticeable loss of performance from the pure fortran version, as the heavy loops are carried out in fortran too. These fortran routines, taken from py_interp, are located in src/py_interp_fortran.F90 inside a module call "routines".

## Install

py_interp can be installed downloading the source and typing python setup.py install, but the recommended way is to use PyPI.

## Usage

Please type python py_interp.py -h for help.

## Standalone executable

It is possible to make py_interp to be a standalone execulable, though it is heavy (14M), using pyinstaller. Use make_executable.sh to build it.
