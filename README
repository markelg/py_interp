>> Markel Garcia Diez 8 Feb 2013

py_interp is a version of p_interp written in python. It uses the f2py interface
to call fortran subroutines extracted from the original p_interp. These are located
in  py_interp_fortran.F90 inside a module call "routines". To make it a python module
type:

f2py -c py_interp_fortran.F90 -m py_interp_fortran

Check https://sysbio.ioc.ee/projects/f2py2e/ for more info about f2py.

Once compiled, it can be imported into python as a regular python package.

from py_interp_fortran import routines as f90

Usage: type python py_interp.py -h for help.

>> Markel Garcia Diez 11 Feb 2014

Now it is possible to make py_interp to be a standalone execulable, though it is heavy (14M), using pyinstaller. Use make_executable.sh to build it.

>> Carlos Blanco 23 Sep 2015 

py_interp is compatible with 2.x and 3.x

>> Carlos Blanco 23 Sep 2015

Created a setup.py file for py_interp
