#!/bin/bash
#
# Uses pyinstaller to build a (heavy) standalone binary of py_interp.
# 
# See http://www.pyinstaller.org/
#
pyinstaller -n ./py_interp.exe -F --hidden-import=netCDF4_utils --hidden-import=netcdftime py_interp.py
