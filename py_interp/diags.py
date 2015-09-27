#
# Diagnostics for py_interp
#
import numpy as np
from py_interp_fortran import routines as f90
from py_interp.fun import massvertint
tr = np.transpose

def compute_MSLP(ivar, inc, onc, bf, plevs):
    ovardata = f90.compute_mslp(tr(bf.pres_field), tr(bf.psfc), tr(bf.hgt), tr(bf.temp), tr(bf.qvapor))
    ovarobj = onc.createVariable(ivar, 'float32', ["Time", "south_north", "west_east"], zlib=True, complevel=4, shuffle=True)
    ovarobj[:] = tr(ovardata)
    ovarobj.FieldType  = 104
    ovarobj.MemoryOrder = "Z"
    ovarobj.description = "Pressure levels"
    ovarobj.units = "pa"
    ovarobj.stagger = "-"
    ovarobj.coordinates = "XLONG XLAT"
    return onc

def compute_GHT(ivar, inc, onc, bf, plevs):
    ovardata = f90.interp(tr(bf.ght), tr(bf.pres_field), plevs, tr(bf.psfc), tr(bf.hgt), tr(bf.temp), tr(bf.qvapor),
					linlog=1, extrapolate=1, geopt=True, missing=1.e36)
    ovarobj = onc.createVariable(ivar, 'float32', ["Time", "num_metgrid_levels", "south_north", "west_east"], zlib=True, complevel=4, shuffle=True)
    ovarobj[:] = tr(ovardata)
    ovarobj.FieldType  = 104
    ovarobj.MemoryOrder = "XZY"
    ovarobj.description = "Geopotential Height"
    ovarobj.units = "m"
    ovarobj.stagger = "-"
    ovarobj.coordinates = "XLONG XLAT"
    return onc

def compute_PRES(ivar, inc, onc, bf, plevs):
    nj, ni = bf.pres_field.shape[2], bf.pres_field.shape[3]
    nplev = len(plevs)
    nt = bf.pres_field.shape[0]
    ovardata = np.repeat(np.nan, nt*nplev*nj*ni).reshape(nt, nplev, nj, ni)
    for n in range(nplev):
        ovardata[:, n, :, :] = plevs[n]
    ovarobj = onc.createVariable(ivar, 'float32', ["Time", "num_metgrid_levels", "south_north", "west_east"], zlib=True, complevel=4, shuffle=True)  
    ovarobj[:] = ovardata        
    ovarobj.FieldType = 104
    ovarobj.MemoryOrder = "XZY"
    ovarobj.description = "Pressure"
    ovarobj.units = "Pa"
    ovarobj.stagger = "-"
    ovarobj.coordinates = "XLONG XLAT"
    return onc
        
def compute_TT(ivar, inc, onc, bf, plevs):
    ovardata = f90.interp(tr(bf.temp), tr(bf.pres_field), plevs, tr(bf.psfc), tr(bf.hgt), tr(bf.temp), tr(bf.qvapor),
            linlog=1, extrapolate=1, geopt=False, missing=1.e36)
    ovarobj = onc.createVariable(ivar, 'float32', ["Time", "num_metgrid_levels", "south_north", "west_east"], zlib=True, complevel=4, shuffle=True)
    ovarobj[:] = tr(ovardata)
    ovarobj.FieldType = 104
    ovarobj.MemoryOrder = "XZY"
    ovarobj.description = "Temperature"
    ovarobj.units = "K"
    ovarobj.stagger = "-"
    ovarobj.coordinates = "XLONG XLAT"
    return onc
        
def compute_RH(ivar, inc, onc, bf, plevs):
    ovardata = f90.interp(tr(bf.rh), tr(bf.pres_field), plevs, tr(bf.psfc), tr(bf.hgt), tr(bf.temp), tr(bf.qvapor),
            linlog=1, extrapolate=1, geopt=False, missing=1.e36)
    ovarobj = onc.createVariable(ivar, 'float32', ["Time", "num_metgrid_levels", "south_north", "west_east"], zlib=True, complevel=4, shuffle=True)  
    ovarobj[:] = tr(ovardata)
    ovarobj.FieldType = 104 ;
    ovarobj.MemoryOrder = "XZY" ;
    ovarobj.description = "Relative Humidity" ;
    ovarobj.units = "%" ;
    ovarobj.stagger = "-" ;
    ovarobj.coordinates = "XLONG XLAT" ;
    return onc

def compute_CLT_OLD(ivar, inc, onc, bf, plevs):
    cldfra = inc.variables["CLDFRA"][:]
    ovardata = f90.clt_sundqvist(tr(cldfra))
    ovarobj = onc.createVariable(ivar, 'float32', ["Time", "south_north", "west_east"], zlib=True, complevel=4, shuffle=True)
    ovarobj[:] =  tr(ovardata)
    ovarobj.FieldType = 104
    ovarobj.MemoryOrder = "XY"
    ovarobj.description = "total cloud fraction Sundqvist"
    ovarobj.units = "1"
    ovarobj.stagger = ""
    ovarobj.coordinates = "XLONG XLAT"
    return onc
        
def compute_CLT(ivar, inc, onc, bf, plevs):
    cldfra = inc.variables["CLDFRA"][:]
    ovardata = f90.clt_maxrand(tr(cldfra))
    ovarobj = onc.createVariable(ivar, 'float32', ["Time", "south_north", "west_east"], zlib=True, complevel=4, shuffle=True)
    ovarobj[:] =  tr(ovardata)
    ovarobj.FieldType = 104
    ovarobj.MemoryOrder = "XY"
    ovarobj.description = "total cloud fraction"
    ovarobj.units = "1"
    ovarobj.stagger = ""
    ovarobj.coordinates = "XLONG XLAT"
    return onc
        
def compute_CLL(ivar, inc, onc, bf, plevs):
    cldfra = inc.variables["CLDFRA"][:]
    ovardata = f90.clt_maxrand_levels(tr(cldfra), tr(bf.pres_field), maxpres=1000, minpres=680)
    ovarobj = onc.createVariable(ivar, 'float32', ["Time", "south_north", "west_east"], zlib=True, complevel=4, shuffle=True)
    ovarobj[:] =  tr(ovardata)
    ovarobj.FieldType = 104
    ovarobj.MemoryOrder = "XY"
    ovarobj.description = "Low cloud fraction (1000-680 hPa)"
    ovarobj.units = "1"
    ovarobj.stagger = ""
    ovarobj.coordinates = "XLONG XLAT"
    return onc
        
def compute_CLM(ivar, inc, onc, bf, plevs):
    cldfra = inc.variables["CLDFRA"][:]
    ovardata = f90.clt_maxrand_levels(tr(cldfra), tr(bf.pres_field), maxpres=680, minpres=440)
    ovarobj = onc.createVariable(ivar, 'float32', ["Time", "south_north", "west_east"], zlib=True, complevel=4, shuffle=True)
    ovarobj[:] =  tr(ovardata)
    ovarobj.FieldType = 104
    ovarobj.MemoryOrder = "XY"
    ovarobj.description = "Medium cloud fraction (680-440 hPa)"
    ovarobj.units = "1"
    ovarobj.stagger = ""
    ovarobj.coordinates = "XLONG XLAT"
    return onc
        
def compute_CLH(ivar, inc, onc, bf, plevs):
    cldfra = inc.variables["CLDFRA"][:]
    ovardata = f90.clt_maxrand_levels(tr(cldfra), tr(bf.pres_field), maxpres=440, minpres=10)
    ovarobj = onc.createVariable(ivar, 'float32', ["Time", "south_north", "west_east"], zlib=True, complevel=4, shuffle=True)
    ovarobj[:] =  tr(ovardata)
    ovarobj.FieldType = 104
    ovarobj.MemoryOrder = "XY"
    ovarobj.description = "High cloud fraction (440-10 hPa)"
    ovarobj.units = "1"
    ovarobj.stagger = ""
    ovarobj.coordinates = "XLONG XLAT"
    return onc

def compute_VIM(ivar, inc, onc, bf, plevs):
    iarr = inc.variables["QVAPOR"][:]
    ovardata =  massvertint(iarr, inc)
    ovarobj = onc.createVariable(ivar, 'float32', ["Time", "south_north", "west_east"], zlib=True, complevel=4, shuffle=True)
    ovarobj[:] =  ovardata
    ovarobj.FieldType = 104
    ovarobj.MemoryOrder = "XY"
    ovarobj.description = "Vertically integrated moisture"
    ovarobj.units = "Kg m-2"
    ovarobj.stagger = ""
    ovarobj.coordinates = "XLONG XLAT" 
    return onc

def compute_VIQC(ivar, inc, onc, bf, plevs):
    iarr = inc.variables["QCLOUD"][:]
    ovardata =  massvertint(iarr, inc)
    ovarobj = onc.createVariable(ivar, 'float32', ["Time", "south_north", "west_east"], zlib=True, complevel=4, shuffle=True)
    ovarobj[:] =  ovardata
    ovarobj.FieldType = 104
    ovarobj.MemoryOrder = "XY"
    ovarobj.description = "Vertically integrated cloud water"
    ovarobj.units = "Kg m-2"
    ovarobj.stagger = ""
    ovarobj.coordinates = "XLONG XLAT" 
    return onc

def compute_VIQI(ivar, inc, onc, bf, plevs):
    iarr = inc.variables["QICE"][:]
    ovardata =  massvertint(iarr, inc)
    ovarobj = onc.createVariable(ivar, 'float32', ["Time", "south_north", "west_east"], zlib=True, complevel=4, shuffle=True)
    ovarobj[:] =  ovardata
    ovarobj.FieldType = 104
    ovarobj.MemoryOrder = "XY"
    ovarobj.description = "Vertically integrated cloud ice"
    ovarobj.units = "Kg m-2"
    ovarobj.stagger = ""
    ovarobj.coordinates = "XLONG XLAT" 
    return onc
