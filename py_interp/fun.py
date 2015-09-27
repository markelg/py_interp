#
# Library with utilities to maniputale netCDF files
#
import logging
import numpy   as np
import netCDF4 as ncdf
import sys, time
from py_interp_fortran import routines as f90

tr = np.transpose # Shorter way to call it
#
# Function to copy netCDF structures.
#
def copy_netcdf_structure(ifile, ofile, variables, dimensions=None, isncobj = False,
    xydims = ["x", "y"], oformat = None, del_gattr=False, limtime2unlim=False):
    logging.debug( "Creating %s netCDF file" % (ofile) )
    if isncobj:
        inc = ifile
    else:
        inc = ncdf.Dataset(ifile,'r')
    if not oformat: oformat = inc.file_format
    onc = ncdf.Dataset(ofile, "w", format = oformat)
    if variables == "all":
        variables = inc.variables.keys()
    # 
    # Copy global attributes and redefine history
    #
    if del_gattr == False:
        for ikey, ivalue in inc.__dict__.items():
            onc.setncattr(ikey, ivalue)
            onc.history = "Created by %s on %s" % (sys.argv[0],time.ctime(time.time()))
    #
    # Copy dimensions
    #
    for dimname, dimobj in inc.dimensions.items():
        if dimensions:
            if dimname not in dimensions:
                continue
        logging.debug( "Setting dimension %s %s" % (dimname, dimobj) )
        if limtime2unlim and dimname == "time":
            onc.createDimension("time", None)
            continue
        if dimobj.isunlimited():
            logging.debug( "Dimension is unlimited" )
            onc.createDimension(dimname, None)
        else:
            onc.createDimension(dimname, len(dimobj))
    #
    # Copy variables specified in the argument
    #
    for ivarname in variables:
        ivarobj = inc.variables[ivarname]
        if len(ivarobj.dimensions) > 2:
            ovarobj = onc.createVariable(ivarname, ivarobj.dtype, ivarobj.dimensions, zlib=True, complevel=4, shuffle=True)
        else:
            ovarobj = onc.createVariable(ivarname, ivarobj.dtype, ivarobj.dimensions)
        for attrname, attrvalue in ivarobj.__dict__.items():
            ovarobj.setncattr(attrname, attrvalue)
        ovarobj[:] = ivarobj[:]
    onc.sync()
    return onc

def copy_n_filter_wrfout(inc, ofile, copyvars):
    #
    # Check if we need the soil layers dimension
    #
    if ("TSLB" in copyvars) or (("SMOIS" in copyvars) or ("SH2O" in copyvars)):
        dims = ["Time", "south_north", "west_east", "DateStrLen", "soil_layers_stag"]
    else:
        dims = ["Time", "south_north", "west_east", "DateStrLen"]
    #
    # Copy the input wrfout
    #
    onc = copy_netcdf_structure(inc, ofile, variables=copyvars, oformat='NETCDF4_CLASSIC',
        dimensions=dims, isncobj = True, xydims = ["west_east", "south_north"])
    return onc

def add_pressure_axis(onc, plevs):
    onc.createDimension("num_metgrid_levels", len(plevs))
    plev = onc.createVariable("PLEV", 'float32', ["num_metgrid_levels"])
    plev[:] = plevs
    plev.FieldType  =104
    plev.MemoryOrder = "Z"
    plev.description = "Pressure levels"
    plev.units = "pa"
    plev.stagger = "-"
    plev.coordinates = "XLONG XLAT"
    onc.sync()
    return onc	

def is_staggered(varobj):
    if varobj.stagger != "":
        return True
    else:
        return False

def de_stagger(ivarobj, ivardata):
    dim_staggered = ivarobj.stagger
    if dim_staggered == "X":
        ovardata = 0.5*(ivardata[:, :, :, 0:-1] + ivardata[:, :, :, 1:]) 
    elif dim_staggered == "Y":
        ovardata = 0.5*(ivardata[:, :, 0:-1, :] + ivardata[:, :, 1:, :]) 
    elif dim_staggered == "Z":
        ovardata = 0.5*(ivardata[:, 0:-1, :, :] + ivardata[:, 1:, :, :])
    return ovardata
    
def interp2plevs(ivar, inc, onc, bf, plevs):
    ivarobj = inc.variables[ivar]
    ivardata = ivarobj[:]
    #
    # Check if the variable is staggered and de-stagger it if necessary
    #
    if is_staggered(ivarobj):
        logging.debug( "Variable %s is staggered, de-staggering" % ivar )
        ivardata = de_stagger(ivarobj, ivardata)
    #
    # Call fortran interpolation routine
    #
    logging.debug( "Calling interpolation routine" )
    ovardata = f90.interp(tr(ivardata), tr(bf.pres_field), plevs, tr(bf.psfc), tr(bf.hgt), tr(bf.temp), tr(bf.qvapor),
                        linlog=1, extrapolate=1, geopt=False, missing=1.e36)
    #
    # Create the output variable and add data and attributes
    #
    ovarobj = onc.createVariable(ivar, ivarobj.dtype, ["Time", "num_metgrid_levels", "south_north", "west_east"], zlib=True, complevel=4, shuffle=True)
    ovarobj[:] = tr(ovardata)
    for attrname, attrvalue in ivarobj.__dict__.items():
        if attrname == "stagger":
            ovarobj.setncattr(attrname, "")
            continue
        ovarobj.setncattr(attrname, attrvalue)
    return onc
#
# Function to compute mass-weighted vertical integrals: Input, an array and a netcdf object.
# As WRF uses an hybrid sigma vertical coordinate, the distance between levels 
# ,measured in the eta coordinate, is proportional to the mass differential:
# dm = -(colmass/g)*deta Where colmass is the total
# weight of the air columns and g 9.8 m2s-2
#
def massvertint(iarr, inc):
    deta = inc.variables["DNW"][0, :]
    #
    # Ugly and unefficient, but numpy.tile and numpy.resize are giving weird
    # rounding problems. After 1 hour trying to fix it I surrender...
    #
    NT, NZ, NJ, NI = iarr.shape
    deta_3d = np.repeat(np.nan, NZ*NJ*NI).reshape((NZ, NJ, NI))
    for j in range(NJ):
        for i in range(NI):
            deta_3d[:, j, i] = deta[:]
    colmass = inc.variables["MU"][:] + inc.variables["MUB"][:]
    intarr = np.sum(iarr*deta_3d, axis=1)
    intarr = -colmass*intarr*(1./(9.8))
    return intarr
    
class BasicFields:
    def __init__(self, inc):
        #
        # Basic constants
        #
        self.Rd = 287.04   # Gas constant for the air (J K-1 Kg-1)
        self.Cp = 7.*self.Rd/2. #  Specific heat (J Kg-1)
        self.RCP = self.Rd/self.Cp
        self.P0 = 100000 # Base pressure for the Poisson equation in Pa
        logging.debug( "Reading basic fields..." )
        self.get_sfc_pressure(inc)
        self.get_terrain_height(inc)
        self.get_pressure(inc)
        self.get_geopotential(inc)
        self.get_temperature(inc)
        self.get_qvapor(inc)
        self.get_relative_humidity(inc)
    #
    # Methods for reading basic fields
    #
    def get_sfc_pressure(self, inc):
        # in Pa
        self.psfc = inc.variables["PSFC"][:]
    def get_terrain_height(self, inc):
        # in m
        self.hgt = inc.variables["HGT"][0, :, :]
    def get_pressure(self, inc):
        # in Pa
        base_pres = inc.variables["PB"][:]
        perturbation_pres =  inc.variables["P"][:]
        self.pres_field = base_pres + perturbation_pres
    def get_geopotential(self, inc):
        # in m2s-2
        base_ght = inc.variables["PHB"][:]
        perturbation_ght =  inc.variables["PH"][:]
        ght = base_ght + perturbation_ght
        # de-stagger
        self.ght = (ght[:, 0:-1, :, :] + ght[:, 1:, :, :])/2
    def get_temperature(self, inc):
        # in k
        base_pot_temp = 300 # not T0 which is other thing!!
        perturbation_pot_temp = inc.variables["T"][:]
        pot_temp = base_pot_temp + perturbation_pot_temp
        self.temp = pot_temp*(self.pres_field/self.P0)**self.RCP
    def get_qvapor(self, inc):
        # in k
        self.qvapor = inc.variables["QVAPOR"][:]
    def get_relative_humidity(self, inc):
        # in 100/1
        data1 = 10.*0.6112*np.exp(17.67*(self.temp - 273.16)/(self.temp - 29.65))
        data2 = 0.622*data1/(0.01 * self.pres_field -  (1.-0.622)*data1)
        self.rh = 100.*self.qvapor/data2
        self.rh = np.where(self.rh > 100., 100., self.rh) 
