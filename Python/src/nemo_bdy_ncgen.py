##
# Creates Nemo Bdy netCDF file ready for population
##
# Written by John Kazimierz Farey, started August 30, 2012
# Port of Matlab code of James Harle
##
from netCDF4 import Dataset 
import datetime
import logging

def CreateBDYNetcdfFile(filename, N,I,J,K,rw,h,orig,fv,grd):
    dimensionNames = ['T', 'I','U','V','E']
    
    # Dimension Lengths
    xb_len = N
    yb_len = 1
    x_len  = I
    y_len  = J
    depth_len = K
    
    # Enter define mode
    ncid = Dataset(filename, 'w', clobber=True, format='NETCDF4')
    
    #define dimensions
    if grd in dimensionNames:
        dimztID = ncid.createDimension('z', xb_len)
    dimxbID = ncid.createDimension('xb',xb_len)
    dimybID = ncid.createDimension('yb',yb_len)
    dimxID  = ncid.createDimension('x', x_len)
    dimyID  = ncid.createDimension('y', y_len)
    dimtcID = ncid.createDimension('time_counter',None)

    #define variable  
    vartcID = ncid.createVariable('time_counter','f4',('time_counter',))
    varlonID = ncid.createVariable('nav_lon','f4',('x','y',))
    varlatID = ncid.createVariable('nav_lat','f4',('x','y',))
    

    if grd in ['E'] :
        varztID = ncid.createVariable('deptht','f4',('xb','yb','z',))
        varmskID = ncid.createVariable('bdy_msk','f4',('x','y',),fill_value=fv)
        varN1pID = ncid.createVariable('N1p','f4',('xb','yb','z','time_counter',),fill_value=fv)
        varN3nID = ncid.createVariable('N3n','f4',('xb','yb','z','time_counter',),fill_value=fv)
        varN5sID = ncid.createVariable('N5s','f4',('xb','yb','z','time_counter',),fill_value=fv)
    elif grd in ['T','I'] :
        varztID = ncid.createVariable('deptht','f4',('xb','yb','z',))
        varmskID = ncid.createVariable('bdy_msk','f4',('x','y',),fill_value=fv)
        vartmpID = ncid.createVariable('votemper','f4',('xb','yb','z','time_counter',),fill_value=fv)
        varsalID = ncid.createVariable('vosaline','f4',('xb','yb','z','time_counter',),fill_value=fv)
        if grd == 'I' :
            varildID = ncid.createVariable('ileadfra','f4',('xb','yb','time_counter',),fill_value=fv)
            variicID = ncid.createVariable('iicethic','f4',('xb','yb','time_counter',),fill_value=fv)
            varisnID = ncid.createVariable('isnowthi','f4',('xb','yb','time_counter',),fill_value=fv)
    elif grd == 'U':
        varztID = ncid.createVariable('depthu','f4',('xb','yb','z',),fill_value=fv)
        varbtuID = ncid.createVariable('vobtcrtx','f4',('xb','yb','time_counter',),fill_value=fv)
        vartouID = ncid.createVariable('vozocrtx','f4',('xb','yb','z','time_counter',),fill_value=fv)
    elif grd == 'V':
        varztID = ncid.createVariable('depthv','f4',('xb','yb','z',))
        varbtvID = ncid.createVariable('vobtcrty','f4',('xb','yb','time_counter',),fill_value=fv)
        vartovID = ncid.createVariable('vomecrty','f4',('xb','yb','z','time_counter',),fill_value=fv)
    elif grd == 'Z':
        varsshID = ncid.createVariable('sossheig','f4',('xb','yb','time_counter',),fill_value=fv)
    else :
        logging.error("Unknow Grid input")
        
    
    varnbiID = ncid.createVariable('nbidta','i4',('xb','yb',))
    varnbjID = ncid.createVariable('nbjdta','i4',('xb','yb',))
    varnbrID = ncid.createVariable('nbrdta','i4',('xb','yb',))
    #Global Attributes
    ncid.file_name = filename
    ncid.creation_date = str(datetime.datetime.now())
    ncid.rim_width = rw
    ncid.history = h
    ncid.institution = 'National Geography Centre, Liverpool, UK'
    
    #Time axis attributes
    vartcID.axis = 'T'
    vartcID.standard_name = 'time'
    vartcID.units = 'seconds since '+orig
    vartcID.title = 'Time'
    vartcID.long_name = 'Time axis'
    vartcID.time_origin = orig
    
    #Longitude axis attributes
    varlonID.axis = 'Longitude'
    varlonID.short_name = 'nav_lon'
    varlonID.units = 'degrees_east'
    varlonID.long_name = 'Longitude'
    
    #Latitude axis attributes
    varlatID.axis = 'Latitude'
    varlatID.short_name = 'nav_lat'
    varlatID.units = 'degrees_east'
    varlatID.long_name = 'Latitude'
    
    #nbidta attributes
    varnbiID.short_name = 'nbidta'
    varnbiID.units = 'unitless'
    varnbiID.long_name = 'Bdy i indices'
    
    #nbjdta attributes
    varnbjID.short_name = 'nbjdta'
    varnbjID.units = 'unitless'
    varnbjID.long_name = 'Bdy j indices'
    
    #nbrdta attributes
    varnbrID.short_name = 'nbrdta'
    varnbrID.units = 'unitless'
    varnbrID.long_name = 'Bdy discrete distance'
    if grd == 'E' :
        varztID.axis = 'Depth'
        varztID.short_name = 'deptht'
        varztID.units = 'm'
        varztID.long_name = 'Depth'
        
        varmskID.short_name = 'bdy_msk'
        varmskID.units = 'unitless'
        varmskID.long_name = 'Unstructed boundary mask'
        
        varN1pID.units = 'mmol/m^3'
        varN1pID.short_name = 'N1p'
        varN1pID.long_name = 'Phosphate'
        varN1pID.grid = 'bdyT'
        
        varN3nID.units = 'mmol/m^3'
        varN3nID.short_name = 'N3n'
        varN3nID.long_name = 'Nitrate'
        varN3nID.grid = 'bdyT'
        
        varN5sID.units = 'mmol/m^3'
        varN5sID.short_name = 'N5s'
        varN5sID.long_name = 'Silicate'
        varN5sID.grid = 'bdyT'
        
    if grd in ['T', 'I'] :
        varztID.axis = 'Depth'
        varztID.short_name = 'deptht'
        varztID.units = 'm'
        varztID.long_name = 'Depth'
        
        varmskID.short_name = 'bdy_msk'
        varmskID.units = 'unitless'
        varmskID.long_name = 'Unstructed boundary mask'
        
        vartmpID.units = 'C'
        vartmpID.short_name = 'votemper'
        vartmpID.long_name = 'Temperature'
        vartmpID.grid = 'bdyT'
        
        varsalID.units = 'PSU'
        varsalID.short_name = 'vosaline'
        varsalID.long_name = 'Salinity'
        varsalID.grid = 'bdyT'
        
        if grd == 'I' :
            varildID.units = '%'
            varildID.short_name = 'ildsconc'
            varildID.long_name = 'Ice lead fraction'
            varildID.grid = 'bdyT'
        
            variicID.units = 'm'
            variicID.short_name = 'iicethic'
            variicID.long_name = 'Ice thickness'
            variicID.grid = 'bdyT'
            
            varisnID.units = 'm'
            varisnID.short_name = 'isnowthi'
            varisnID.long_name = 'Snow thickness'
            varisnID.grid = 'bdyT'
    elif grd == 'U' :
        varztID.axis = 'Depth'
        varztID.short_name = 'depthu'
        varztID.units = 'm'
        varztID.long_name = 'Depth'
        
        varbtuID.units = 'm/s'
        varbtuID.short_name = 'vobtcrtx'
        varbtuID.long_name = 'Barotropic zonal Current'
        varbtuID.grid = 'bdyU'
        
        vartouID.units = 'm/s'
        vartouID.short_name = 'vozocrtx'
        vartouID.long_name = 'Zonal Current'
        vartouID.grid = 'bdyU'
        
    elif grd == 'V':
        varztID.axis = 'Depth'
        varztID.short_name = 'depthv'
        varztID.units = 'm'
        varztID.long_name = 'Depth'
        
        varbtvID.units = 'm/s'
        varbtvID.short_name = 'vobtcrty'
        varbtvID.long_name = 'Barotropic meridional Current'
        varbtvID.grid = 'bdyV'
        
        vartovID.units = 'm/s'
        vartovID.short_name = 'vomecrty'
        vartovID.long_name = 'Meridional Current'
        vartovID.grid = 'bdyV'
        
    elif grd == 'Z':
        varsshID.units = 'm'
        varsshID.short_name = 'sossheig'
        varsshID.long_name = 'Sea Surface Height'
        varsshID.grid = 'bdyT'
        
    else :
        logging.error('Unknown Grid')
        
    ncid.close()

              
    


