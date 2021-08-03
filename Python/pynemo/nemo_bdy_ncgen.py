'''
Creates Nemo Bdy netCDF file ready for population

Written by John Kazimierz Farey, started August 30, 2012
Port of Matlab code of James Harle
'''
# pylint: disable=E1103
# pylint: disable=no-name-in-module
from netCDF4 import Dataset
import datetime
import logging

def CreateBDYNetcdfFile(filename, N, I, J, K, rw, h, orig, fv, calendar, grd):
    """ This method creates a template of bdy netcdf files. A common for
    T, I, U, V, E grid types.

    Keyword arguments:
    filename -- name of output bdy file name.
    N -- Number of border points.
    I -- Grid width.
    J -- Grid Height.
    K -- Grid Depth.
    rw -- Rim Width.
    h -- Meta data information.
    orig -- unit origin (calendar).
    fv -- Fill Value.
    calendar -- time calendar.
    grd -- grid types. ['T'=3D Temperature/Salinity,
                        'I'=3D Temperature/Salinity + 2D ice,
                        'U'=3D u-vel,
                        'V'=3D v-vel,
                        'E'=3D biogeochem (phosphate, nitrate, silocate),
                        'Z'=2D sea surface height]
    """

    gridNames = ['T', 'I', 'U', 'V', 'E', 'Z'] # All possible grids

    # Dimension Lengths
    xb_len = N
    yb_len = 1
    x_len = I
    y_len = J
    depth_len = K

    # Enter define mode
    ncid = Dataset(filename, 'w', clobber=True, format='NETCDF4')

    #define dimensions
    if grd in gridNames and grd != 'Z': # i.e grid NOT barotropic (Z) 
        dimztID = ncid.createDimension('z', depth_len)
    else:
        logging.error('Grid tpye not known')
    dimxbID = ncid.createDimension('xb', xb_len)
    dimybID = ncid.createDimension('yb', yb_len)
    dimxID = ncid.createDimension('x', x_len)
    dimyID = ncid.createDimension('y', y_len)
    dimtcID = ncid.createDimension('time_counter', None)

    #define variable
    vartcID = ncid.createVariable('time_counter', 'f4', ('time_counter', ))
    varlonID = ncid.createVariable('nav_lon', 'f4', ('y', 'x', ))
    varlatID = ncid.createVariable('nav_lat', 'f4', ('y', 'x', ))

    if grd in ['E']:
        varztID = ncid.createVariable('gdept', 'f4', ('z', 'yb', 'xb', ))
        vardzID = ncid.createVariable('e3t', 'f4', ('z', 'yb', 'xb', ))
        varmskID = ncid.createVariable('bdy_msk', 'f4', ('y', 'x', ), fill_value=fv)
        varN1pID = ncid.createVariable('N1p', 'f4', ('time_counter', 'z', 'yb', 'xb', ),
                                       fill_value=fv)
        varN3nID = ncid.createVariable('N3n', 'f4', ('time_counter', 'z', 'yb', 'xb', ),
                                       fill_value=fv)
        varN5sID = ncid.createVariable('N5s', 'f4', ('time_counter', 'z', 'yb', 'xb', ),
                                       fill_value=fv)
    elif grd in ['T', 'I']:
        varztID = ncid.createVariable('gdept', 'f4', ('z', 'yb', 'xb', ))
        vardzID = ncid.createVariable('e3t', 'f4', ('z', 'yb', 'xb', ))
        varmskID = ncid.createVariable('bdy_msk', 'f4', ('y', 'x', ), fill_value=fv)
        vartmpID = ncid.createVariable('votemper', 'f4', ('time_counter', 'z', 'yb', 'xb', ),
                                       fill_value=fv)
        varsalID = ncid.createVariable('vosaline', 'f4', ('time_counter', 'z', 'yb', 'xb', ),
                                       fill_value=fv)
        varchnID = ncid.createVariable('CHN', 'f4', ('time_counter', 'z', 'yb', 'xb', ),
                                       fill_value=fv)
        varchdID = ncid.createVariable('CHD', 'f4', ('time_counter', 'z', 'yb', 'xb', ),
                                       fill_value=fv)
        varphnID = ncid.createVariable('PHN', 'f4', ('time_counter', 'z', 'yb', 'xb', ),
                                       fill_value=fv)
        varphdID = ncid.createVariable('PHD', 'f4', ('time_counter', 'z', 'yb', 'xb', ),
                                       fill_value=fv)
        varzmiID = ncid.createVariable('ZMI', 'f4', ('time_counter', 'z', 'yb', 'xb', ),
                                       fill_value=fv)
        varzmeID = ncid.createVariable('ZME', 'f4', ('time_counter', 'z', 'yb', 'xb', ),
                                       fill_value=fv)
        vardinID = ncid.createVariable('DIN', 'f4', ('time_counter', 'z', 'yb', 'xb', ),
                                       fill_value=fv)
        varsilID = ncid.createVariable('SIL', 'f4', ('time_counter', 'z', 'yb', 'xb', ),
                                       fill_value=fv)
        varferID = ncid.createVariable('FER', 'f4', ('time_counter', 'z', 'yb', 'xb', ),
                                       fill_value=fv)
        vardetID = ncid.createVariable('DET', 'f4', ('time_counter', 'z', 'yb', 'xb', ),
                                       fill_value=fv)
        varpdsID = ncid.createVariable('PDS', 'f4', ('time_counter', 'z', 'yb', 'xb', ),
                                       fill_value=fv)
        vardtcID = ncid.createVariable('DTC', 'f4', ('time_counter', 'z', 'yb', 'xb', ),
                                       fill_value=fv)
        vardicID = ncid.createVariable('DIC', 'f4', ('time_counter', 'z', 'yb', 'xb', ),
                                       fill_value=fv)
        varalkID = ncid.createVariable('ALK', 'f4', ('time_counter', 'z', 'yb', 'xb', ),
                                       fill_value=fv)
        varoxyID = ncid.createVariable('OXY', 'f4', ('time_counter', 'z', 'yb', 'xb', ),
                                       fill_value=fv)
        if grd == 'I':
            varildID = ncid.createVariable('ileadfra', 'f4', ('time_counter', 'yb', 'xb',),
                                           fill_value=fv)
            variicID = ncid.createVariable('iicethic', 'f4', ('time_counter', 'yb', 'xb',),
                                           fill_value=fv)
            varisnID = ncid.createVariable('isnowthi', 'f4', ('time_counter', 'yb', 'xb',),
                                           fill_value=fv)
    elif grd == 'U':
        varztID = ncid.createVariable('gdepu', 'f4', ('z', 'yb', 'xb', ), fill_value=fv)
        vardzID = ncid.createVariable('e3u', 'f4', ('z', 'yb', 'xb', ), fill_value=fv)
        varbtuID = ncid.createVariable('vobtcrtx', 'f4', ('time_counter', 'yb', 'xb', ),
                                       fill_value=fv)
        vartouID = ncid.createVariable('vozocrtx', 'f4', ('time_counter', 'z', 'yb', 'xb', ),
                                       fill_value=fv)
    elif grd == 'V':
        varztID = ncid.createVariable('gdepv', 'f4', ('z', 'yb', 'xb', ))
        vardzID = ncid.createVariable('e3v', 'f4', ('z', 'yb', 'xb', ))
        varbtvID = ncid.createVariable('vobtcrty', 'f4', ('time_counter', 'yb', 'xb', ),
                                       fill_value=fv)
        vartovID = ncid.createVariable('vomecrty', 'f4', ('time_counter', 'z', 'yb', 'xb',),
                                       fill_value=fv)
    elif grd == 'Z':
        varsshID = ncid.createVariable('sossheig', 'f4', ('time_counter', 'yb', 'xb', ),
                                       fill_value=fv)
        varmskID = ncid.createVariable('bdy_msk', 'f4', ('y', 'x', ), fill_value=fv)
    else:
        logging.error("Unknow Grid input")


    varnbiID = ncid.createVariable('nbidta', 'i4', ('yb', 'xb', ))
    varnbjID = ncid.createVariable('nbjdta', 'i4', ('yb', 'xb', ))
    varnbrID = ncid.createVariable('nbrdta', 'i4', ('yb', 'xb', ))
    #Global Attributes
    ncid.file_name = filename
    ncid.creation_date = str(datetime.datetime.now())
    ncid.rim_width = rw
    ncid.history = h
    ncid.institution = 'National Oceanography Centre, Livepool, U.K.'

    #Time axis attributes
    vartcID.axis = 'T'
    vartcID.standard_name = 'time'
    vartcID.units = 'seconds since '+orig
    vartcID.title = 'Time'
    vartcID.long_name = 'Time axis'
    vartcID.time_origin = orig
    vartcID.calendar = calendar

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
    if grd == 'E':
        varztID.axis = 'Depth'
        varztID.short_name = 'gdept'
        varztID.units = 'm'
        varztID.long_name = 'Depth'

        vardzID.axis = 'Depth'
        vardzID.short_name = 'e3t'
        vardzID.units = 'm'
        vardzID.long_name = 'Depth'

        varmskID.short_name = 'bdy_msk'
        varmskID.units = 'unitless'
        varmskID.long_name = 'Structured boundary mask'

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

    if grd in ['T', 'I']:
        varztID.axis = 'Depth'
        varztID.short_name = 'gdept'
        varztID.units = 'm'
        varztID.long_name = 'Depth'

        vardzID.axis = 'Depth'
        vardzID.short_name = 'e3t'
        vardzID.units = 'm'
        vardzID.long_name = 'Depth'

        varmskID.short_name = 'bdy_msk'
        varmskID.units = 'unitless'
        varmskID.long_name = 'Structured boundary mask'

        vartmpID.units = 'C'
        vartmpID.short_name = 'votemper'
        vartmpID.long_name = 'Temperature'
        vartmpID.grid = 'bdyT'

        varsalID.units = 'PSU'
        varsalID.short_name = 'vosaline'
        varsalID.long_name = 'Salinity'
        varsalID.grid = 'bdyT'

        varchnID.units = 'mg Chl/m3'
        varchnID.short_name = 'CHN'
        varchnID.long_name = 'CHN'
        varchnID.grid = 'bdyT'

        varchdID.units = 'mg Chl/m3'
        varchdID.short_name = 'CHD'
        varchdID.long_name = 'CHD'
        varchdID.grid = 'bdyT'

        varphnID.units = 'mmolN/m3'
        varphnID.short_name = 'PHN'
        varphnID.long_name = 'PHN'
        varphnID.grid = 'bdyT'

        varphdID.units = 'mmolN/m3'
        varphdID.short_name = 'PHD'
        varphdID.long_name = 'PHD'
        varphdID.grid = 'bdyT'

        varzmiID.units = 'mmolN/m3'
        varzmiID.short_name = 'ZMI'
        varzmiID.long_name = 'ZMI'
        varzmiID.grid = 'bdyT'

        varzmeID.units = 'mmolN/m3'
        varzmeID.short_name = 'ZME'
        varzmeID.long_name = 'ZME'
        varzmeID.grid = 'bdyT'

        vardinID.units = 'mmolN/m3'
        vardinID.short_name = 'DIN'
        vardinID.long_name = 'DIN'
        vardinID.grid = 'bdyT'

        varsilID.units = 'mmolSi/m3'
        varsilID.short_name = 'SIL'
        varsilID.long_name = 'SIL'
        varsilID.grid = 'bdyT'

        varferID.units = 'mmolFe/m3'
        varferID.short_name = 'FER'
        varferID.long_name = 'FER'
        varferID.grid = 'bdyT'

        vardetID.units = 'mmolN/m3'
        vardetID.short_name = 'DET'
        vardetID.long_name = 'DET'
        vardetID.grid = 'bdyT'

        varpdsID.units = 'mmolSi/m3'
        varpdsID.short_name = 'PDS'
        varpdsID.long_name = 'PDS'
        varpdsID.grid = 'bdyT'

        vardtcID.units = 'mmolC/m3'
        vardtcID.short_name = 'DTC'
        vardtcID.long_name = 'DTC'
        vardtcID.grid = 'bdyT'

        varalkID.units = 'meq/m3'
        varalkID.short_name = 'ALK'
        varalkID.long_name = 'ALK'
        varalkID.grid = 'bdyT'

        varoxyID.units = 'mmolO2/m3'
        varoxyID.short_name = 'OXY'
        varoxyID.long_name = 'OXY'
        varoxyID.grid = 'bdyT'

        if grd == 'I':
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
    elif grd == 'U':
        varztID.axis = 'Depth'
        varztID.short_name = 'depthu'
        varztID.units = 'm'
        varztID.long_name = 'Depth'

        vardzID.axis = 'Depth'
        vardzID.short_name = 'e3u'
        vardzID.units = 'm'
        vardzID.long_name = 'Depth'

        varbtuID.units = 'm/s'
        varbtuID.short_name = 'vobtcrtx'
        varbtuID.long_name = 'Thickness-weighted depth-averaged zonal Current'
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

        vardzID.axis = 'Depth'
        vardzID.short_name = 'e3v'
        vardzID.units = 'm'
        vardzID.long_name = 'Depth'

        varbtvID.units = 'm/s'
        varbtvID.short_name = 'vobtcrty'
        varbtvID.long_name = 'Thickness-weighted depth-averaged meridional Current'
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

        varmskID.short_name = 'bdy_msk'
        varmskID.units = 'unitless'
        varmskID.long_name = 'Structured boundary mask'

    else:
        logging.error('Unknown Grid')

    ncid.close()

