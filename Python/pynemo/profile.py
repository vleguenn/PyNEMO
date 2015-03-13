###
## This is the equivalent to the matlab of same name.
## The driver thast makes things happen. bit messy atm
## Note: Individual matlab vars are grouped into relevant dictionaries
## for clarity/ease of use
##
## Algorithm:
## 1) Interpolate the data from the global grid to local grid
## 2) then interpolate the local grid data time data to input time resolution
## Written by John Kazimierz Farey, Sep 2012
## Port of Matlab code of the Cardinal James Harle
###


# pylint: disable=E1103
# pylint: disable=no-name-in-module
from time import clock
from calendar import monthrange

import nemo_bdy_setup as setup
import nemo_bdy_msk_c as msk
import nemo_bdy_gen_c as gen_grid
import nemo_coord_gen_pop as coord
import nemo_bdy_zgrv2 as zgrv
import nemo_bdy_src_time as source_time

import nemo_bdy_source_coord as source_coord
import nemo_bdy_dst_coord as dst_coord
import nemo_bdy_ice
import nemo_bdy_extr_tm3
from tide import nemo_bdy_tide3
from tide import nemo_bdy_tide_ncgen

from netCDF4 import Dataset
import numpy as np
from scipy.interpolate import interp1d

import nemo_bdy_ncgen
import nemo_bdy_ncpop
from gui.nemo_bdy_mask import Mask as Mask_File
#import pickle

import logging
from netcdftime import datetime
import copy

def go(setup_filepath=0, mask_gui=False):
    #Logger
    logger = logging.getLogger(__name__)
    logger.info('START')

    start = clock()
    SourceCoord = source_coord.SourceCoord()
    DstCoord = dst_coord.DstCoord()

    logger.info(clock() - start)
    start = clock()
    Setup = setup.Setup(setup_filepath) # default settings file
    settings = Setup.settings

    logger.info(clock() - start)
    ice = settings['ice']
    logger.info('ice = %s', ice)

    logger.info('Done Setup')
    # default file, region settingas

    # Use mask from the file if exists otherwise use default mask from nemo_bdy_msk_c
    start = clock()
    if mask_gui:
        #Open the gui to create a mask
        #open_settings_window()
        Setup.refresh()
    else:
        try:
            if Setup.bool_settings['mask_file'] and settings['mask_file'] is not None:        
                Mask =  Mask_File(settings['bathy'], settings['mask_file'])
                bdy_msk = Mask.data
            elif Setup.bool_settings['mask_file']:
                logger.error("Mask file is not given")
                return
            else:        
                logger.warning("Using default mask with bathymetry!!!!")
                Mask = msk.Mask(settings['bathy'], med=1, blk=1, hud=1, bal=1, v2=1, custom_areas=None)
                bdy_msk = Mask.bdy_msk
        except:
            return

    logger.info(clock() - start)
    logger.info('Done Mask')

    
    DstCoord.bdy_msk = bdy_msk == 1

    start = clock()
    logger.info('start bdy_t')
    Grid_T = gen_grid.Boundary(bdy_msk, settings, 't')

    logger.info(clock() - start)
    start = clock()
    logger.info('start bdy_u')
    Grid_U = gen_grid.Boundary(bdy_msk, settings, 'u')
    logger.info('start bdy_v')

    logger.info(clock() - start)
    start = clock()
    Grid_V = gen_grid.Boundary(bdy_msk, settings, 'v')
    logger.info('start bdy_f')

    logger.info(clock() - start)
    start = clock()

    Grid_F = gen_grid.Boundary(bdy_msk, settings, 'f')
    logger.info('done  bdy t,u,v,f')

    logger.info(clock() - start)
    start = clock()
    if ice:
        Grid_ice = nemo_bdy_ice.BoundaryIce()
        Grid_ice.grid_type = 't'
        Grid_ice.bdy_r = Grid_T.bdy_r

    # Be careful this stays as refs not copies
    '''
    bdy_ind = {}
    bdy_ind['t'] = [Grid_T.bdy_i, Grid_T.bdy_r]

    bdy_ind['u'] = [Grid_U.bdy_i, Grid_U.bdy_r]
    bdy_ind['v'] = [Grid_V.bdy_i, Grid_V.bdy_r] 
    bdy_ind['f'] = [Grid_F.bdy_i, Grid_F.bdy_r]

    '''
    bdy_ind = {'t': Grid_T, 'u': Grid_U, 'v': Grid_V, 'f': Grid_F}

    for k in bdy_ind.keys():
        logger.info('bdy_ind %s %s %s', k, bdy_ind[k].bdy_i.shape, bdy_ind[k].bdy_r.shape)

    start = clock()
    co_set = coord.Coord(settings['dst_dir']+'/toreador.nc', bdy_ind)
    logger.info('done coord gen')
    logger.info(clock() - start)
    start = clock()
    logger.info(settings['dst_hgr'])
    co_set.populate(settings['dst_hgr'])
    logger.info('done coord pop')

    logger.info(clock() - start)

    # may need to rethink grid info
    # tracer 3d frs over rw
    # tracer 2d frs over rw (i.e. ice)
    # dyn 2d over 1st rim of T grid (i.e. ssh)
    # dyn 2d over 1st rim
    # dyn 2d frs over rw
    # dyn 3d over 1st rim
    # dyn 3d frs over rw

    Grid_U.bdy_i = Grid_U.bdy_i[Grid_U.bdy_r == 0, :]
    Grid_V.bdy_i = Grid_V.bdy_i[Grid_V.bdy_r == 0, :]
    # Ought to look at this
    bdy_r = Grid_T.bdy_r

    # number of points
    num_bdy = {}
    num_bdy['t'] = len(Grid_T.bdy_i[:, 0])
    num_bdy['z'] = len(Grid_T.bdy_i[bdy_r == 0, 0])
    num_bdy['u'] = len(Grid_U.bdy_i[:, 0])
    num_bdy['v'] = len(Grid_V.bdy_i[:, 0])

    ######
    # Gather grid info
    ######
    logger.info('gather grid info')
    start = clock()
    # how grid broken up, long/lats
    nc = Dataset(settings['src_zgr'], 'r')
    SourceCoord.zt = nc.variables['gdept_0'][:]
    nc.close()

    logger.info(clock() - start)
    logger.info('generating depth info')
    start = clock()
    z = zgrv.Depth(Grid_T.bdy_i, Grid_U.bdy_i, Grid_V.bdy_i, settings)

    logger.info(clock() - start)
    zp = z.zpoints
    zk = zp.keys()
    logger.info(zk)
    for zk in zp:
        logger.info('%s  is %s and a %s  max: %s', zk, zp[zk].shape, np.sum(zp[zk][:]),
                    np.amax(zp[zk][:10, :2], axis=0))

    # FIX? Consider renaming to bdy_depths

    DstCoord.depths = {'t': {}, 'u': {}, 'v': {}}
    #temp = ['t', 'u', 'v']
    #temp2 = ['bdy_hbat', 'bdy_dz', 'bdy_z']

    for g in 't', 'u', 'v':
        DstCoord.depths[g]['bdy_hbat'] = np.nanmax(zp['w' + g], axis=0)
        DstCoord.depths[g]['bdy_dz'] = np.diff(zp['w' + g], axis=0)
        DstCoord.depths[g]['bdy_z'] = zp[g]

    # Might want to stake up some heretics now
    logger.info('horizontal grid info')
    start = clock()
    # Horizontal grid info
    nc = Dataset(settings['src_hgr'], 'r')
    SourceCoord.lon = nc.variables['glamt'][:, :]
    SourceCoord.lat = nc.variables['gphit'][:, :]
    nc.close()

    logger.info(clock() - start)
    DstCoord.lonlat = {'t': {}, 'u': {}, 'v': {}}

    start = clock()
    nc = Dataset(settings['dst_hgr'], 'r')
    """
    DstCoord.lonlat['t']['lon'] = nc.variables['glamt'][0,:,:]
    DstCoord.lonlat['t']['lat'] = nc.variables['gphit'][0,:,:]
    DstCoord.lonlat['u']['lon'] = nc.variables['glamu'][0,:,:]
    DstCoord.lonlat['u']['lat'] = nc.variables['gphiu'][0,:,:]
    DstCoord.lonlat['v']['lon'] = nc.variables['glamv'][0,:,:]
    DstCoord.lonlat['v']['lat'] = nc.variables['gphiv'][0,:,:]
    """
    logger.info('read and assign netcdf data, lon lat')
    # Read and assign netcdf data
    for g in 't', 'u', 'v':
        DstCoord.lonlat[g]['lon'] = nc.variables['glam' + g][0, :, :]
        DstCoord.lonlat[g]['lat'] = nc.variables['gphi' + g][0, :, :]
    nc.close()

    logger.info(clock() - start)
    # Set of
    DstCoord.bdy_lonlat = {'t': {}, 'u': {}, 'v': {}, 'z': {}}
    for g in 't', 'u', 'v', 'z':
        for l in 'lon', 'lat':
            DstCoord.bdy_lonlat[g][l] = np.zeros(num_bdy[g])


    for g in 't', 'u', 'v':
        logger.info('range is %s', num_bdy[g])
        for i in range(num_bdy[g]):
            x = bdy_ind[g].bdy_i[i, 1]
            y = bdy_ind[g].bdy_i[i, 0]
            DstCoord.bdy_lonlat[g]['lon'][i] = DstCoord.lonlat[g]['lon'][x, y]
            DstCoord.bdy_lonlat[g]['lat'][i] = DstCoord.lonlat[g]['lat'][x, y]

        DstCoord.lonlat[g]['lon'][DstCoord.lonlat[g]['lon'] > 180] -= 360

    DstCoord.bdy_lonlat['z']['lon'] = DstCoord.bdy_lonlat['t']['lon'][bdy_r == 0]
    DstCoord.bdy_lonlat['z']['lat'] = DstCoord.bdy_lonlat['t']['lat'][bdy_r == 0]

    logger.info(DstCoord.bdy_lonlat['t']['lon'].shape)

    # TEMP
    #################
    ## Set Constants ##
    #################
    logger.debug('src time')
    start = clock()
    # Create a Time handler object
    SourceTime = source_time.SourceTime(settings['src_dir'])
    acc = settings['src_time_adj']
#    acc = -3168000 - 86400 # Straight out Matlab code. Question this.
#    acc = 0
    # Extract source data on dst grid
    # Assign time info to Boundary T, U, V objects
    # FIX MY STRUCTURING
    Grid_T.source_time = SourceTime.get_source_time('t', acc)
    Grid_U.source_time = SourceTime.get_source_time('u', acc)
    Grid_V.source_time = SourceTime.get_source_time('v', acc)
    if ice:
        Grid_ice.source_time = SourceTime.get_source_time('i', acc)

    logger.debug(clock() - start)
    """
    # List source filenames

    for grid in [Grid_T, Grid_U, Grid_V]:
        grid.dir_list = grid.get_dir_list(Setup.settings['src_dir'])
    # This might not be best way
    if Setup.settings['ice']:
        Grid_ice.dir_list = Grid_T.get_dir_list(Setup.settings['src_dir'], ice=True)
    """

    # Set constants before years loop

    # Watch out for the eternal Fortran C war
    # subdomain size
    j, i = DstCoord.lonlat['t']['lon'].shape
    k = DstCoord.depths['t']['bdy_z'].shape[0]

    unit_string = 'seconds since %d-01-01 00:00:00' %settings['base_year']
    unit_origin = '%d-01-01 00:00:00' %settings['base_year']
    """
    ################
    # Extract source data on dst grid
    ################

    # extract source date/time info from source filenames prior to entering main loop

    source_time = {}
    source_time['t'] = nemo_bdy_src_time
    """



    # Example tide extraction TODO complete loop over constituents read in from namelist
#    for Grid in [Grid_T,Grid_U,Grid_V]:
        #extract_write_tidal_data(Setup, DstCoord, Grid, num_bdy, settings["clname"])
    cosz, sinz, cosu, sinu, cosv, sinv = nemo_bdy_tide3.nemo_bdy_tpx7p2_rot(Setup, DstCoord,
                                                                            Grid_T, Grid_U, Grid_V,
                                                                            settings["clname"])
    write_tidal_data(Setup, DstCoord, Grid_T, num_bdy, settings["clname"], cosz, sinz)
    write_tidal_data(Setup, DstCoord, Grid_U, num_bdy, settings["clname"], cosu, sinu)
    write_tidal_data(Setup, DstCoord, Grid_V, num_bdy, settings["clname"], cosv, sinv)


#        extract_tide_z = nemo_bdy_tide2.Extract(Setup.settings,DstCoord,Grid_T)
#        extract_tide_z.extract_con('m2')
#        extract_tide_u = nemo_bdy_tide2.Extract(Setup.settings,DstCoord,Grid_U)
#        extract_tide_u.extract_con('m2')
#        extract_tide_v = nemo_bdy_tide2.Extract(Setup.settings,DstCoord,Grid_V)
#        extract_tide_v.extract_con('m2')



    # Enter Years Loop
    # Loop not yet implemented
    years = [1979]
    months = [11]



    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    #print '\n- Old Callington Ratsmith should never be trusted -\n'
    #print 'DstCoord: \n', dir(DstCoord), '\nSourceCoord: \n', dir(SourceCoord)

    #ret1 = {'t': Grid_T, 'u':Grid_U, 'v':Grid_V, 'f':Grid_F}
    # return settings, ret1, SourceCoord, DstCoord

    for year in years:
        for month in months:

            if Setup.settings['tra']:
                extract_t = nemo_bdy_extr_tm3.Extract(Setup.settings, SourceCoord, DstCoord,
                                                      Grid_T, ['votemper','vosaline'],
                                                      years=year, months=month) 
                extract_t.extract_month(year, month)
                #Get Date as a Number used in interpolation
                time_counter = np.zeros([len(extract_t.sc_time)])
                for i in range(0,len(extract_t.sc_time)):
                    time_counter[i] = extract_t.convert_date_to_destination_num(extract_t.sc_time[i]['date'])

                date_start = datetime(year, month, 1, 0, 0, 0)
                date_end = datetime(year, month, monthrange(year, month)[1], 24, 60, 60)
                time_num_start = extract_t.convert_date_to_destination_num(date_start)
                time_num_end = extract_t.convert_date_to_destination_num(date_end)

                if monthrange(year,month)[1] == 29:
                    interp1fn = interp1d(np.arange(0,len(extract_t.sc_time),1),time_counter)
                    time_counter = interp1fn(np.append(np.arange(0,
                                                                 len(extract_t.sc_time)-2+0.2,
                                                                 0.2),
                                                       np.arange(len(extract_t.sc_time)-2+1/6.0,
                                                                 len(extract_t.sc_time)-1+1/6.0,
                                                                 1/6.0)))      #added 0.2 to include boundary              
                else:
                    interp1fn = interp1d(np.arange(0, len(extract_t.sc_time), 1), time_counter)
                    time_counter = interp1fn(np.arange(0, len(extract_t.sc_time)-1+0.2, 0.2))      #added 0.2 to include boundary   

                ft = np.where(((time_counter >= time_num_start) & (time_counter <= time_num_end)))    
                time_counter = time_counter[ft]

                interpolate_data(extract_t, year, month, ft) # interpolate the data to daily period in a month                                                              
                        
                output_filename_t = Setup.settings['dst_dir']+Setup.settings['fn']+ \
                                    '_bdyT_y'+str(year)+'m'+str(month)+'.nc'
                logger.info('Outputing T file:%s', output_filename_t)
                grid_id = 'T'
                if Setup.settings['ice']:
                    grid_id = 'I'
                nemo_bdy_ncgen.CreateBDYNetcdfFile(output_filename_t, num_bdy['t'], 
                                                   DstCoord.lonlat['t']['lon'].shape[1],
                                                   DstCoord.lonlat['t']['lon'].shape[0],
                                                   DstCoord.depths['t']['bdy_z'].shape[0],
                                                   Setup.settings['rimwidth'],
                                                   Setup.settings['dst_metainfo'],
                                                   unit_origin, Setup.settings['fv'],
                                                   Setup.settings['dst_calendar'],
                                                   grid_id)
                #Replace nan values with fillvalues
                nanindex = np.isnan(extract_t.d_bdy['votemper'][year]['data'])
                extract_t.d_bdy['votemper'][year]['data'][nanindex] = Setup.settings['fv']    
                nanindex = np.isnan(extract_t.d_bdy['vosaline'][year]['data'])      
                extract_t.d_bdy['vosaline'][year]['data'][nanindex] = Setup.settings['fv']
                
                nemo_bdy_ncpop.WriteDataToFile(output_filename_t, 'votemper',
                                               extract_t.d_bdy['votemper'][year]['data'])
                nemo_bdy_ncpop.WriteDataToFile(output_filename_t, 'vosaline',
                                               extract_t.d_bdy['vosaline'][year]['data'])
                    
                nemo_bdy_ncpop.WriteDataToFile(output_filename_t, 'bdy_msk', DstCoord.bdy_msk) 
                nemo_bdy_ncpop.WriteDataToFile(output_filename_t, 'deptht',
                                               DstCoord.depths['t']['bdy_z']) 
                nemo_bdy_ncpop.WriteDataToFile(output_filename_t, 'nav_lon',
                                               DstCoord.lonlat['t']['lon']) 
                nemo_bdy_ncpop.WriteDataToFile(output_filename_t, 'nav_lat',
                                               DstCoord.lonlat['t']['lat']) 
                nemo_bdy_ncpop.WriteDataToFile(output_filename_t, 'nbidta', Grid_T.bdy_i[:,0]+1)          #Added to match the index of matlab                       
                nemo_bdy_ncpop.WriteDataToFile(output_filename_t, 'nbjdta', Grid_T.bdy_i[:,1]+1)          #Added to match the index of matlab
                nemo_bdy_ncpop.WriteDataToFile(output_filename_t, 'nbrdta', Grid_T.bdy_r[:]+1)            #Added to match the index of matlab
                nemo_bdy_ncpop.WriteDataToFile(output_filename_t, 'time_counter',time_counter)
 
                if Setup.settings['ice']:
                    extract_write_ice_data(Setup, SourceCoord, DstCoord, Grid_ice, year, month,
                                           ft, num_bdy, time_counter, unit_origin, logger)               
                #-------------------------bt---------------------------------------------------------------------------
                extract_write_bt_data(Setup, SourceCoord, DstCoord, Grid_T, year, month,
                                      ft, num_bdy, time_counter, unit_origin, logger)  
            #Eco

            #U, V
            if Setup.settings['dyn2d'] or Setup.settings['dyn3d'] :
                #---------------------------------------------U----------------------------------------------------------                
                class Grid_2(object):pass
                Grid_U2 = Grid_2()
                Grid_U2.grid_type = 'u'    
                Grid_U2.bdy_i = Grid_U.bdy_i
                Grid_U2.maxI = DstCoord.lonlat['t']['lon'].shape[1]
                Grid_U2.maxJ = DstCoord.lonlat['t']['lon'].shape[0]    
                Grid_U2.fname_2 = Grid_V.source_time
                extract_write_u_data(Setup, SourceCoord, DstCoord, Grid_U, Grid_U2, year, month,
                                     ft, num_bdy, time_counter, unit_origin, logger)                
                #----------------------------------------------V--------------------------------               
                Grid_V2 = Grid_2()
                Grid_V2.grid_type = 'v'
                Grid_V2.bdy_i = Grid_V.bdy_i
                Grid_V2.maxI = DstCoord.lonlat['t']['lon'].shape[1]
                Grid_V2.maxJ = DstCoord.lonlat['t']['lon'].shape[0]    
                Grid_V2.fname_2 = Grid_V.source_time   
                Grid_UV = Grid_2()
                Grid_UV.bdy_i = Grid_V.bdy_i
                Grid_UV.grid_type = 'v'
                Grid_UV.bdy_r = Grid_V.bdy_r
                Grid_UV.source_time = Grid_U.source_time
                extract_write_v_data(Setup, SourceCoord, DstCoord, Grid_UV, Grid_V2, year, month,
                                     ft, num_bdy, time_counter, unit_origin, logger)                 
                
def extract_write_ice_data(Setup, SourceCoord, DstCoord, grid_ice, year, month, ft, num_bdy,
                           time_counter, unit_origin, logger):
    """ writes the ice data to file"""
    source_coord_var_copy = copy.deepcopy(SourceCoord)
    source_coord_var_copy.zt = np.zeros([1, 1])
    dst_coord_var_copy = copy.deepcopy(DstCoord)
    dst_coord_var_copy.depths['t']['bdy_z'] = np.zeros([1, 1])   
    extract_ice = nemo_bdy_extr_tm3.Extract(Setup.settings, source_coord_var_copy, 
                                            dst_coord_var_copy, grid_ice,
                                            ['iicethic', 'ileadfra', 'isnowthi'],
                                            years=year, months=month)
    extract_ice.extract_month(year,month)
    interpolate_data(extract_ice, year, month,ft) #interpolate the data to daily period in a month
                        
    output_filename_t = Setup.settings['dst_dir']+Setup.settings['fn']+'_bdyT_y'+str(year)+ \
                        'm'+str(month)+'.nc'
    logger.info('Outputing ICE values to T file:%s', output_filename_t)                        

    nemo_bdy_ncpop.WriteDataToFile(output_filename_t, 'iicethic',
                                   extract_ice.d_bdy['iicethic'][year]['data'])
    nemo_bdy_ncpop.WriteDataToFile(output_filename_t, 'ileadfra',
                                   extract_ice.d_bdy['ileadfra'][year]['data'])
    nemo_bdy_ncpop.WriteDataToFile(output_filename_t, 'isnowthi',
                                   extract_ice.d_bdy['isnowthi'][year]['data'])
                               
def extract_write_bt_data(Setup, SourceCoord, DstCoord, Grid_T, year, month, 
                          ft, num_bdy, time_counter, unit_origin, logger):
    """Writes BT data file"""
    SourceCoordLatLon = copy.deepcopy(SourceCoord)
    SourceCoordLatLon.zt = np.zeros([1,1])
    DstCoordLatLon = copy.deepcopy(DstCoord)
    DstCoordLatLon.depths['t']['bdy_z'] = np.zeros([1,1])
    Grid_T.grid_type = 'z'
    extract_t = nemo_bdy_extr_tm3.Extract(Setup.settings,
                                          SourceCoordLatLon,
                                          DstCoord,
                                          Grid_T,
                                          ['sossheig'],
                                          years=year,
                                          months=month)
    extract_t.extract_month(year,month)
    interpolate_data(extract_t,year,month,ft) #interpolate the data to daily period in a month
    output_filename_bt = Setup.settings['dst_dir']+Setup.settings['fn']+'_bt_bdyT_y'+str(year)+'m'+str(month)+'.nc'
    logger.info('Outputting bt_T file: %s', output_filename_bt)
    nemo_bdy_ncgen.CreateBDYNetcdfFile(output_filename_bt, num_bdy['z'], 
                                       DstCoord.lonlat['t']['lon'].shape[1],
                                       DstCoord.lonlat['t']['lon'].shape[0], 
                                       1, Setup.settings['rimwidth'],
                                       Setup.settings['dst_metainfo'],
                                       unit_origin, Setup.settings['fv'],
                                       Setup.settings['dst_calendar'],'Z')             
    nemo_bdy_ncpop.WriteDataToFile(output_filename_bt, 'sossheig', extract_t.d_bdy['sossheig'][year]['data'])
    nemo_bdy_ncpop.WriteDataToFile(output_filename_bt, 'nav_lon', DstCoord.lonlat['t']['lon']) 
    nemo_bdy_ncpop.WriteDataToFile(output_filename_bt, 'nav_lat', DstCoord.lonlat['t']['lat']) 
    nemo_bdy_ncpop.WriteDataToFile(output_filename_bt, 'nbidta',
                                   Grid_T.bdy_i[np.where(Grid_T.bdy_r==0),0]+1)                                 
    nemo_bdy_ncpop.WriteDataToFile(output_filename_bt, 'nbjdta',
                                   Grid_T.bdy_i[np.where(Grid_T.bdy_r==0),1]+1)
    nemo_bdy_ncpop.WriteDataToFile(output_filename_bt, 'nbrdta',
                                   Grid_T.bdy_r[np.where(Grid_T.bdy_r==0)]+1)
    nemo_bdy_ncpop.WriteDataToFile(output_filename_bt, 'time_counter',time_counter) 
                
def extract_write_u_data(setup_var,source_coord_var,dst_coord_var,grid_u,grid_u2,year,month,ft,num_bdy,time_counter,unit_origin,logger):
    """ Extracts the U component of velocity data, interpolates the data for a month and writes  data to netCDF file"""
    extract_u = nemo_bdy_extr_tm3.Extract(setup_var.settings, source_coord_var, dst_coord_var,
                                          grid_u, ['vozocrtx','vomecrty'],Grid_2=grid_u2,
                                          years=year, months=month)
    extract_u.extract_month(year, month)
    interpolate_data(extract_u, year, month, ft)
    nanindex = np.isnan(extract_u.d_bdy['vozocrtx'][year]['data'])
    tmp_vozocrtx = extract_u.d_bdy['vozocrtx'][year]['data']
    tmp_vozocrtx[nanindex] = setup_var.settings['fv']
    output_filename_u = setup_var.settings['dst_dir'] + setup_var.settings['fn'] + '_bdyU_y' + \
                        str(year) + 'm' + str(month) + '.nc'
    logger.info('Outputting U file: %s', output_filename_u)
    nemo_bdy_ncgen.CreateBDYNetcdfFile(output_filename_u, num_bdy['u'],
                                       dst_coord_var.lonlat['u']['lon'].shape[1],
                                       dst_coord_var.lonlat['u']['lon'].shape[0],
                                       dst_coord_var.depths['u']['bdy_z'].shape[0],
                                       setup_var.settings['rimwidth'],
                                       setup_var.settings['dst_metainfo'],
                                       unit_origin,
                                       setup_var.settings['fv'],
                                       setup_var.settings['dst_calendar'],
                                       'U')
    nemo_bdy_ncpop.WriteDataToFile(output_filename_u, 'vozocrtx', tmp_vozocrtx)
    nemo_bdy_ncpop.WriteDataToFile(output_filename_u, 'depthu', dst_coord_var.depths['u']['bdy_z'])
    # TODO: write vobtcrtx
    tmp_vozocrtx[nanindex] = float('NaN')
    # ft size has some problem
    tmp_tile = np.tile(dst_coord_var.depths['u']['bdy_dz'], [ft[0].size, 1, 1, 1])
    tmp_vobtcrtx = np.nansum(np.reshape(tmp_vozocrtx[:, :-1, :], tmp_tile.shape) * tmp_tile, 2) / np.nansum(tmp_tile, 2)  
    tmp_vobtcrtx[np.isnan(tmp_vobtcrtx)] = setup_var.settings['fv']                                                                      
    nemo_bdy_ncpop.WriteDataToFile(output_filename_u, 'vobtcrtx', tmp_vobtcrtx)               
    nemo_bdy_ncpop.WriteDataToFile(output_filename_u, 'nav_lon', dst_coord_var.lonlat['u']['lon'])
    nemo_bdy_ncpop.WriteDataToFile(output_filename_u, 'nav_lat', dst_coord_var.lonlat['u']['lat'])
    nemo_bdy_ncpop.WriteDataToFile(output_filename_u, 'nbidta', grid_u.bdy_i[:, 0]+1)                
    nemo_bdy_ncpop.WriteDataToFile(output_filename_u, 'nbjdta', grid_u.bdy_i[:, 1]+1)                
    nemo_bdy_ncpop.WriteDataToFile(output_filename_u, 'nbrdta', np.ones([1, num_bdy['u']]))
    nemo_bdy_ncpop.WriteDataToFile(output_filename_u, 'time_counter', time_counter)                
    #-------------------------------------------End U-------------------------------    
    
def extract_write_v_data(setup_var, source_coord_var, dst_coord_var, grid_v, grid_v2, year, month,
                         ft, num_bdy, time_counter, unit_origin, logger):
    """ Extracts the V component of velocity data, interpolates the data for a month and writes  data to netCDF file"""
    extract_v = nemo_bdy_extr_tm3.Extract(setup_var.settings, source_coord_var, dst_coord_var,
                                          grid_v, ['vozocrtx','vomecrty'], Grid_2=grid_v2,
                                          years=year, months=month)
    extract_v.extract_month(year, month) #return values in vozocrtx instead of vomecrty
    interpolate_data(extract_v, year, month, ft)  
    nanindex = np.isnan(extract_v.d_bdy['vozocrtx'][year]['data'])
    tmp_vozocrty = extract_v.d_bdy['vozocrtx'][year]['data']
    tmp_vozocrty[nanindex] = setup_var.settings['fv']
                    
    output_filename_v = setup_var.settings['dst_dir']+setup_var.settings['fn']+'_bdyV_y'+str(year)+ \
                        'm'+str(month)+'.nc'
    logger.info('Outputting V file: %s', output_filename_v)
    nemo_bdy_ncgen.CreateBDYNetcdfFile(output_filename_v, num_bdy['v'], 
                                       dst_coord_var.lonlat['v']['lon'].shape[1],
                                       dst_coord_var.lonlat['v']['lon'].shape[0],
                                       dst_coord_var.depths['v']['bdy_z'].shape[0],
                                       setup_var.settings['rimwidth'],
                                       setup_var.settings['dst_metainfo'],
                                       unit_origin,
                                       setup_var.settings['fv'],
                                       setup_var.settings['dst_calendar'],
                                       'V')
    nemo_bdy_ncpop.WriteDataToFile(output_filename_v, 'vomecrty', 
                                   extract_v.d_bdy['vozocrtx'][year]['data'])
    nemo_bdy_ncpop.WriteDataToFile(output_filename_v, 'depthv', dst_coord_var.depths['v']['bdy_z']) 
    tmp_vozocrty[nanindex] = float('NaN')
    #ft size has some problem
    tmp_tile = np.tile(dst_coord_var.depths['v']['bdy_dz'],[ft[0].size,1,1,1])
    tmp_vobtcrty = np.nansum(np.reshape(tmp_vozocrty[:,:-1,:], tmp_tile.shape)*tmp_tile,2) \
                   /np.nansum(tmp_tile,2)  
    tmp_vobtcrty[np.isnan(tmp_vobtcrty)] = setup_var.settings['fv']
    nemo_bdy_ncpop.WriteDataToFile(output_filename_v, 'vobtcrty', tmp_vobtcrty) 
                                                  
    nemo_bdy_ncpop.WriteDataToFile(output_filename_v, 'nav_lon', dst_coord_var.lonlat['u']['lon'])
    nemo_bdy_ncpop.WriteDataToFile(output_filename_v, 'nav_lat', dst_coord_var.lonlat['u']['lat'])
    nemo_bdy_ncpop.WriteDataToFile(output_filename_v, 'nbidta', grid_v.bdy_i[:,0]+1)                
    nemo_bdy_ncpop.WriteDataToFile(output_filename_v, 'nbjdta', grid_v.bdy_i[:,1]+1)                
    nemo_bdy_ncpop.WriteDataToFile(output_filename_v, 'nbrdta', np.ones([1,num_bdy['v']]))
    nemo_bdy_ncpop.WriteDataToFile(output_filename_v, 'time_counter', time_counter)   
                
                
def write_tidal_data(setup_var, dst_coord_var, grid, num_bdy, tide_cons, cos_g, sin_g):
    """ This method writes the tidal data to netcdf file"""
    indx =0
    for tide_con in tide_cons:
        const_name = setup_var.settings['clname'][tide_con]
        const_name = const_name.replace("'","").upper() 

        if grid.grid_type == 't' :
            output_filename_tide = setup_var.settings['dst_dir']+setup_var.settings['fn']+"_bdytide_rotT_"+const_name+"_grid_"+grid.grid_type.upper()+".nc"
            nemo_bdy_tide_ncgen.CreateBDYTideNetcdfFile(output_filename_tide, num_bdy['z'], dst_coord_var.lonlat['t']['lon'].shape[1], dst_coord_var.lonlat['t']['lon'].shape[0], 
                                                        "tidal elevation components for:"+tide_con, setup_var.settings['fv'],'T')
            nemo_bdy_ncpop.WriteDataToFile(output_filename_tide, 'z1', cos_g[indx])
            nemo_bdy_ncpop.WriteDataToFile(output_filename_tide, 'z2', sin_g[indx])  
            nemo_bdy_ncpop.WriteDataToFile(output_filename_tide, 'bdy_msk', dst_coord_var.bdy_msk) 
            nemo_bdy_ncpop.WriteDataToFile(output_filename_tide, 'nav_lon', dst_coord_var.lonlat['t']['lon']) 
            nemo_bdy_ncpop.WriteDataToFile(output_filename_tide, 'nav_lat', dst_coord_var.lonlat['t']['lat']) 
            nemo_bdy_ncpop.WriteDataToFile(output_filename_tide, 'nbidta', grid.bdy_i[np.where(grid.bdy_r==0), 0]+1)                                 
            nemo_bdy_ncpop.WriteDataToFile(output_filename_tide, 'nbjdta', grid.bdy_i[np.where(grid.bdy_r==0), 1]+1)
            nemo_bdy_ncpop.WriteDataToFile(output_filename_tide, 'nbrdta', grid.bdy_r[np.where(grid.bdy_r==0)]+1)
        elif grid.grid_type == 'u':
            output_filename_tide = setup_var.settings['dst_dir']+setup_var.settings['fn']+"_bdytide_rotT_"+const_name+"_grid_"+grid.grid_type.upper()+".nc"
            nemo_bdy_tide_ncgen.CreateBDYTideNetcdfFile(output_filename_tide, num_bdy['u'], dst_coord_var.lonlat['t']['lon'].shape[1], dst_coord_var.lonlat['t']['lon'].shape[0], 
                                                        "tidal east velocity components for:"+tide_con, setup_var.settings['fv'],'U')
            nemo_bdy_ncpop.WriteDataToFile(output_filename_tide, 'u1', cos_g[indx])
            nemo_bdy_ncpop.WriteDataToFile(output_filename_tide, 'u2', sin_g[indx])              
            nemo_bdy_ncpop.WriteDataToFile(output_filename_tide, 'nav_lon', dst_coord_var.lonlat['u']['lon'])
            nemo_bdy_ncpop.WriteDataToFile(output_filename_tide, 'nav_lat', dst_coord_var.lonlat['u']['lat'])
            nemo_bdy_ncpop.WriteDataToFile(output_filename_tide, 'nbidta', grid.bdy_i[:, 0]+1)                
            nemo_bdy_ncpop.WriteDataToFile(output_filename_tide, 'nbjdta', grid.bdy_i[:, 1]+1)                
            nemo_bdy_ncpop.WriteDataToFile(output_filename_tide, 'nbrdta', np.ones([1, num_bdy['u']]))      
        elif grid.grid_type == 'v':
            output_filename_tide = setup_var.settings['dst_dir']+setup_var.settings['fn']+"_bdytide_rotT_"+const_name+"_grid_"+grid.grid_type.upper()+".nc"            
            nemo_bdy_tide_ncgen.CreateBDYTideNetcdfFile(output_filename_tide, num_bdy['v'], dst_coord_var.lonlat['t']['lon'].shape[1], dst_coord_var.lonlat['t']['lon'].shape[0], 
                                                        "tidal north velocity components for:"+tide_con, setup_var.settings['fv'],'V')
            nemo_bdy_ncpop.WriteDataToFile(output_filename_tide, 'v1', cos_g[indx])
            nemo_bdy_ncpop.WriteDataToFile(output_filename_tide, 'v2', sin_g[indx])                          
            nemo_bdy_ncpop.WriteDataToFile(output_filename_tide,'nav_lon',dst_coord_var.lonlat['v']['lon'])
            nemo_bdy_ncpop.WriteDataToFile(output_filename_tide,'nav_lat',dst_coord_var.lonlat['v']['lat'])
            nemo_bdy_ncpop.WriteDataToFile(output_filename_tide,'nbidta',grid.bdy_i[:,0]+1)                
            nemo_bdy_ncpop.WriteDataToFile(output_filename_tide,'nbjdta',grid.bdy_i[:,1]+1)                
            nemo_bdy_ncpop.WriteDataToFile(output_filename_tide,'nbrdta',np.ones([1,num_bdy['v']]))   
        indx=indx+1
                                                                                  
                     
        
def interpolate_data(extract, year, month, time_indexes):
    """This method does a 1D interpolation of daily mean data of a month"""
    
    for variable_name in extract.d_bdy:
        if monthrange(year,month)[1] == 29: #leap year
            lon = np.squeeze(np.asarray(np.mat(np.arange(0,len(extract.d_bdy[variable_name][year]['data'][:,0,0]),1)).transpose()))                        
            lat = extract.d_bdy[variable_name][year]['data'].transpose(0, 2, 1)
            interp1fn = interp1d(lon,lat,axis=0)
            extra_axis = np.append(np.arange(0,len(extract.d_bdy[variable_name][year]['data'][:,0,0])-2+0.2,0.2),
                                           np.arange(len(extract.d_bdy[variable_name][year]['data'][:,0,0])-2+1/6.0,
                                                     len(extract.d_bdy[variable_name][year]['data'][:,0,0])-1+1/6.0,1/6.0))
            extract.d_bdy[variable_name][year]['data'] = interp1fn(extra_axis)
        else:
            lon = np.squeeze(np.asarray(np.mat(np.arange(0,len(extract.d_bdy[variable_name][year]['data'][:,0,0]),1)).transpose()))
            lat = extract.d_bdy[variable_name][year]['data'].transpose(0, 2, 1)
            interp1fn = interp1d(lon,lat,axis=0)
            extract.d_bdy[variable_name][year]['data'] = interp1fn(np.squeeze(np.asarray(np.mat(np.arange(0,len(extract.d_bdy[variable_name][year]['data'][:,0,0])-1+0.2,0.2)).transpose())))      #added 0.2 to include boundary
        tmp = np.squeeze(extract.d_bdy[variable_name][year]['data'][time_indexes,:,:]) 
        if len(tmp.shape) == 2:            
            tmp = tmp.reshape(tmp.shape+(1L,))        
        extract.d_bdy[variable_name][year]['data'] =  tmp.transpose(0,2,1)
        
                  
