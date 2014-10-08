###
## This is the equivalent to the matlab of same name. 
## The driver thast makes things happen. bit messy atm
## Note: Individual matlab vars are grouped into relevant dictionaries
## for clarity/ease of use
##
## Written by John Kazimierz Farey, Sep 2012
## Port of Matlab code of the Cardinal James Harle
###



from time import time, clock
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

from netCDF4 import Dataset
import numpy as np
from scipy.interpolate import interp1d
#import pickle

import logging
from netcdftime import date2num, datetime
logging.basicConfig(level=logging.INFO)
    
def go(): 
    #Logger
    logger = logging.getLogger(__name__)
    logger.info('START')   
    start = clock()
    SourceCoord = source_coord.SourceCoord()
    DstCoord = dst_coord.DstCoord()
    
    logger.info( clock() - start)
    start = clock()
    Setup = setup.Setup(0) # default settings file
    settings = Setup.settings

    
    logger.info( clock() - start)
    ice = settings['ice']
    logger.info('ice = %s', ice)

    logger.info( 'Done Setup')
    # default file, region settingas
    
    #This is to create mask. a future improvement will be to use GUI.
    start = clock()
    Mask = msk.Mask(0, med=1, blk=1, hud=1, bal=1, v2=1, custom_areas=None)

    logger.info( clock() - start)
    logger.info( 'Done Mask')
        
    bdy_msk = Mask.bdy_msk
    DstCoord.bdy_msk = bdy_msk == 1
    
    start = clock()
    logger.info( 'start bdy_t')
    Grid_T = gen_grid.Boundary(bdy_msk, settings, 't')
    
    logger.info( clock() - start)
    start = clock()
    logger.info( 'start bdy_u')
    Grid_U = gen_grid.Boundary(bdy_msk, settings, 'u')
    logger.info( 'start bdy_v')
    
    logger.info( clock() - start)
    start = clock()
    Grid_V = gen_grid.Boundary(bdy_msk, settings, 'v')
    logger.info('start bdy_f')
    
    logger.info( clock() - start)
    start = clock()
    
    Grid_F = gen_grid.Boundary(bdy_msk, settings, 'f')
    logger.info( 'done  bdy t,u,v,f')

    logger.info( clock() - start)
    start = clock()
    if ice:
        Grid_ice = nemo_bdy_ice.BoundaryIce()

    

   
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
        logger.info( 'bdy_ind %s %s %s', k, bdy_ind[k].bdy_i.shape, bdy_ind[k].bdy_r.shape)
    
    start = clock()
    co_set = coord.Coord(None, bdy_ind)
    logger.info( 'done coord gen')
    logger.info( clock() - start)
    start = clock()
    logger.info( settings['dst_hgr'])
    co_set.populate(settings['dst_hgr'])
    logger.info( 'done coord pop')
    
    logger.info( clock() - start)
    
    
    Grid_U.bdy_i = Grid_U.bdy_i[Grid_U.bdy_r == 1, :]
    Grid_V.bdy_i = Grid_V.bdy_i[Grid_V.bdy_r == 1, :]
    # Ought to look at this
    bdy_r = Grid_T.bdy_r

    # number of points
    num_bdy = {}
    num_bdy['t'] = len(Grid_T.bdy_i[:,0])
    num_bdy['z'] = len(Grid_T.bdy_i[bdy_r == 1, 0])
    num_bdy['u'] = len(Grid_U.bdy_i[:,0])
    num_bdy['v'] = len(Grid_V.bdy_i[:,0])

    ######
    # Gather grid info
    ######
    logger.info( 'gather grid info')
    start = clock()
    # how grid broken up, long/lats
    nc = Dataset(settings['src_zgr'], 'r')
    SourceCoord.zt = nc.variables['gdept_0'][:]
    nc.close()

    logger.info( clock() - start)
    logger.info( 'generating depth info')
    start = clock()
    z = zgrv.Depth(Grid_T.bdy_i,Grid_U.bdy_i,Grid_V.bdy_i, settings)

    logger.info( clock() - start)
    zp =  z.zpoints
    zk = zp.keys()
    logger.info( zk)
    for zk in zp:
        print zk, ' is ', zp[zk].shape, ' and a ', np.sum(zp[zk][:]), ' max: ', np.amax(zp[zk][:10,:2], axis=0)
        
    # FIX? Consider renaming to bdy_depths
    
    DstCoord.depths = {'t': {}, 'u': {}, 'v': {}}
    #temp = ['t', 'u', 'v']
    #temp2 = ['bdy_hbat', 'bdy_dz', 'bdy_z']

    for g in 't', 'u', 'v':
        DstCoord.depths[g]['bdy_hbat'] = np.nanmax(zp['w' + g], axis=0)
        DstCoord.depths[g]['bdy_dz'] = np.diff(zp['w' + g], axis=0)
        DstCoord.depths[g]['bdy_z'] = zp[g]

    # Might want to stake up some heretics now
    logger.info( 'horizontal grid info')
    start = clock()
    # Horizontal grid info
    nc = Dataset(settings['src_hgr'], 'r')
    SourceCoord.lon = nc.variables['glamt'][:,:]
    SourceCoord.lat = nc.variables['gphit'][:,:]
    nc.close()

    logger.info( clock() - start)
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
    logger.info( 'read and assign netcdf data, lon lat')
    # Read and assign netcdf data
    for g in 't', 'u', 'v':
        DstCoord.lonlat[g]['lon'] = nc.variables['glam' + g][0,:,:]
        DstCoord.lonlat[g]['lat'] = nc.variables['gphi' + g][0,:,:]
    nc.close()

    
    logger.info( clock() - start)
    # Set of 
    DstCoord.bdy_lonlat = {'t': {}, 'u': {}, 'v': {}, 'z': {}}
    for g in 't', 'u', 'v', 'z':
        for l in 'lon', 'lat':
            DstCoord.bdy_lonlat[g][l] = np.zeros(num_bdy[g])


    for g in 't', 'u', 'v':
        logger.info( 'range is %s', num_bdy[g])
        for i in range(num_bdy[g]):
            x = bdy_ind[g].bdy_i[i,1]
            y = bdy_ind[g].bdy_i[i,0]
            DstCoord.bdy_lonlat[g]['lon'][i] = DstCoord.lonlat[g]['lon'][x,y]
            DstCoord.bdy_lonlat[g]['lat'][i] = DstCoord.lonlat[g]['lat'][x,y]

        DstCoord.lonlat[g]['lon'][DstCoord.lonlat[g]['lon'] > 180] -= 360

    DstCoord.bdy_lonlat['z']['lon'] = DstCoord.bdy_lonlat['t']['lon'][bdy_r == 1]
    DstCoord.bdy_lonlat['z']['lat'] = DstCoord.bdy_lonlat['t']['lat'][bdy_r == 1]

    logger.info( DstCoord.bdy_lonlat['t']['lon'].shape)


    # TEMP


     #################
    ## Set Constants ##
     #################
    print 'src time'
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

    print clock() - start
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

            if Setup.settings['tra'] :
                extract = nemo_bdy_extr_tm3.Extract(Setup.settings,SourceCoord,DstCoord,Grid_T,['votemper','vosaline'],
                                        years=year,months=month) 
                extract.extract_month(year, month)
                #Get Date as a Number used in interpolation
                time_counter=np.zeros([len(extract.sc_time)])
                for i in range(0,len(extract.sc_time)):
                    time_counter[i] = extract.convert_date_to_destination_num(extract.sc_time[i]['date'])
                
                date_start = datetime(year, month, 1, 0, 0, 0)
                date_end = datetime(year, month, monthrange(year,month)[1],24, 60,60)
                time_num_start = extract.convert_date_to_destination_num(date_start)
                time_num_end = extract.convert_date_to_destination_num(date_end)
                    
                if monthrange(year,month)[1] == 29:
                    interp1fn = interp1d(np.arange(0,len(extract.sc_time),1),time_counter) 
                    time_counter = interp1fn(np.append(np.arange(0,len(extract.sc_time)-2+0.2,0.2),
                                                       np.arange(len(extract.sc_time)-2+1/6.0,len(extract.sc_time)-1+1/6.0,1/6.0)))      #added 0.2 to include boundary              
                else:
                    interp1fn = interp1d(np.arange(0,len(extract.sc_time),1),time_counter) 
                    time_counter = interp1fn(np.arange(0,len(extract.sc_time)-1+0.2,0.2))      #added 0.2 to include boundary   
                 
                ft = np.where(((time_counter >= time_num_start) - (time_counter <= time_num_end)).all())    
                time_counter = time_counter[ft]
                    
                for variable_name in extract.d_bdy:
                    if monthrange(year,month)[1] == 29:
                        x = np.squeeze(np.asarray(np.mat(np.arange(0,len(extract.d_bdy[variable_name][year]['data'][1,1,:]),1)).transpose()))                        
                        y = extract.d_bdy[variable_name][year]['data'].transpose(2, 0, 1)
                        interp1fn = interp1d(x,y,axis=0)
                        extra_axis = np.append(np.arange(0,len(extract.d_bdy[variable_name][year]['data'][1,1,:])-2+0.2,0.2),
                                                       np.arange(len(extract.d_bdy[variable_name][year]['data'][1,1,:])-2+1/6.0,
                                                                 len(extract.d_bdy[variable_name][year]['data'][1,1,:])-1+1/6.0,1/6.0))
                        extract.d_bdy[variable_name][year]['data'] = interp1fn(extra_axis)
                    else:
                        x = np.squeeze(np.asarray(np.mat(np.arange(0,len(extract.d_bdy[variable_name][year]['data'][1,1,:]),1)).transpose()))
                        y = extract.d_bdy[variable_name][year]['data'].transpose(2, 0, 1)
                        interp1fn = interp1d(x,y,axis=0)
                        extract.d_bdy[variable_name][year]['data'] = interp1fn(np.mat(np.arange(0,len(extract.d_bdy[variable_name][year]['data'][1,1,:])-1+0.2,0.2)).transpose())      #added 0.2 to include boundary
                    #extract.d_bdy[variable_name][year]['data'] = extract.d_bdy[variable_name][year]['data'][ft,:,:].transpose([1,2,0]) #not working
                        
                
            if Setup.settings['dyn2d'] or Setup.settings['dyn3d'] :
                extract_u = nemo_bdy_extr_tm3.Extract(Setup.settings,SourceCoord,DstCoord,Grid_U,['vozocrtx'], 
                                        years=year,months=month)
                extract_u.extract_month(year,month)
                extract_v = nemo_bdy_extr_tm3.Extract(Setup.settings,SourceCoord,DstCoord,Grid_V,['vomecrty'],
                                        years=year,months=month)
                extract_v.extract_month(year,month)
                
                
            
                  

go()
