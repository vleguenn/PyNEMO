###
## This is the equivalent to the matlab of same name. 
## The driver thast makes things happen. bit messy atm
## Note: Individual matlab vars are grouped into relevant dictionaries
## for clarity/ease of use
##
## Written by John Kazimierz Farey, Sep 2012
## Port of Matlab code of the Cardinal James Harle
###


from timeit import Timer

import nemo_bdy_setup as setup
import nemo_bdy_msk_c as msk
import nemo_bdy_gen_c as gen_grid
import nemo_coord_gen_pop as coord
import nemo_bdy_zgrv2 as zgrv
import nemo_bdy_src_time as source_time

import nemo_bdy_source_coord as source_coord
import nemo_bdy_dst_coord as dst_coord
import nemo_bdy_ice

from netCDF4 import Dataset
import numpy as np
#import pickle

def go(): 
    print 'START'   
    
    SourceCoord = source_coord.SourceCoord()
    DstCoord = dst_coord.DstCoord()
    
    Setup = setup.Setup(0) # default settings file

    ice = Setup.settings['ice']
    print 'ice = ', ice

    print 'Done Setup'
    # default file, region settingas
    
    Mask = msk.Mask(0, med=1, blk=1, hud=1, bal=1, v2=1, reduced=0)

    print 'Done Mask'
    bdy_msk = Mask.bdy_msk
    DstCoord.bdy_msk = bdy_msk == 1
    
    print 'start bdy_t'
    Grid_T = gen_grid.Boundary(bdy_msk, Setup, 't')
    
    print 'start bdy_u'
    Grid_U = gen_grid.Boundary(bdy_msk, Setup, 'u')
    print 'start bdy_v'
    
    Grid_V = gen_grid.Boundary(bdy_msk, Setup, 'v')
    print 'start bdy_f'
    
    
    Grid_F = gen_grid.Boundary(bdy_msk, Setup, 'f')
    print 'done  bdy t,u,v,f'

    if ice:
        Grid_ice = nemo_bdy_ice.BoundaryIce()

    

   
    # Be careful this stays as refs not copies
    bdy_ind = {}
    bdy_ind['t'] = [Grid_T.bdy_i, Grid_T.bdy_r] 
    
    bdy_ind['u'] = [Grid_U.bdy_i, Grid_U.bdy_r] 
    bdy_ind['v'] = [Grid_V.bdy_i, Grid_V.bdy_r] 
    bdy_ind['f'] = [Grid_F.bdy_i, Grid_F.bdy_r]

    
    for k in bdy_ind.keys():
        print 'bdy_ind ', k, bdy_ind[k][0].shape, bdy_ind[k][1].shape
    
    return Setup, bdy_ind
    co_set = coord.Coord(None, bdy_ind)
    print 'done coord gen'
    co_set.populate(Setup.settings['dst_hgr'])
    print 'done coord pop'
    
    
    
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
    
    # how grid broken up, long/lats
    nc = Dataset(Setup.settings['src_zgr'], 'r')
    SourceCoord.zt = nc.variables['gdept_0'][:]
    nc.close()


    z = zgrv.Depth(Grid_T.bdy_i,Grid_U.bdy_i,Grid_V.bdy_i,Setup)

    zp =  z.zpoints
    zk = zp.keys()
    print zk
    for zk in zp:
        print zk, ' is ', zp[zk].shape, ' and a ', np.sum(zp[zk][:]), ' max: ', np.amax(zp[zk][:10,:2], axis=0)
        
    # FIX? Consider renaming to bdy_depths
    
    DstCoord.depths = {'t': {}, 'u': {}, 'v': {}}
    #temp = ['t', 'u', 'v']
    #temp2 = ['bdy_hbat', 'bdy_dz', 'bdy_z']

    for g in 't', 'u', 'v':
        DstCoord.depths[g]['bdy_hbat'] = np.amax(zp['w' + g], axis=0)
        DstCoord.depths[g]['bdy_dz'] = np.diff(zp['w' + g], axis=0)
        DstCoord.depths[g]['bdy_z'] = zp[g]

    # Might want to stake up some heretics now

    # Horizontal grid info
    nc = Dataset(Setup.settings['src_hgr'], 'r')
    SourceCoord.lon = nc.variables['glamt'][:,:]
    SourceCoord.lat = nc.variables['gphit'][:,:]
    nc.close()

    DstCoord.lonlat = {'t': {}, 'u': {}, 'v': {}}

    nc = Dataset(Setup.settings['dst_hgr'], 'r')
    """
    DstCoord.lonlat['t']['lon'] = nc.variables['glamt'][0,:,:]
    DstCoord.lonlat['t']['lat'] = nc.variables['gphit'][0,:,:]
    DstCoord.lonlat['u']['lon'] = nc.variables['glamu'][0,:,:]
    DstCoord.lonlat['u']['lat'] = nc.variables['gphiu'][0,:,:]
    DstCoord.lonlat['v']['lon'] = nc.variables['glamv'][0,:,:]
    DstCoord.lonlat['v']['lat'] = nc.variables['gphiv'][0,:,:]
    """

    # Read and assign netcdf data
    for g in 't', 'u', 'v':
        DstCoord.lonlat[g]['lon'] = nc.variables['glam' + g][0,:,:]
        DstCoord.lonlat[g]['lat'] = nc.variables['gphi' + g][0,:,:]
    nc.close()

    
    # Set of 
    DstCoord.bdy_lonlat = {'t': {}, 'u': {}, 'v': {}, 'z': {}}
    for g in 't', 'u', 'v', 'z':
        for l in 'lon', 'lat':
            DstCoord.bdy_lonlat[g][l] = np.zeros(num_bdy[g])


    for g in 't', 'u', 'v':
        print 'range is ', num_bdy[g]
        for i in range(num_bdy[g]):
            x = bdy_ind[g][0][i,1]
            y = bdy_ind[g][0][i,0]
            DstCoord.bdy_lonlat[g]['lon'][i] = DstCoord.lonlat[g]['lon'][x,y]
            DstCoord.bdy_lonlat[g]['lat'][i] = DstCoord.lonlat[g]['lat'][x,y]

        DstCoord.lonlat[g]['lon'][DstCoord.lonlat[g]['lon'] > 180] -= 360

    DstCoord.bdy_lonlat['z']['lon'] = DstCoord.bdy_lonlat['t']['lon'][bdy_r == 1]
    DstCoord.bdy_lonlat['z']['lat'] = DstCoord.bdy_lonlat['t']['lat'][bdy_r == 1]

    print DstCoord.bdy_lonlat['t']['lon'].shape


    # TEMP


     #################
    ## Set Constants ##
     #################
    
    # Create a Time handler object
    SourceTime = source_time.SourceTime(Setup.settings['src_dir'])
    acc = -3168000 - 86400 # Straight out Matlab code. Question this. 

    # Extract source data on dst grid
    # Assign time info to Boundary T, U, V objects
    # FIX MY STRUCTURING
    Grid_T.source_time, Grid_T.dir_list = SourceTime.get_source_time('t', acc)
    Grid_U.source_time, Grid_U.dir_list = SourceTime.get_source_time('u', acc)
    Grid_V.source_time, Grid_V.dir_list = SourceTime.get_source_time('v', acc)
    if ice:
        Grid_ice.source_time, Grid_ice.dir_list = SourceTime.get_source_time('i', acc)

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

    unit_string = 'seconds since %d-01-01 00:00:00' %Setup.settings['base_year']
    unit_origin = '%d-01-01 00:00:00' %Setup.settings['base_year']
    
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
    year = 1979
    month = 11



    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    print '\n- Old Callington Ratsmith should never be trusted -\n'
    print 'DstCoord: \n', dir(DstCoord), '\nSourceCoord: \n', dir(SourceCoord)

    ret1 = {'t': Grid_T, 'u':Grid_U, 'v':Grid_V, 'f':Grid_F}
    return Setup, ret1, SourceCoord, DstCoord

#go()
