# ===================================================================
# The contents of this file are dedicated to the public domain.  To
# the extent that dedication to the public domain is not available,
# everyone is granted a worldwide, perpetual, royalty-free,
# non-exclusive license to exercise all rights associated with the
# contents of this file for any purpose whatsoever.
# No rights are reserved.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# ===================================================================


'''
###
## This is the equivalent to the matlab of same name.
## The driver thast makes things happen. bit messy atm
## Note: Individual matlab vars are grouped into relevant dictionaries
## for clarity/ease of use
## Algorithm:
## 1) Interpolate the data from the global grid to local grid
## 2) then interpolate the local grid data time data to input time resolution
## Written by John Kazimierz Farey, Sep 2012
## Port of Matlab code of the Cardinal James Harle
## @author James Harle
## @author John Kazimierz Farey
## $Last commit on:$
###
'''

# pylint: disable=E1103
# pylint: disable=no-name-in-module

#External imports
from time import clock
from calendar import monthrange
import numpy as np
from scipy.interpolate import interp1d
import logging
from pynemo.reader.factory import GetFile
from pynemo.reader.factory import GetReader
from netcdftime import datetime
import copy

#local imports
from pynemo import nemo_bdy_setup as setup
from pynemo import nemo_bdy_gen_c as gen_grid
from pynemo import nemo_coord_gen_pop as coord
from pynemo import nemo_bdy_zgrv2 as zgrv
from pynemo import nemo_bdy_src_time as source_time
from pynemo.reader import factory

from pynemo import nemo_bdy_source_coord as source_coord
from pynemo import nemo_bdy_dst_coord as dst_coord
from pynemo import nemo_bdy_ice
from pynemo import nemo_bdy_extr_tm3
from pynemo.tide import nemo_bdy_tide3
from pynemo.tide import nemo_bdy_tide_ncgen
from pynemo import pynemo_settings_editor
from pynemo.utils import Constants

from pynemo import nemo_bdy_ncgen
from pynemo import nemo_bdy_ncpop
from pynemo.gui.nemo_bdy_mask import Mask as Mask_File
from PyQt4.QtGui import QMessageBox
#import pickle

class Grid(object):
    """ A Grid object that stores bdy grid information """    
    def __init__(self):
        self.bdy_i = None
        self.bdy_r = None
        self.grid_type = None
        self.max_i = None
        self.max_j = None
        self.fname_2 = None
        self.source_time = None

logger = logging.getLogger(__name__)
def process_bdy(setup_filepath=0, mask_gui=False):
    """ Main entry to the processing of the bdy 
    Keyword arguments:
    setup_filepath -- file path to bdy file
    mask_gui -- whether gui to select the mask file needs to be poped up
    """
    #Logger
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

    start = clock()
    bdy_msk = _get_mask(Setup, mask_gui)

    logger.info(clock() - start)
    logger.info('Done Mask')

    DstCoord.bdy_msk = bdy_msk == 1

    start = clock()
    logger.info('start bdy_t')
    grid_t = gen_grid.Boundary(bdy_msk, settings, 't')

    logger.info(clock() - start)
    start = clock()
    logger.info('start bdy_u')
    grid_u = gen_grid.Boundary(bdy_msk, settings, 'u')
    logger.info('start bdy_v')

    logger.info(clock() - start)
    start = clock()
    grid_v = gen_grid.Boundary(bdy_msk, settings, 'v')
    logger.info('start bdy_f')

    logger.info(clock() - start)
    start = clock()

    grid_f = gen_grid.Boundary(bdy_msk, settings, 'f')
    logger.info('done  bdy t,u,v,f')

    logger.info(clock() - start)
    start = clock()
    if ice:
        grid_ice = nemo_bdy_ice.BoundaryIce()
        grid_ice.grid_type = 't'
        grid_ice.bdy_r = grid_t.bdy_r

    bdy_ind = {'t': grid_t, 'u': grid_u, 'v': grid_v, 'f': grid_f}

    for k in bdy_ind.keys():
        logger.info('bdy_ind %s %s %s', k, bdy_ind[k].bdy_i.shape, bdy_ind[k].bdy_r.shape)

    start = clock()
    co_set = coord.Coord(settings['dst_dir']+'/coordinates.bdy.nc', bdy_ind)
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

    grid_u.bdy_i = grid_u.bdy_i[grid_u.bdy_r == 0, :]
    grid_v.bdy_i = grid_v.bdy_i[grid_v.bdy_r == 0, :]
    # Ought to look at this
    bdy_r = grid_t.bdy_r

    # number of points
    num_bdy = {}
    num_bdy['t'] = len(grid_t.bdy_i[:, 0])
    num_bdy['z'] = len(grid_t.bdy_i[bdy_r == 0, 0])
    num_bdy['u'] = len(grid_u.bdy_i[:, 0])
    num_bdy['v'] = len(grid_v.bdy_i[:, 0])

    ######
    # Gather grid info
    ######
    logger.info('gather grid info')
    start = clock()
    # how grid broken up, long/lats
    nc = GetFile(settings['src_zgr'])
    SourceCoord.zt = nc['gdept_0'][:]
    nc.close()

    logger.info(clock() - start)
    logger.info('generating depth info')
    start = clock()
    z = zgrv.Depth(grid_t.bdy_i, grid_u.bdy_i, grid_v.bdy_i, settings)

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
    nc = GetFile(settings['src_hgr'])
    SourceCoord.lon = nc['glamt'][:,:]
    SourceCoord.lat = nc['gphit'][:,:]
    try: # if they are masked array convert them to normal arrays
        SourceCoord.lon = SourceCoord.lon.filled()
    except:
        pass
    try:
        SourceCoord.lat = SourceCoord.lat.filled()
    except:
        pass
        
    nc.close()

    logger.info(clock() - start)
    DstCoord.lonlat = {'t': {}, 'u': {}, 'v': {}}

    start = clock()
    nc = GetFile(settings['dst_hgr'])

    logger.info('read and assign netcdf data, lon lat')
    # Read and assign netcdf data
    for g in 't', 'u', 'v':
        DstCoord.lonlat[g]['lon'] = nc['glam' + g][0, :, :]
        DstCoord.lonlat[g]['lat'] = nc['gphi' + g][0, :, :]
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
#    SourceTime = source_time.SourceTime(settings['src_dir'])
    acc = settings['src_time_adj']
#    acc = -3168000 - 86400 # Straight out Matlab code. Question this.
#    acc = 0
    # Extract source data on dst grid
    # Assign time info to Boundary T, U, V objects
    # FIX MY STRUCTURING
#    grid_t.source_time = SourceTime.get_source_time('t', acc)
#    grid_u.source_time = SourceTime.get_source_time('u', acc)
#    grid_v.source_time = SourceTime.get_source_time('v', acc)
#    grid_t.source_time = factory.GetRepository(settings['src_dir'], 't', acc).grid_source_data['t']
#    grid_u.source_time = factory.GetRepository(settings['src_dir'], 'u', acc).grid_source_data['u']
#    grid_v.source_time = factory.GetRepository(settings['src_dir'], 'v', acc).grid_source_data['v']
    reader = factory.GetReader(settings['src_dir'],acc)
    grid_t.source_time = reader['t']
    grid_u.source_time = reader['u']
    grid_v.source_time = reader['v']
    
    if ice:
#        grid_ice.source_time = SourceTime.get_source_time('i', acc)
#        grid_ice.source_time = factory.GetRepository(settings['src_dir'], 'i', acc).grid_source_data['i']
        grid_ice.source_time = reader['i']

    logger.debug(clock() - start)
    """
    # List source filenames

    for grid in [grid_t, grid_u, grid_v]:
        grid.dir_list = grid.get_dir_list(Setup.settings['src_dir'])
    # This might not be best way
    if Setup.settings['ice']:
        Grid_ice.dir_list = grid_t.get_dir_list(Setup.settings['src_dir'], ice=True)
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
    if settings['tide']:
        cosz, sinz, cosu, sinu, cosv, sinv = nemo_bdy_tide3.nemo_bdy_tpx7p2_rot(Setup, DstCoord,
                                                                                grid_t, grid_u, grid_v,
                                                                                settings["clname"])
        write_tidal_data(Setup, DstCoord, grid_t, num_bdy, settings["clname"], cosz, sinz)
        write_tidal_data(Setup, DstCoord, grid_u, num_bdy, settings["clname"], cosu, sinu)
        write_tidal_data(Setup, DstCoord, grid_v, num_bdy, settings["clname"], cosv, sinv)

    # Set the Year and month range
    start_year = settings["year_000"]
    end_year = settings["year_end"]
    if start_year > end_year:
        logging.error("Please check the nn_year_000 and nn_year_end values in input bdy file")
        return
    
    years = range(start_year, end_year+1)
    months = []
    if end_year - start_year >= 1:
        months = range(1,13)
    else:
        start_month = settings["month_000"]
        end_month = settings["month_end"]
        if end_month > 12 or start_month < 1:
            logging.error("Please check the nn_month_000 and nn_month_end values in input bdy file")
            return
        months = range(start_month, end_month+1)
    #enter the loop for each year and month extraction
    for year in years:
        for month in months:

            if Setup.settings['tra']:
                extract_t = nemo_bdy_extr_tm3.Extract(Setup.settings, SourceCoord, DstCoord,
                                                      grid_t, ['votemper', 'vosaline'],
                                                      years=year, months=month)
                extract_t.extract_month(year, month)
                #Get Date as a Number used in interpolation
                time_counter = np.zeros([len(extract_t.sc_time.time_counter)])
                for i in range(0, len(extract_t.sc_time.time_counter)):
                    time_counter[i] = extract_t.convert_date_to_destination_num(extract_t.sc_time.date_counter[i])

                date_start = datetime(year, month, 1, 0, 0, 0)
                date_end = datetime(year, month, monthrange(year, month)[1], 24, 60, 60)
                time_num_start = extract_t.convert_date_to_destination_num(date_start)
                time_num_end = extract_t.convert_date_to_destination_num(date_end)

                if monthrange(year, month)[1] == 29:
                    interp1fn = interp1d(np.arange(0, len(extract_t.sc_time.time_counter), 1), time_counter)
                    time_counter = interp1fn(np.append(np.arange(0,
                                                                 len(extract_t.sc_time.time_counter)-2+0.2,
                                                                 0.2),
                                                       np.arange(len(extract_t.sc_time.time_counter)-2+1/6.0,
                                                                 len(extract_t.sc_time.time_counter)-1+1/6.0,
                                                                 1/6.0)))
                                            #added 0.2 to include boundary
                else:
                    interp1fn = interp1d(np.arange(0, len(extract_t.sc_time.time_counter), 1), time_counter)
                    #added 0.2 to include boundary
                    time_counter = interp1fn(np.arange(0, len(extract_t.sc_time.time_counter)-1+0.2, 0.2))

                ft = np.where(((time_counter >= time_num_start) & (time_counter <= time_num_end)))
                time_counter = time_counter[ft]
                # interpolate the data to daily period in a month
                interpolate_data(extract_t, year, month, ft)

                output_filename_t = Setup.settings['dst_dir']+Setup.settings['fn']+ \
                                    '_bdyT_y'+str(year)+'m'+"%02d" % month+'.nc'
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

                nemo_bdy_ncpop.write_data_to_file(output_filename_t, 'votemper',
                                                  extract_t.d_bdy['votemper'][year]['data'])
                nemo_bdy_ncpop.write_data_to_file(output_filename_t, 'vosaline',
                                                  extract_t.d_bdy['vosaline'][year]['data'])

                nemo_bdy_ncpop.write_data_to_file(output_filename_t, 'bdy_msk', DstCoord.bdy_msk)
                nemo_bdy_ncpop.write_data_to_file(output_filename_t, 'deptht',
                                                  DstCoord.depths['t']['bdy_z'])
                nemo_bdy_ncpop.write_data_to_file(output_filename_t, 'nav_lon',
                                                  DstCoord.lonlat['t']['lon'])
                nemo_bdy_ncpop.write_data_to_file(output_filename_t, 'nav_lat',
                                                  DstCoord.lonlat['t']['lat'])
                #Added +1 to index nbidta, nbjdta, nbrdta to match the index of matlab
                nemo_bdy_ncpop.write_data_to_file(output_filename_t, 'nbidta',
                                                  grid_t.bdy_i[:, 0]+1)
                nemo_bdy_ncpop.write_data_to_file(output_filename_t, 'nbjdta',
                                                  grid_t.bdy_i[:, 1]+1)
                nemo_bdy_ncpop.write_data_to_file(output_filename_t, 'nbrdta',
                                                  grid_t.bdy_r[:]+1)
                nemo_bdy_ncpop.write_data_to_file(output_filename_t, 'time_counter', time_counter)

                if Setup.settings['ice']:
                    extract_write_ice_data(Setup, SourceCoord, DstCoord, grid_ice, year, month,
                                           ft, num_bdy, time_counter, unit_origin)
                #-------------------------bt------------------------------------------------------
                extract_write_bt_data(Setup, SourceCoord, DstCoord, grid_t, year, month,
                                      ft, num_bdy, time_counter, unit_origin)
            #Eco

            #U, V
            if Setup.settings['dyn2d'] or Setup.settings['dyn3d']:
                #---------------------------------------------U-----------------------------------
                grid_u2 = Grid()
                grid_u2.grid_type = 'u'
                grid_u2.bdy_i = grid_u.bdy_i
                grid_u2.max_i = DstCoord.lonlat['t']['lon'].shape[1]
                grid_u2.max_j = DstCoord.lonlat['t']['lon'].shape[0]
                grid_u2.fname_2 = grid_v.source_time
                extract_write_u_data(Setup, SourceCoord, DstCoord, grid_u, grid_u2, year, month,
                                     ft, num_bdy, time_counter, unit_origin)
                #----------------------------------------------V----------------------------------
                grid_v2 = Grid()
                grid_v2.grid_type = 'v'
                grid_v2.bdy_i = grid_v.bdy_i
                grid_v2.max_i = DstCoord.lonlat['t']['lon'].shape[1]
                grid_v2.max_j = DstCoord.lonlat['t']['lon'].shape[0]
                grid_v2.fname_2 = grid_v.source_time
                grid_uv = Grid()
                grid_uv.bdy_i = grid_v.bdy_i
                grid_uv.grid_type = 'v'
                grid_uv.bdy_r = grid_v.bdy_r
                grid_uv.source_time = grid_u.source_time
                extract_write_v_data(Setup, SourceCoord, DstCoord, grid_uv, grid_v2, year, month,
                                     ft, num_bdy, time_counter, unit_origin)


def _get_mask(Setup, mask_gui):
    """ This method reads the mask information from the netcdf file or opens a gui
    to create a mask depending on the mask_gui input. return the mask data. The default mask
    data is using bathymetry and applying a 1px halo
    Keyword arguments:
    Setup -- settings for bdy
    mask_gui -- boolean to open mask gui.
    """
    bdy_msk = None
    if mask_gui:
        #Open the gui to create a mask
        _, mask = pynemo_settings_editor.open_settings_dialog(Setup)
        bdy_msk = mask.data
        Setup.refresh()
	logger.info('Using GUI defined mask')
    else:
        try:
            #mask filename and mask file flag is set
            if Setup.bool_settings['mask_file'] and Setup.settings['mask_file'] is not None:
                mask = Mask_File(Setup.settings['bathy'], Setup.settings['mask_file'])
                bdy_msk = mask.data
		logger.info('Using input mask file')
            elif Setup.bool_settings['mask_file']:
                logger.error("Mask file is not given")
                return
            else: #no mask file specified then use default 1px halo mask
                logger.warning("Using default mask with bathymetry!!!!")
                mask = Mask_File(Setup.settings['bathy'])
                mask.apply_border_mask(Constants.DEFAULT_MASK_PIXELS)
                bdy_msk = mask.data
        except:
            return
    if np.amin(bdy_msk) == 0:
        # Mask is not set throw a warning message and set border to 1px.
        logger.warning("Setting the mask to 1px border")
        QMessageBox.warning(None,"pyNEMO", "Mask is not set, setting a 1 pixel border mask")
        if bdy_msk is not None and 1 < bdy_msk.shape[0] and 1 < bdy_msk.shape[1]:
            tmp = np.ones(bdy_msk.shape, dtype=bool)
            tmp[1:-1, 1:-1] = False
            bdy_msk[tmp] = -1
            
    return bdy_msk
    
def extract_write_ice_data(Setup, SourceCoord, DstCoord, grid_ice, year, month, ft, num_bdy,
                           time_counter, unit_origin):
    """ writes the ice data to file"""
    source_coord_var_copy = copy.deepcopy(SourceCoord)
    source_coord_var_copy.zt = np.zeros([1, 1])
    dst_coord_var_copy = copy.deepcopy(DstCoord)
    dst_coord_var_copy.depths['t']['bdy_z'] = np.zeros([1, 1])
    extract_ice = nemo_bdy_extr_tm3.Extract(Setup.settings, source_coord_var_copy,
                                            dst_coord_var_copy, grid_ice,
                                            ['iicethic', 'ileadfra', 'isnowthi'],
                                            years=year, months=month)
    extract_ice.extract_month(year, month)
    interpolate_data(extract_ice, year, month, ft) #interpolate the data to daily period in a month

    output_filename_t = Setup.settings['dst_dir']+Setup.settings['fn']+'_bdyT_y'+str(year)+ \
                        'm'+"%02d" % month+'.nc'
    logger.info('Outputing ICE values to T file:%s', output_filename_t)
    nemo_bdy_ncpop.write_data_to_file(output_filename_t, 'iicethic',
                                      extract_ice.d_bdy['iicethic'][year]['data'])
    nemo_bdy_ncpop.write_data_to_file(output_filename_t, 'ileadfra',
                                      extract_ice.d_bdy['ileadfra'][year]['data'])
    nemo_bdy_ncpop.write_data_to_file(output_filename_t, 'isnowthi',
                                      extract_ice.d_bdy['isnowthi'][year]['data'])

def extract_write_bt_data(Setup, SourceCoord, DstCoord, grid_t, year, month,
                          ft, num_bdy, time_counter, unit_origin):
    """Writes BT data file"""
    SourceCoordLatLon = copy.deepcopy(SourceCoord)
    SourceCoordLatLon.zt = np.zeros([1, 1])
    DstCoordLatLon = copy.deepcopy(DstCoord)
    DstCoordLatLon.depths['t']['bdy_z'] = np.zeros([1, 1])
    grid_t.grid_type = 'z'
    extract_t = nemo_bdy_extr_tm3.Extract(Setup.settings,
                                          SourceCoordLatLon,
                                          DstCoord,
                                          grid_t,
                                          ['sossheig'],
                                          years=year,
                                          months=month)
    extract_t.extract_month(year, month)
    interpolate_data(extract_t, year, month, ft) #interpolate the data to daily period in a month
    output_filename_bt = Setup.settings['dst_dir']+ \
                         Setup.settings['fn']+ \
                         '_bt_bdyT_y'+str(year)+'m'+"%02d" % month+'.nc'
    logger.info('Outputting bt_T file: %s', output_filename_bt)
    nemo_bdy_ncgen.CreateBDYNetcdfFile(output_filename_bt, num_bdy['z'],
                                       DstCoord.lonlat['t']['lon'].shape[1],
                                       DstCoord.lonlat['t']['lon'].shape[0],
                                       1, Setup.settings['rimwidth'],
                                       Setup.settings['dst_metainfo'],
                                       unit_origin, Setup.settings['fv'],
                                       Setup.settings['dst_calendar'], 'Z')
    nemo_bdy_ncpop.write_data_to_file(output_filename_bt, 'sossheig',
                                      extract_t.d_bdy['sossheig'][year]['data'])
    nemo_bdy_ncpop.write_data_to_file(output_filename_bt, 'bdy_msk', DstCoord.bdy_msk)
    nemo_bdy_ncpop.write_data_to_file(output_filename_bt, 'nav_lon', DstCoord.lonlat['t']['lon'])
    nemo_bdy_ncpop.write_data_to_file(output_filename_bt, 'nav_lat', DstCoord.lonlat['t']['lat'])
    nemo_bdy_ncpop.write_data_to_file(output_filename_bt, 'nbidta',
                                      grid_t.bdy_i[np.where(grid_t.bdy_r == 0), 0]+1)
    nemo_bdy_ncpop.write_data_to_file(output_filename_bt, 'nbjdta',
                                      grid_t.bdy_i[np.where(grid_t.bdy_r == 0), 1]+1)
    nemo_bdy_ncpop.write_data_to_file(output_filename_bt, 'nbrdta',
                                      grid_t.bdy_r[np.where(grid_t.bdy_r == 0)]+1)
    nemo_bdy_ncpop.write_data_to_file(output_filename_bt, 'time_counter', time_counter)

def extract_write_u_data(setup_var, source_coord_var, dst_coord_var, grid_u, grid_u2, year,
                         month, ft, num_bdy, time_counter, unit_origin):
    """ Extracts the U component of velocity data, interpolates the data for a month and writes
      data to netCDF file"""
    extract_u = nemo_bdy_extr_tm3.Extract(setup_var.settings, source_coord_var, dst_coord_var,
                                          grid_u, ['vozocrtx', 'vomecrty'], Grid_2=grid_u2,
                                          years=year, months=month)
    extract_u.extract_month(year, month)
    interpolate_data(extract_u, year, month, ft)
    nanindex = np.isnan(extract_u.d_bdy['vozocrtx'][year]['data'])
    tmp_vozocrtx = extract_u.d_bdy['vozocrtx'][year]['data']
    tmp_vozocrtx[nanindex] = setup_var.settings['fv']
    output_filename_u = setup_var.settings['dst_dir'] + setup_var.settings['fn'] + '_bdyU_y' + \
                        str(year) + 'm' + "%02d" % month + '.nc'
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
    nemo_bdy_ncpop.write_data_to_file(output_filename_u, 'vozocrtx', tmp_vozocrtx)
    nemo_bdy_ncpop.write_data_to_file(output_filename_u, 'depthu',
                                      dst_coord_var.depths['u']['bdy_z'])
    # TODO: write vobtcrtx
    tmp_vozocrtx[nanindex] = float('NaN')
    # ft size has some problem
    tmp_tile = np.tile(dst_coord_var.depths['u']['bdy_dz'], [ft[0].size, 1, 1, 1])
    tmp_vobtcrtx = np.nansum(np.reshape(tmp_vozocrtx[:, :-1, :],
                                        tmp_tile.shape) * tmp_tile, 2) / np.nansum(tmp_tile, 2)
    tmp_vobtcrtx[np.isnan(tmp_vobtcrtx)] = setup_var.settings['fv']
    nemo_bdy_ncpop.write_data_to_file(output_filename_u, 'vobtcrtx', tmp_vobtcrtx)
    nemo_bdy_ncpop.write_data_to_file(output_filename_u, 'nav_lon',
                                      dst_coord_var.lonlat['u']['lon'])
    nemo_bdy_ncpop.write_data_to_file(output_filename_u, 'nav_lat',
                                      dst_coord_var.lonlat['u']['lat'])
    nemo_bdy_ncpop.write_data_to_file(output_filename_u, 'nbidta', grid_u.bdy_i[:, 0]+1)
    nemo_bdy_ncpop.write_data_to_file(output_filename_u, 'nbjdta', grid_u.bdy_i[:, 1]+1)
    nemo_bdy_ncpop.write_data_to_file(output_filename_u, 'nbrdta', np.ones([1, num_bdy['u']]))
    nemo_bdy_ncpop.write_data_to_file(output_filename_u, 'time_counter', time_counter)
    #-------------------------------------------End U-------------------------------

def extract_write_v_data(setup_var, source_coord_var, dst_coord_var, grid_v, grid_v2, year, month,
                         ft, num_bdy, time_counter, unit_origin):
    """ Extracts the V component of velocity data, interpolates the data for a month and writes
      data to netCDF file"""
    extract_v = nemo_bdy_extr_tm3.Extract(setup_var.settings, source_coord_var, dst_coord_var,
                                          grid_v, ['vozocrtx', 'vomecrty'], Grid_2=grid_v2,
                                          years=year, months=month)
    extract_v.extract_month(year, month) #return values in vozocrtx instead of vomecrty
    interpolate_data(extract_v, year, month, ft)
    nanindex = np.isnan(extract_v.d_bdy['vozocrtx'][year]['data'])
    tmp_vozocrty = extract_v.d_bdy['vozocrtx'][year]['data']
    tmp_vozocrty[nanindex] = setup_var.settings['fv']

    output_filename_v = setup_var.settings['dst_dir']+setup_var.settings['fn']+ \
                        '_bdyV_y'+str(year)+ 'm'+"%02d" % month+'.nc'
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
    nemo_bdy_ncpop.write_data_to_file(output_filename_v, 'vomecrty',
                                      extract_v.d_bdy['vozocrtx'][year]['data'])
    nemo_bdy_ncpop.write_data_to_file(output_filename_v, 'depthv',
                                      dst_coord_var.depths['v']['bdy_z'])
    tmp_vozocrty[nanindex] = float('NaN')
    #ft size has some problem
    tmp_tile = np.tile(dst_coord_var.depths['v']['bdy_dz'], [ft[0].size, 1, 1, 1])
    tmp_vobtcrty = np.nansum(np.reshape(tmp_vozocrty[:, :-1, :], tmp_tile.shape)*tmp_tile, 2) \
                   /np.nansum(tmp_tile, 2)
    tmp_vobtcrty[np.isnan(tmp_vobtcrty)] = setup_var.settings['fv']
    nemo_bdy_ncpop.write_data_to_file(output_filename_v, 'vobtcrty', tmp_vobtcrty)
    nemo_bdy_ncpop.write_data_to_file(output_filename_v, 'nav_lon',
                                      dst_coord_var.lonlat['u']['lon'])
    nemo_bdy_ncpop.write_data_to_file(output_filename_v, 'nav_lat',
                                      dst_coord_var.lonlat['u']['lat'])
    nemo_bdy_ncpop.write_data_to_file(output_filename_v, 'nbidta', grid_v.bdy_i[:, 0]+1)
    nemo_bdy_ncpop.write_data_to_file(output_filename_v, 'nbjdta', grid_v.bdy_i[:, 1]+1)
    nemo_bdy_ncpop.write_data_to_file(output_filename_v, 'nbrdta', np.ones([1, num_bdy['v']]))
    nemo_bdy_ncpop.write_data_to_file(output_filename_v, 'time_counter', time_counter)

def write_tidal_data(setup_var, dst_coord_var, grid, num_bdy, tide_cons, cos_g, sin_g):
    """ This method writes the tidal data to netcdf file"""
    indx = 0
    for tide_con in tide_cons:
        const_name = setup_var.settings['clname'][tide_con]
        const_name = const_name.replace("'", "").upper()

        if grid.grid_type == 't':
            output_filename_tide = setup_var.settings['dst_dir']+setup_var.settings['fn']+ \
                                   "_bdytide_rotT_"+const_name+"_grid_"+ \
                                   grid.grid_type.upper()+".nc"
            nemo_bdy_tide_ncgen.CreateBDYTideNetcdfFile(output_filename_tide, num_bdy['z'],
                                                        dst_coord_var.lonlat['t']['lon'].shape[1],
                                                        dst_coord_var.lonlat['t']['lon'].shape[0],
                                                        "tidal elevation components for:"+tide_con,
                                                        setup_var.settings['fv'], 'T')
            nemo_bdy_ncpop.write_data_to_file(output_filename_tide, 'z1', cos_g[indx])
            nemo_bdy_ncpop.write_data_to_file(output_filename_tide, 'z2', sin_g[indx])
            nemo_bdy_ncpop.write_data_to_file(output_filename_tide, 'bdy_msk',
                                              dst_coord_var.bdy_msk)
            nemo_bdy_ncpop.write_data_to_file(output_filename_tide, 'nav_lon',
                                              dst_coord_var.lonlat['t']['lon'])
            nemo_bdy_ncpop.write_data_to_file(output_filename_tide, 'nav_lat',
                                              dst_coord_var.lonlat['t']['lat'])
            nemo_bdy_ncpop.write_data_to_file(output_filename_tide, 'nbidta',
                                              grid.bdy_i[np.where(grid.bdy_r == 0), 0]+1)
            nemo_bdy_ncpop.write_data_to_file(output_filename_tide, 'nbjdta',
                                              grid.bdy_i[np.where(grid.bdy_r == 0), 1]+1)
            nemo_bdy_ncpop.write_data_to_file(output_filename_tide, 'nbrdta',
                                              grid.bdy_r[np.where(grid.bdy_r == 0)]+1)
        elif grid.grid_type == 'u':
            output_filename_tide = setup_var.settings['dst_dir']+setup_var.settings['fn']+ \
                                   "_bdytide_rotT_"+const_name+"_grid_"+ \
                                   grid.grid_type.upper()+".nc"
            nemo_bdy_tide_ncgen.CreateBDYTideNetcdfFile(output_filename_tide, num_bdy['u'],
                                                        dst_coord_var.lonlat['t']['lon'].shape[1],
                                                        dst_coord_var.lonlat['t']['lon'].shape[0],
                                                        "tidal east velocity components for:" \
                                                        +tide_con,
                                                        setup_var.settings['fv'], 'U')
            nemo_bdy_ncpop.write_data_to_file(output_filename_tide, 'u1', cos_g[indx])
            nemo_bdy_ncpop.write_data_to_file(output_filename_tide, 'u2', sin_g[indx])
            nemo_bdy_ncpop.write_data_to_file(output_filename_tide, 'nav_lon',
                                              dst_coord_var.lonlat['u']['lon'])
            nemo_bdy_ncpop.write_data_to_file(output_filename_tide, 'nav_lat',
                                              dst_coord_var.lonlat['u']['lat'])
            nemo_bdy_ncpop.write_data_to_file(output_filename_tide, 'nbidta',
                                              grid.bdy_i[:, 0]+1)
            nemo_bdy_ncpop.write_data_to_file(output_filename_tide, 'nbjdta',
                                              grid.bdy_i[:, 1]+1)
            nemo_bdy_ncpop.write_data_to_file(output_filename_tide, 'nbrdta',
                                              np.ones([1, num_bdy['u']]))
        elif grid.grid_type == 'v':
            output_filename_tide = setup_var.settings['dst_dir']+setup_var.settings['fn']+ \
                                   "_bdytide_rotT_"+const_name+"_grid_"+ \
                                   grid.grid_type.upper()+".nc"
            nemo_bdy_tide_ncgen.CreateBDYTideNetcdfFile(output_filename_tide, num_bdy['v'],
                                                        dst_coord_var.lonlat['t']['lon'].shape[1],
                                                        dst_coord_var.lonlat['t']['lon'].shape[0],
                                                        "tidal north velocity components for:"+ \
                                                        tide_con, setup_var.settings['fv'], 'V')
            nemo_bdy_ncpop.write_data_to_file(output_filename_tide, 'v1', cos_g[indx])
            nemo_bdy_ncpop.write_data_to_file(output_filename_tide, 'v2', sin_g[indx])
            nemo_bdy_ncpop.write_data_to_file(output_filename_tide, 'nav_lon',
                                              dst_coord_var.lonlat['v']['lon'])
            nemo_bdy_ncpop.write_data_to_file(output_filename_tide, 'nav_lat',
                                              dst_coord_var.lonlat['v']['lat'])
            nemo_bdy_ncpop.write_data_to_file(output_filename_tide, 'nbidta',
                                              grid.bdy_i[:, 0]+1)
            nemo_bdy_ncpop.write_data_to_file(output_filename_tide, 'nbjdta',
                                              grid.bdy_i[:, 1]+1)
            nemo_bdy_ncpop.write_data_to_file(output_filename_tide, 'nbrdta',
                                              np.ones([1, num_bdy['v']]))
        indx = indx+1

        
def interpolate_data(extract, year, month, time_indexes):
    """This method does a 1D interpolation of daily mean data of a month"""

    for variable_name in extract.d_bdy:
        lon_len = len(extract.d_bdy[variable_name][year]['data'][:, 0, 0])
        lon = np.squeeze(np.asarray(np.mat(np.arange(0, lon_len, 1)).transpose()))
        lat = extract.d_bdy[variable_name][year]['data'].transpose(0, 2, 1)
        if monthrange(year, month)[1] == 29: #leap year
            interp1fn = interp1d(lon, lat, axis=0)
            extra_axis = np.append(np.arange(0, lon_len-2+0.2, 0.2),
                                   np.arange(lon_len-2+1/6.0, lon_len-1+1/6.0, 1/6.0))
            extract.d_bdy[variable_name][year]['data'] = interp1fn(extra_axis)
        else:
            interp1fn = interp1d(lon, lat, axis=0)
            extract.d_bdy[variable_name][year]['data'] = interp1fn(np.squeeze(np.asarray(np.mat(np.arange(0, lon_len-1+0.2, 0.2)).transpose())))      #added 0.2 to include boundary
        tmp = np.squeeze(extract.d_bdy[variable_name][year]['data'][time_indexes, :, :])
        if len(tmp.shape) == 2:
            tmp = tmp.reshape(tmp.shape+(1L, ))
        extract.d_bdy[variable_name][year]['data'] = tmp.transpose(0, 2, 1)

