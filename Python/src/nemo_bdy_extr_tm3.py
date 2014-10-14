#
#
# Written by John Kazimierz Farey, Sep 2012
# Port of Matlab code of James Harle
#
#

'''
Count var is not implemented, I suggest handling multi entry files
via splitting up in src_time module

Initial arguments:
setup: dictionary
SourceCoord : source coordinates object
DstCoord: destination coordinates object
Grid: t, u, v or f. Decides how to run
var_nam: netcdf file variable names
Grid_2: Rotation grid, defaults to None if unneeded
fnames_2: source_time dictionary of secondary grid for rotation
years: year of data to be extracted
months: month of data to be extracted
'''


from calendar import monthrange, isleap
from datetime import datetime

import numpy as np
import scipy.spatial as sp
from netCDF4 import Dataset
from netcdftime import utime
#from netCDF4 import date2num, num2date

import nemo_bdy_grid_angle
from nemo_bdy_lib import sub2ind

import copy # DEBUG ONLY- allows multiple runs without corruption

class Extract:


                # SourceCoord, DstCoord, grid, Setup dict, 
    def __init__(self, setup, SourceCoord, DstCoord, Grid, var_nam, 
                 Grid_2=None, fnames_2=None, years=[1979], months=[11]):

        self.g_type = Grid.grid_type
        self.settings = setup
        SC = copy.deepcopy(SourceCoord)
        DC = copy.deepcopy(DstCoord)
        if self.g_type == 't':  # this should be removed as bdy_r needs to be defined for all grids
            bdy_r = copy.deepcopy(Grid.bdy_r)
        else:
            bdy_r = Grid.bdy_r[Grid.bdy_r == 0]

        self.fnames_2 = fnames_2 # second source times dict
        sc_time = Grid.source_time
        self.var_nam = var_nam
        sc_z = SC.zt[0]
        sc_z_len = len(sc_z)
        
        # We can infer rot_str and rot_dir
        if Grid_2:
            bdy_vel = Grid_2.bdy_i
            self.rot_str = Grid_2.grid_type
            if self.rot_str == 'u':
                self.rot_dir = 'i'
            elif self.rot_str == 'v':
                self.rot_dir = 'j'
            else:
                raise ValueError('Invalid rotation grid type: %s' 
                                 %self.rot_str)

        dst_lon = DC.bdy_lonlat[self.g_type]['lon']
        dst_lat = DC.bdy_lonlat[self.g_type]['lat']
        try:
            dst_dep = DC.depths[self.g_type]['bdy_z']
        except KeyError:
            dst_dep = np.zeros([1])
        self.isslab = len(dst_dep) == 1

        # ??? Should this be read from settings?
        wei_121 = np.array([0.5, 0.25, 0.25])

        #years = [1979]
        #months = [11] #12
        
        SC.lon = SC.lon.squeeze()
        SC.lat = SC.lat.squeeze()
        
        # Make function of grid resolution
        fr = 0.2

        # infer key_vec
        if Grid_2 is not None:
            self.key_vec = True
        else: 
            self.key_vec = False

        num_bdy = len(dst_lon)
        self.nvar = len(self.var_nam)
        if self.key_vec:
            self.nvar = self.nvar / 2 
 
        print 'nvar: ', self.nvar
        print 'key vec: ', self.key_vec

        # Find subset of source data set required to produce bdy points

        # Temp vars. Just clearer way of doing lots than all on one line
        ive = SC.lon < np.amax(dst_lon)
        got = SC.lon > np.amin(dst_lon)
        ivegot = np.logical_and(ive, got)
        two = SC.lat > np.amin(dst_lat)
        legs = SC.lat < np.amax(dst_lat)
        twolegs = np.logical_and(two, legs)

        dejavu = np.logical_and(ivegot, twolegs)
        dejavu = np.argwhere(dejavu.flatten(0) != 0)#.flatten(1)

        sub_j, sub_i = np.unravel_index(dejavu, SC.lon.shape)

        # Find I/J range
        imin = np.maximum(np.amin(sub_i) - 2, 0)
        imax = np.minimum(np.amax(sub_i) + 2, len(SC.lon[0,:]) - 1) + 1
        jmin = np.maximum(np.amin(sub_j) - 2, 0)
        jmax = np.minimum(np.amax(sub_j) + 2, len(SC.lon[:,0]) - 1) + 1

        print ' \n imin: %d\n imax: %d\n jmin: %d\n jmax: %d\n' %(imin,imax,jmin,jmax)
        # SC.lon: big matrix from an nc, dst_lon - bdy_lon_t
        SC.lon = SC.lon[jmin:jmax, imin:imax]
        SC.lat = SC.lat[jmin:jmax, imin:imax]

        # Init of gsin* and gcos* at first call for rotation of vectors
        if self.key_vec:
            T_GridAngles = nemo_bdy_grid_angle.GridAngle(
                       self.settings['src_hgr'], imin, imax, jmin, jmax, 't')
            RotStr_GridAngles = nemo_bdy_grid_angle.GridAngle(
                         self.settings['dst_hgr'], 0, imax, 0, jmax, self.rot_str)
            
            self.gcos = T_GridAngles.cosval
            self.gsin = T_GridAngles.sinval
            dst_gcos = RotStr_GridAngles.cosval
            dst_gsin = RotStr_GridAngles.sinval
            
            # BEWARE THE ETERNAL FORTRAN C CONFLICT
            dst_gcos = np.tile(np.array([dst_gcos]), (sc_z_len, 1, 1))
            dst_gsin = np.tile(np.array([dst_gsin]), (sc_z_len, 1, 1))

            # Retain only boundary points rotation information
            tmp_gcos = np.zeros((sc_z_len, bdy_vel.shape[0]))
            tmp_gsin = np.zeros((sc_z_len, bdy_vel.shape[0]))
            for p in range(bdy_vel.shape[0]):
                tmp_gcos[:,p] = dst_gcos[:, bdy_vel[p,1], bdy_vel[p,0]]
                tmp_gsin[:,p] = dst_gsin[:, bdy_vel[p,1], bdy_vel[p,0]]

            self.dst_gcos = tmp_gcos
            self.dst_gsin = tmp_gsin


        # Determine size of source data subset  
        if not self.isslab:  
            dst_len_z = len(dst_dep[:,0])
        else:
            dst_len_z = 0
        source_dims = SC.lon.shape

        # Find nearest neighbour on the source grid to each dst bdy point
        # Ann Query substitute
        source_tree = sp.cKDTree(zip(SC.lon.ravel(order='F'), 
                                 SC.lat.ravel(order='F')))
        dst_pts = zip(dst_lon[:].ravel(order='F'), dst_lat[:].ravel(order='F'))
        nn_dist, nn_id = source_tree.query(dst_pts, k=1)

        # Find surrounding points
        j_sp, i_sp = np.unravel_index(nn_id, source_dims, order='F')
        j_sp = np.vstack((j_sp, j_sp + 1, j_sp - 1))
        j_sp = np.vstack((j_sp, j_sp, j_sp))
        i_sp = np.vstack((i_sp, i_sp, i_sp))
        i_sp = np.vstack((i_sp, i_sp + 1, i_sp - 1))

        # FIG NOT IMPLEMENTED

        # Index out of bounds error check not implemented
        
        # Determine 9 nearest neighbours based on distance
        ind = sub2ind(source_dims, i_sp, j_sp)
        ind_rv = np.ravel(ind, order='F')
        sc_lon_rv = np.ravel(SC.lon, order='F')
        sc_lat_rv = np.ravel(SC.lat, order='F')
        sc_lon_ind = sc_lon_rv[ind_rv]
        
        diff_lon = sc_lon_ind - np.repeat(dst_lon,9).T
        diff_lon = diff_lon.reshape(ind.shape, order='F')
        out = np.abs(diff_lon) > 180
        diff_lon[out] = -np.sign(diff_lon[out]) * (360 - np.abs(diff_lon[out]))
        
        dst_lat_rep = np.repeat(dst_lat.T, 9)
        diff_lon_rv = np.ravel(diff_lon, order='F')
        dist_merid = diff_lon_rv * np.cos(dst_lat_rep * np.pi / 180)
        dist_zonal = sc_lat_rv[ind_rv] - dst_lat_rep
        
        dist_tot = np.power((np.power(dist_merid, 2) + 
                            np.power(dist_zonal, 2)), 0.5)
        dist_tot = dist_tot.reshape(ind.shape, order='F').T
        # Get sort inds, and sort
        dist_ind = np.argsort(dist_tot, axis=1, kind='mergesort')
        dist_tot = dist_tot[np.arange(dist_tot.shape[0])[:,None], dist_ind]

        # Shuffle ind to reflect ascending dist of source and dst points
        ind = ind.T
        for p in range(ind.shape[0]):
            ind[p,:] = ind[p,dist_ind[p,:]]

        if self.key_vec:
            self.gcos = self.gcos.flatten(1)[ind].reshape(ind.shape, order='F')
            self.gsin = self.gsin.flatten(1)[ind].reshape(ind.shape, order='F')
        
        sc_ind = {}
        sc_ind['ind'] = ind
        sc_ind['imin'], sc_ind['imax'] = imin, imax
        sc_ind['jmin'], sc_ind['jmax'] = jmin, jmax

        # Fig not implemented

        # Determine 1-2-1 filter indices
        id_121 = np.zeros((num_bdy, 3), dtype=np.int64)
        for r in range(int(np.amax(bdy_r))+1):
            print r            
            r_id = bdy_r != r
            rr_id = bdy_r == r
            tmp_lon = dst_lon.copy()
            tmp_lon[r_id] = -9999
            tmp_lat = dst_lat.copy()
            tmp_lat[r_id] = -9999
            source_tree = sp.cKDTree(zip(tmp_lon.ravel(order='F'),
                                     tmp_lat.ravel(order='F')))
            dst_pts = zip(dst_lon[rr_id].ravel(order='F'),
                          dst_lat[rr_id].ravel(order='F'))
            junk, an_id = source_tree.query(dst_pts, k=3, 
                                            distance_upper_bound=fr)    
            id_121[rr_id,:] = an_id
            id_121[id_121 == len(dst_lon)] = 0

        
        reptile = np.tile(id_121[:,0], 3).reshape(id_121.shape, order='F')
        id_121 = id_121 + reptile * (id_121 == 0)

        rep_dims = (id_121.shape[0], id_121.shape[1], sc_z_len)
        # These tran/tiles work like matlab. Tested with same Data.
        id_121 = id_121.repeat(sc_z_len).reshape(rep_dims).transpose(2,0,1)
        reptile = np.arange(sc_z_len).repeat(num_bdy).reshape(sc_z_len, 
                                                              num_bdy)
        reptile = reptile.repeat(3).reshape(num_bdy, 3, sc_z_len, 
                                            order='F').transpose(2,0,1)
  
        id_121 = sub2ind((sc_z_len, num_bdy), id_121, reptile)
        tmp_filt = wei_121.repeat(num_bdy).reshape(num_bdy, len(wei_121), 
                                                   order='F')
        tmp_filt = tmp_filt.repeat(sc_z_len).reshape(num_bdy, len(wei_121),
                                                     sc_z_len).transpose(2,0,1)

        # Fig not implemented
        
        if self.isslab != 1:
            # Determine vertical weights for the linear interpolation 
            # onto Dst grid
            # Allocate vertical index array
            dst_dep_rv = dst_dep.ravel(order='F')
            z_ind = np.zeros((num_bdy * dst_len_z, 2), dtype=np.int64)
            source_tree = sp.cKDTree(zip(sc_z.ravel(order='F')))
            junk, nn_id = source_tree.query(zip(dst_dep_rv), k=1)
            
            # WORKAROUND: the tree query returns out of range val when
            # dst_dep point is NaN, causing ref problems later.
            nn_id[nn_id == sc_z_len] = 0
            sc_z[nn_id]
            # Find next adjacent point in the vertical
            z_ind[:,0] = nn_id
            z_ind[sc_z[nn_id] > dst_dep_rv[:], 1] = nn_id[sc_z[nn_id] > 
                                                          dst_dep_rv[:]] - 1
            z_ind[sc_z[nn_id] <= dst_dep_rv[:], 1] = nn_id[sc_z[nn_id] <= 
                                                           dst_dep_rv[:]] + 1
            # Adjust out of range values
            z_ind[z_ind == -1] = 0
            z_ind[z_ind == sc_z_len] = sc_z_len - 1

            # Create weightings array
            sc_z[z_ind]
            z_dist = np.abs(sc_z[z_ind] - dst_dep.T.repeat(2).reshape(len(dst_dep_rv),2))
            rat = np.sum(z_dist, axis=1)
            z_dist = 1 - (z_dist / rat.repeat(2).reshape(len(rat), 2))

            # Update z_ind for the dst array dims and vector indexing
            # Replicating this part of matlab is difficult without causing
            # a Memory Error. This workaround may be +/- brilliant
            # In theory it maximises memory efficiency
            z_ind[:,:] += (np.arange(0, (num_bdy) * sc_z_len, sc_z_len)
                           [np.arange(num_bdy).repeat(2*dst_len_z)].reshape(
                                                                z_ind.shape))
        else:
            z_ind = np.zeros([1,1])
            z_dist = np.zeros([1,1])
        # End self.isslab
        
        # Set instance attributes
        self.first = True
        self.z_ind = z_ind
        self.z_dist = z_dist
        self.sc_ind = sc_ind
        self.dst_dep = dst_dep
        self.num_bdy = num_bdy
        self.id_121 = id_121
        if not self.isslab:
            self.bdy_z = DC.depths[self.g_type]['bdy_hbat']
        else:
            self.bdy_z = np.zeros([1])
        self.sc_z_len = sc_z_len
        self.sc_time = sc_time
        self.tmp_filt = tmp_filt
        self.dist_tot = dist_tot

        self.d_bdy = {}
        for v in range(self.nvar):
            self.d_bdy[self.var_nam[v]] = {}
            
            
        #    for year in years:
        #        self.d_bdy[self.var_nam[v]][year] = {'data': None, 'date': {}}

        
        # This loop should be removed. Once the module is initialised 
        # necessary attributes, the extract_month function should be called
        # externally 
        #
        # Extract source data level by level and interpolate onto bc points
        # Determine stretch factors for calendar out based on calendar in
        #for year in years:
        #    for month in months:
        #        self.extract_month(year, month)
        
# # # # # # # # # # # # # # # # # # # # #

     # # # # # # # # # 
    # # Functions # # 
   # # # # # # # # #    

    def extract_month(self, year, month):
        print '\n extract_month function called \n'
        # Check year entry exists in d_bdy, if not create it.
        for v in range(self.nvar):
            try:
                self.d_bdy[self.var_nam[v]][year]
            except KeyError:        
                self.d_bdy[self.var_nam[v]][year] = {'data': None, 'date': {}}
        
        i_run = np.arange(self.sc_ind['imin'], self.sc_ind['imax']) 
        j_run = np.arange(self.sc_ind['jmin'], self.sc_ind['jmax'])
        extended_i = np.arange(self.sc_ind['imin'] - 1, self.sc_ind['imax'])
        extended_j = np.arange(self.sc_ind['jmin'] - 1, self.sc_ind['jmax'])
        ind = self.sc_ind['ind']
        sc_time = self.sc_time
        sc_z_len = self.sc_z_len
 
        # define src/dst cals
        sf, ed = self.cal_trans(sc_time[0]['calendar'], 
                                self.settings['dst_calendar'], year, month)
        DstCal = utime('seconds since %d-1-1' %year, 
                       self.settings['dst_calendar'])
        dst_start = DstCal.date2num(datetime(year, month, 1))
        dst_end = DstCal.date2num(datetime(year, month, ed, 23, 59, 59))

        self.S_cal = utime(sc_time[0]['units'], sc_time[0]['calendar'])
#        self.D_cal = utime('seconds since %d-1-1' %year, 
#                           self.settings['dst_calendar'])
        self.D_cal = utime('seconds since %d-1-1' %self.settings['base_year'], 
                           self.settings['dst_calendar'])

        for date in sc_time:
            date['date_num'] = DstCal.date2num(date['date']) * sf

        # Get first and last date within range, init to cover entire range
        first_date = 0
        last_date = len(sc_time) - 1 
        rev_seq = range(len(sc_time))
        rev_seq.reverse()
        for date in rev_seq:
            if sc_time[date]['date_num'] < dst_start:
                first_date = date
                break
        for date in range(len(sc_time)):
            if sc_time[date]['date_num'] > dst_end:
                last_date = date
                break

        print 'first/last dates: ', first_date, last_date

        if self.first:
            nc_3 = Dataset(self.settings['src_msk'], 'r')
            varid_3 = nc_3.variables['tmask']
            t_mask = varid_3[:1, :sc_z_len, j_run, i_run]
            if self.key_vec:
                varid_3 = nc_3.variables['umask']
                u_mask = varid_3[:1, :sc_z_len, j_run, extended_i]
                varid_3 = nc_3.variables['vmask']
                v_mask = varid_3[:1, :sc_z_len, extended_j, i_run]
            nc_3.close()

        # Identify missing values and scale factors if defined
        meta_data = []
        meta_range = self.nvar
        if self.key_vec:
            meta_range += 1 
        for v in range(meta_range):
            meta_data.append({})
            for x in 'mv', 'sf', 'os', 'fv':
                meta_data[v][x] = np.ones((self.nvar, 1)) * np.NaN

        for v in range(self.nvar):
            meta_data[v] = self._get_meta_data(sc_time[first_date]['fname'], 
                                               self.var_nam[v], meta_data[v])

        if self.key_vec:
            n = self.nvar
            meta_data[n] = (self._get_meta_data(
                             self.fnames_2[first_date]['fname'], 
                             self.var_nam[n], meta_data[n]))

        # Loop over identified files
        for f in range(first_date, last_date + 1):
            sc_array = [None, None]
            sc_alt_arr = [None, None]
            print 'opening nc file: %s' %sc_time[f]['fname']
            nc = Dataset(sc_time[f]['fname'], 'r')
            # If extracting vector quantities open second file
            if self.key_vec:
                nc_2 = Dataset(self.fnames_2[f]['fname'], 'r')

            # Counters not implemented

            sc_bdy = np.zeros((len(self.var_nam), sc_z_len, ind.shape[0], 
                              ind.shape[1]))

            # Loop over time entries from file f
            for vn in range(self.nvar):
                # Extract sub-region of data
                print 'var_nam = %s' %self.var_nam[vn]
                varid = nc.variables[self.var_nam[vn]]
                # If extracting vector quantities open second var
                if self.key_vec:
                    varid_2 = nc_2.variables[self.var_nam[vn + 1]]

                # Extract 3D scalar variables
                if not self.isslab and not self.key_vec:
                    print ' 3D source array '
                    sc_array[0] = varid[:1 , :sc_z_len, j_run, i_run]
                # Extract 3D vector variables
                elif self.key_vec:
                    # For u vels take i-1
                    sc_alt_arr[0] = varid[:1, :sc_z_len, j_run, extended_i]
                    # For v vels take j-1
                    sc_alt_arr[1] = varid_2[:1, :sc_z_len, extended_j, i_run]
                # Extract 2D scalar vars
                else:
                    print ' 2D source array '
                    sc_array[0] = varid[:1, j_run, i_run].reshape([1,1,j_run.size,i_run.size])

                # Average vector vars onto T-grid
                if self.key_vec:
                    # First make sure land points have a zero val
                    sc_alt_arr[0] *= u_mask
                    sc_alt_arr[1] *= v_mask
                    # Average from to T-grid assuming C-grid stagger
                    sc_array[0] = 0.5 * (sc_alt_arr[0][:,:,:,:-1] + 
                                         sc_alt_arr[0][:,:,:,1:])
                    sc_array[1] = 0.5 * (sc_alt_arr[1][:,:,:-1,:] +
                                         sc_alt_arr[1][:,:,1:,:])

                # Set land points to NaN and adjust with any scaling
                # Factor offset
                # Note using isnan/sum is relatively fast, but less than 
                # bottleneck external lib
                print 'SC ARRAY MIN MAX : ', np.nanmin(sc_array[0]), np.nanmax(sc_array[0])
                sc_array[0][t_mask == 0] = np.NaN
                print 'SC ARRAY MIN MAX : ', np.nanmin(sc_array[0]), np.nanmax(sc_array[0])
                
                if self.key_vec:
                    sc_array[1][t_mask == 0] = np.NaN
                    if not np.isnan(np.sum(meta_data[vn + 1]['sf'])):
                        sc_array[1] *= meta_data[vn + 1]['sf']
                    if not np.isnan(np.sum(meta_data[vn + 1]['os'])):
                        sc_array[1] += meta_data[vn + 1]['os']

                # Now collapse the extracted data to an array 
                # containing only nearest neighbours to dest bdy points
                # Loop over the depth axis
                for dep in range(sc_z_len):
                    tmp_arr = [None, None]
                    # Consider squeezing
                    tmp_arr[0] = sc_array[0][0,dep,:,:].flatten(1) #[:,:,dep]
                    if not self.key_vec:
                        sc_bdy[vn, dep, :, :] = self._flat_ref(tmp_arr[0], ind)
                    else:
                        tmp_arr[1] = sc_array[1][0,dep,:,:].flatten(1) #[:,:,dep]
                        # Include in the collapse the rotation from the
                        # grid to real zonal direction, ie ij -> e
                        sc_bdy[vn, dep, :] = (tmp_arr[0][ind[:]] * self.gcos - 
                                              tmp_arr[1][ind[:]] * self.gsin)
                        # Include... meridinal direction, ie ij -> n
                        sc_bdy[vn+1, dep, :] = (tmp_arr[1][ind[:]] * self.gcos +
                                                tmp_arr[0][ind[:]] * self.gsin)
                
                # End depths loop
                print ' END DEPTHS LOOP '
            # End Looping over vars
            print ' END VAR LOOP '
            # ! Skip sc_bdy permutation
            
            x = sc_array[0]
            y = np.isnan(x)
            z = np.invert(np.isnan(x))
            x[y] = 0
            print 'nans: ', np.sum(y[:])
            #x = x[np.invert(y)]
            print x.shape, np.sum(x[z], dtype=np.float64), np.amin(x), np.amax(x)

            # Calculate weightings to be used in interpolation from
            # source data to dest bdy pts. Only need do once.
            if self.first:
                # identify valid pts
                data_ind = np.invert(np.isnan(sc_bdy[0,:,:,:]))
                # dist_tot is currently 2D so extend along depth
                # axis to allow single array calc later, also remove
                # any invalid pts using our eldritch data_ind
                print 'DIST TOT ZEROS BEFORE ', np.sum(self.dist_tot == 0)
                self.dist_tot = (np.repeat(self.dist_tot, sc_z_len).reshape(
                            self.dist_tot.shape[0], 
                            self.dist_tot.shape[1], sc_z_len)).transpose(2,0,1)
                self.dist_tot *= data_ind
                print 'DIST TOT ZEROS ', np.sum(self.dist_tot == 0)
                
                print 'DIST IND ZEROS ', np.sum(data_ind == 0)

                # Identify problem pts due to grid discontinuities 
                # using dists >  lat
                over_dist = np.sum(self.dist_tot[:] > 4)
                if over_dist > 0:
                    raise RuntimeError('''Distance between source location
                                          and new boundary points is greater
                                          than 4 degrees of lon/lat''')

                # Calculate guassian weighting with correlation dist
                r0 = self.settings['r0']
                dist_wei = (1 / np.power(r0 * (2 * np.pi), 0.5) * 
                            np.exp(-0.5 * np.power(self.dist_tot / r0, 2)))
                # Calculate sum of weightings
                dist_fac = np.sum(dist_wei * data_ind, 2)
                # identify loc where all sc pts are land
                nan_ind = np.sum(data_ind, 2) == 0
                print 'NAN IND : ', np.sum(nan_ind)
                
                # Calc max zlevel to which data available on sc grid
                data_ind = np.sum(nan_ind == 0, 0) - 1
                # set land val to level 1 otherwise indexing problems
                # may occur- should not affect later results because
                # land is masked in weightings array
                data_ind[data_ind == -1] = 0
                # transform depth levels at each bdy pt to vector
                # index that can be used to speed up calcs
                data_ind += np.arange(0, sc_z_len * self.num_bdy, sc_z_len)

                # ? Attribute only used on first run so clear. 
                del self.dist_tot
                
            # weighted averaged onto new horizontal grid
            for vn in range(self.nvar):
                print ' sc_bdy ', np.nanmin(sc_bdy), np.nanmax(sc_bdy)
                dst_bdy = (np.nansum(sc_bdy[vn,:,:,:] * dist_wei, 2) / 
                           dist_fac)
                print ' dst_bdy ', np.nanmin(dst_bdy), np.nanmax(dst_bdy)
                # Quick check to see we have not got bad values
                if np.sum(dst_bdy == np.inf) > 0:
                    raise RuntimeError('''Bad values found after 
                                          weighted averaging''')
                # weight vector array and rotate onto dest grid
                if self.key_vec:
                    # [:,:,:,vn+1]
                    dst_bdy_2 = (np.nansum(sc_bdy[vn+1,:,:,:] * dist_wei, 2) / 
                                 dist_fac)
                    print 'time to to rot and rep '
                    print np.nanmin(dst_bdy), np.nanmax(dst_bdy)
                    print self.rot_str, 'en to %s' %self.rot_dir, dst_bdy.shape
                    dst_bdy = self.rot_rep(dst_bdy, dst_bdy_2, self.rot_str,
                                      'en to %s' %self.rot_dir, self.dst_gcos, self.dst_gsin)
                    print np.nanmin(dst_bdy), np.nanmax(dst_bdy)
                # Apply 1-2-1 filter along bdy pts using NN ind self.id_121
                if self.first:
                    tmp_valid = np.invert(np.isnan(
                                            dst_bdy.flatten(1)[self.id_121]))
                    # Finished first run operations
                    self.first = False

                dst_bdy = (np.nansum(dst_bdy.flatten(1)[self.id_121] * 
                           self.tmp_filt, 2) / np.sum(self.tmp_filt * 
                           tmp_valid, 2))
                # Set land pts to zero

                print ' pre dst_bdy[nan_ind] ', np.nanmin(dst_bdy), np.nanmax(dst_bdy)
                dst_bdy[nan_ind] = 0
                print ' post dst_bdy ', np.nanmin(dst_bdy), np.nanmax(dst_bdy)
                # Remove any data on dst grid that is in land
                dst_bdy[:,np.isnan(self.bdy_z)] = 0
                print ' 3 dst_bdy ', np.nanmin(dst_bdy), np.nanmax(dst_bdy)

                # If we have depth dimension
                if not self.isslab:
                    # If all else fails fill down using deepest pt
                    dst_bdy = dst_bdy.flatten(1)
                    dst_bdy += ((dst_bdy == 0) * 
                                dst_bdy[data_ind].repeat(sc_z_len)) 
                    # Weighted averaged on new vertical grid
                    dst_bdy = (dst_bdy[self.z_ind[:,0]] * self.z_dist[:,0] + 
                               dst_bdy[self.z_ind[:,1]] * self.z_dist[:,1])
                    data_out = dst_bdy.reshape(self.dst_dep.shape, order='F')

                    # If z-level replace data below bed !!! make stat
                    # of this as could be problematic
                    ind_z = self.bdy_z.repeat(len(self.dst_dep))
                    ind_z = ind_z.reshape(len(self.dst_dep), 
                                          len(self.bdy_z), order='F')
                    ind_z -= self.dst_dep
                    ind_z = ind_z < 0
                    
                    data_out[ind_z] = np.NaN
                else:
                    data_out = dst_bdy
                    data_out[np.isnan(self.bdy_z)] = np.NaN
                entry = self.d_bdy[self.var_nam[vn]][year]
                if entry['data'] == None:
                    # Create entry with singleton 3rd dimension
                    entry['data'] = np.array([data_out])
                else:
                    entry['data'] = np.concatenate((entry['data'], 
                                                   np.array([data_out])))
                entry['date'] = sc_time[f] #count skipped
            
            nc.close()
            if self.key_vec:
                nc_2.close()
        
        # Need stats on fill pts in z and horiz + missing pts...
    # end month
#end year
# End great loop of crawling chaos
            
                
# # # # # # # # # # # # # # # # #                        


    # Allows reference of two equal sized but misshapen arrays
    # equivalent to Matlab alpha(beta(:)) 
    def _flat_ref(self, alpha, beta):
        return alpha.flatten(1)[beta.flatten(1)].reshape(
                                                   beta.shape, order='F')
    
    # THIS FUNCTION MAY BE BROKEN
    def rot_rep(self, pxin, pyin, cd_type, cd_todo, gcos, gsin):
        if cd_todo.lower() in ['en to i', 'ij to e']:
            print 'en to i or ij to e'
            x,y = pxin, pyin
        elif cd_todo.lower() in ['en to j', 'ij to n']:
            print 'en to j or ij to n'
            x,y = pyin, pxin
        else:
            raise SyntaxError('rot_rep cd_todo %s is invalid' %cd_todo)
        # cd_type = superfluous?
        return x * gcos + y * gsin

    # Convert numeric date from source to dest
    def convert_date(self, date):
        val = self.S_cal.num2date(date)
        return self.D_cal.date2num(val)
    
    def convert_date_to_destination_num(self, date):
        return self.D_cal.date2num(date)

    # Translate Calendars
    def cal_trans(self, source, dest, year, month):
        vals = {'gregorian': 365. + isleap(year), 'noleap': 
                365., '360_day': 360.}
        if source not in vals.keys():
            raise ValueError('Unknown source calendar type: %s' %source)
        # Get month length
        if dest == '360_day':
            ed = 30
        else:
            ed = monthrange(year, month)[1]
        # Calculate scale factor
        sf = vals[source] / vals[dest]
        
        return sf, ed

    # BEWARE FORTRAN V C ordering
    # Replicates and tiles an array
    def _trig_reptile(self, trig, size):
        trig = np.transpose(trig, (2, 1, 0)) # Matlab 2 0 1
        return np.tile(trig, (1, size, 1)) # Matlab size 1 1

    def _get_meta_data(self, fname, var, source_dic):
        nc = Dataset(fname, 'r')
        source_att = nc.variables[var].ncattrs()

        for i in range(len(source_att)):
            varid = nc.variables[var]
            if source_att[i] == 'missing_value':
                source_dic['mv'] = varid.missing_value[:]
            elif source_att[i] == 'scale_factor':
                source_dic['sf'] = varid.scale_factor[:]
            elif source_att[i] == 'add_offset':
                source_dic['os'] = varid.add_offset[:]
            elif source_att[i] == '_FillValue':
                source_dic['fv'] = varid._FillValue #[:]
     
        nc.close()
        return source_dic       




















 
