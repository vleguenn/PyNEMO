'''
Rewritting the break depth implementation from matlab version

@author: Mr. Srikanth Nagella
'''
# pylint: disable=E1103
# pylint: disable=no-name-in-module
import numpy as np
import math
import copy
import pyproj

#import matplotlib.pyplot as plt
import scipy.ndimage as ndimage
from seawater.extras import dist
def gcoms_break_depth(bathy):
    """ This creates a mask for the break depth using histograms """
    ocean_depth = bathy[...]
    ocean_depth = ocean_depth[ocean_depth > 0]

    depth_bin = 10.0
    depth_max = np.max(ocean_depth)
    num_bin = int(math.floor(depth_max/depth_bin))

    # Compute the histogram of depth values over the whole Domain
    depth_vec = (np.arange(1,num_bin+1)+0.5)*depth_bin    
    histbat, dummy = np.histogram(ocean_depth, num_bin)
    max_hist_value_index = np.argmax(histbat)
    #z_smo = (max_hist_value_index * depth_bin)/2.0
    z_smo = 100.0
    nsmo = math.floor(z_smo/depth_bin)
#    print nsmo, z_smo, histbat.dtype
    #print max_hist_value_index
#    plt.subplot(211)
#    plt.hist(ocean_depth, num_bin) 
    #plt.show()
    hist_smooth = ndimage.uniform_filter(histbat.astype(float), int(nsmo)*2+1, mode='nearest')
#    print histbat.shape, dummy.shape
#    plt.subplot(212)
#    plt.bar(dummy[:-1],hist_smooth)
    
#    plt.show()
#    print histbat
#    print hist_smooth
    kshelf = -1
    kbreak = -1
    kplain = -1
    kfloor = -1
    histfloor = 0.0
    for depth_bin_index in range(0,num_bin-1):
        if kshelf == -1:
            if hist_smooth[depth_bin_index] > hist_smooth[depth_bin_index+1]:
                kshelf = depth_bin_index;
        elif kbreak == -1:
            if hist_smooth[depth_bin_index] < hist_smooth[depth_bin_index+1]:
                kbreak = depth_bin_index
        elif kplain == -1:
            if hist_smooth[depth_bin_index] > hist_smooth[depth_bin_index+1]:
                kplain = depth_bin_index
                histfloor = hist_smooth[depth_bin_index]
    
    depth_shelf = depth_vec[kshelf]
    depth_break = depth_vec[kbreak]
    depth_plain = depth_vec[kplain]
#    print kshelf,kbreak,kplain
#    print 'Approximate depths: shelf=%sm, break=%sm, plain=%sm' % (depth_shelf,depth_break,depth_plain)
    h_max = math.floor(depth_break/100)*100
    return depth_shelf, h_max
    
def gcoms_boundary_masks(bathy,ov,lv):
    tmp = np.pad(bathy, (1, 1), 'constant', constant_values=(np.NaN, np.NaN))
    tmp[tmp==ov] = np.NaN
    
    tmp1 = tmp[1:-1, :-2] + tmp[1:-1, 2:] + tmp[:-2, 1:-1] + tmp[2:, 1:-1]

    ob = np.logical_and(np.logical_and(np.isnan(tmp1), bathy != ov) , bathy != lv)
    
    tmp = np.pad(bathy, (1, 1), 'constant', constant_values=(-1,-1))
    tmp[tmp==lv] = np.NaN
    
    tmp1 = tmp[1:-1, :-2] + tmp[1:-1, 2:] + tmp[:-2, 1:-1] + tmp[2:, 1:-1]

    lb = np.logical_and(np.logical_and(np.isnan(tmp1), bathy!=ov), bathy!=lv)
    return ob, lb


def polcoms_select_domain(bathy, lat, lon, roi):
    """ This calculates the shelf break
    :param bathy: This is the input bathymetry data
    :param lat: Latittude array
    :param lon: Longitude array  
    :type bathy: numpy array
    :type lat: numpy array
    :type lon: numpy array
    :return: returns the depth_shelf, h_max
    :rtype: numpy arrays
    
    :Example:
    """
    dr = 200
    dy = 0.1
    dx = 0.1
    
    #create a copy of bathy
    bathy_copy = bathy.copy()
#    bathy[bathy>=0] = 0;
#    bathy = bathy*-1
    global_ind = bathy_copy*np.NaN
    r = np.ceil(dr/(np.pi/180*6400)/dy)
    
    tmp = bathy_copy[roi[2]-r:roi[3]+r,roi[0]-r:roi[1]+r]
    lat = lat[roi[2]-r:roi[3]+r,roi[0]-r:roi[1]+r]
    lon = lon[roi[2]-r:roi[3]+r,roi[0]-r:roi[1]+r]
    
    nanind = np.isnan(tmp) 
    tmp[nanind] = -1
    dummy, lb = gcoms_boundary_masks(tmp, -1,0)
    Zshelf, Hmax = gcoms_break_depth(tmp)
    tmp[tmp>Zshelf] = -1
    tmp[np.logical_and(np.logical_and(tmp!=0, np.logical_not(np.isnan(tmp))), tmp!=-1)] = 1
    
    ob, dummy = gcoms_boundary_masks(tmp, -1, 0)
    
    lat_ob = np.ravel(lat,order='F')[np.ravel(ob,order='F')]
    lon_ob = np.ravel(lon,order='F')[np.ravel(ob,order='F')]
    
    len_lat = len(lat)
    len_lon = len(lon[0,:])
    lat_lon_index = np.nonzero( np.logical_and(lat == lat_ob[0], lon == lon_ob[0]))    
    for idx in range(0, len(lat_ob)):        
        lat_lon_index = np.nonzero( np.logical_and(lat == lat_ob[idx], lon == lon_ob[idx]))
        lat_slice = slice(max(lat_lon_index[0]-r,0),min(lat_lon_index[0]+r+1,len_lat))
        lon_slice = slice(max(lat_lon_index[1]-r,0),min(lat_lon_index[1]+r+1,len_lon))   
        lat_pts = lat[lat_slice, lon_slice]
        lon_pts = lon[lat_slice, lon_slice]
        lat_pts_shape = lat_pts.shape
        lat_pts = np.ravel(lat_pts)
        lon_pts = np.ravel(lon_pts)
        # NOTE: seawater package calculates the distance from point to the next point in the array
        # that is the reason to insert reference point before every point
        #lat_pts = np.insert(lat_pts,range(0,len(lat_pts)), lat_ob[idx])
        #lon_pts = np.insert(lon_pts,range(0,len(lon_pts)), lon_ob[idx])
        #distance_pts = seawater.dist(lat_pts, lon_pts)
        #distances repeat themselves so only pick every alternative distance
        #distance_pts = distance_pts[0][::2]
        
        #Using pyproj
        geod = pyproj.Geod(ellps='WGS84')
        dummy,dummy, distance_pts = geod.inv(len(lon_pts)*[lon_ob[idx]],len(lat_pts)*[lat_ob[idx]], lon_pts, lat_pts)
        distance_pts=distance_pts/1000.0
                         
        distance_pts = np.reshape(distance_pts, lat_pts_shape)
        distance_pts[distance_pts>dr] = np.NaN
        distance_pts[np.logical_not(np.isnan(distance_pts))] = 1
        tmp1 = tmp[lat_slice, lon_slice]
        tmp1[np.logical_and(tmp1==-1, distance_pts==1)] = 1
        tmp[lat_slice, lon_slice] = tmp1
        
    lat_lb = lat[lb]
    lon_lb = lon[lb]
    
    for idx in range(0, len(lat_lb)):
        lat_lon_index = np.nonzero( np.logical_and(lat == lat_lb[idx], lon == lon_lb[idx]))
        lat_slice = slice(max(lat_lon_index[0]-r,0),min(lat_lon_index[0]+r,len_lat))
        lon_slice = slice(max(lat_lon_index[1]-r,0),min(lat_lon_index[1]+r,len_lon))   
        lat_pts = lat[lat_slice, lon_slice]
        lon_pts = lon[lat_slice, lon_slice]
        lat_pts_shape = lat_pts.shape        
        lat_pts = np.ravel(lat_pts)
        lon_pts = np.ravel(lon_pts)
        # NOTE: seawater package calculates the distance from point to the next point in the array
        # that is the reason to insert reference point before every point
        #lat_pts = np.insert(lat_pts,range(0,len(lat_pts)), lat_lb[idx])
        #lon_pts = np.insert(lon_pts,range(0,len(lon_pts)), lon_lb[idx])
        #distance_pts = seawater.dist(lat_pts, lon_pts)
        #distances repeat themselves so only pick every alternative distance
        #distance_pts = distance_pts[0][::2]
        
        #Using pyproj
        geod = pyproj.Geod(ellps='WGS84')
        dummy,dummy, distance_pts = geod.inv(len(lon_pts)*[lon_lb[idx]],len(lat_pts)*[lat_lb[idx]], lon_pts, lat_pts) 
        distance_pts=distance_pts/1000.0
        
        distance_pts = np.reshape(distance_pts, lat_pts_shape)
        distance_pts[distance_pts>dr] = np.NaN
        distance_pts[np.logical_not(np.isnan(distance_pts))] = 1
        tmp1 = tmp[lat_slice, lon_slice]
        tmp1[np.logical_and(tmp1==-1, distance_pts==1)] = 1
        tmp[lat_slice, lon_slice] = tmp1        
         
    #Only select largest sub region 
    tmp[nanind] = np.NaN
    ret_val = np.ones(bathy.shape)
    ret_val[roi[2]-r:roi[3]+r,roi[0]-r:roi[1]+r] = tmp
    return ret_val == 1
    #in
         
    

    