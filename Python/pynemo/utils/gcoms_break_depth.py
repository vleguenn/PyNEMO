'''
Rewritting the break depth implementation from matlab version

@author: Mr. Srikanth Nagella
'''
# pylint: disable=E1103
# pylint: disable=no-name-in-module
import numpy as np
import math
#import matplotlib.pyplot as plt
import scipy.ndimage as ndimage
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
    return lb
