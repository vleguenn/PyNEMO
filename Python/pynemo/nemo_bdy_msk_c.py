#################################################################
#                    Nemo Boundary Mask                         #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                               #
#  This class is initialised with a netcdf file location, the   #
#  bathymetry data is interpreted into a boundary mask array.   #
#  Optional arguments modify the mask by removing regions. The  #
#  'custom_areas' argument allows the removal of any additional #
#  generic areas of the grid, and is specified as nested lists  #
#  of coordinates that match x/y blocks to remove- the same way #
#  as the predefined regions.                                   #
#                                                               #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#  Written by John Kazimierz Farey, August 2012                 #
#  Port of Matlab code by James Harle                           #
#################################################################

# Attributes: bdy_msk

# reduced option not implemented

from netCDF4 import Dataset
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import logging

class Mask:

    def __init__(self, ncfname, med=True, blk=True, hud=True, bal=True,
                 v2=True, custom_areas=None):
        #Logging
        self.logger = logging.getLogger(__name__)
        # Debug
        #self.logger.info('Mask init, ' 'med', med, 'blk', blk, 'hud', hud, 'bal ',bal, 'v2 ',v2)
        self.logger.info( '''Mask initialising with remove med: %s, blk: %s, hud: %s,
                 bal: %s, v2: %s''', med, blk, hud, bal, v2)

        # Read from default location
        if ncfname == 0:
            ncfname = '/Users/jdha/Projects/pynemo_data/NEMO_bdy_tools/grid_C/NNA_R12_bathy_meter_bench.nc'

        # Open and read bathymetry netcdf file
        nc = Dataset(ncfname, 'r')
        bdy_msk = np.asarray(nc.variables['Bathymetry'][:,:])

        # Normalise all elements > 0 to 1, otherwise to 0
        bdy_msk = np.around((bdy_msk + .5).clip(0,1))      
        bdy_msk[bdy_msk ==1] = -1          
        tmp = bdy_msk[3:398,3:348]
        tmp[tmp == -1] = 1
        bdy_msk[3:398,3:348] = tmp
        
        tmp = bdy_msk[0:59,280:350]
        tmp[tmp == 1] = -1
        bdy_msk[0:59,280:350] = tmp
        
        # Debug line
        self.logger.debug ('bdy_msk FINAL: Total: %s , | 1s: %s | -1s: %s | 0s: %s',  
                          np.sum(bdy_msk[:,:]),np.sum(bdy_msk[:]==1),np.sum(bdy_msk[:]==-1),np.sum(bdy_msk[:] == 0))

        # Set our attribute
        self.bdy_msk = bdy_msk
        # Close NetCDF file
        nc.close()
            

# # # # # # #
# Functions #
# # # # # # #

    # Adds or removes an area from grid
    # bathy = mask
    # mod = list of lists of region coords, formatted as the constants
    # mode = defaults to remove region (1s -> -1s). mode='add' reverses
    def _mod_area(self, bathy, mod, mode=''):
        for co in range(len(mod)):
            area = bathy[mod[co][0]:mod[co][1], mod[co][2]:mod[co][3]]
            if mode == 'add':
                area = np.abs(area)
                bathy[mod[co][0]:mod[co][1], mod[co][2]:mod[co][3]] = area
            else:
                area = area[:] == 1
                bathy[mod[co][0]:mod[co][1], mod[co][2]:mod[co][3]][area] = -1
        
        return bathy

    # Removes regions of mask as specified
    def _rmv_regions(self, bathy, med, blk, hud, bal, v2):
        if med:
            bathy = self._mod_area(bathy, self._MED)
        if blk:
            bathy = self._mod_area(bathy, self._BLK)
        if hud:
            bathy = self._mod_area(bathy, self._HUD)
        if bal:
            bathy = self._mod_area(bathy, self._BAL)
            if v2:
                bathy = self._v2_norwegian(bathy)

        return bathy


    # Adds a bit of norwegian back in
    def _v2_norwegian(self, bathy):
        bathy = self._mod_area(bathy, [[479,610, 1159,1227]], mode='add')
        bathy = self._mod_area(bathy, [[494,503, 1225,1227]])
        bathy[535,1226] = -1

        return bathy

    # If baltic and hudson bay removed
    def _rmv_bal_hud(self, bathy):
        bathy = self._mod_area(bathy, [[763,None, None,None]])
        bathy = self._mod_area(bathy, [[763,783, 379,420]], 'add')
        
        # add back a bit of norwegian
        bathy = self._mod_area(bathy, [[763,849, 749,1200]], 'add')
        bathy[848,1149:1151] = -1

        return bathy


    def _plot_map(self,area):
        self.logger.info('Processing Surface Plot of the data ...')
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        Y = np.arange(1, area.shape[0]+1, 1)
        X = np.arange(1, area.shape[1]+1, 1)
        X, Y = np.meshgrid(X, Y)
        Z = area
        self.logger.debug ('Data shape: %s x %s',area.shape[0],area.shape[1])
        surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0, antialiased=False)

        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

        fig.colorbar(surf, shrink=0.5, aspect=5)

        plt.show()
        return
# Below code for plotting projections not currently in use
"""
    # Map plotting unfinised
    def _plot_map(self):

        plotter.pcolor(bdy_msk[500:931,1000:1500])
        plotter.show()
        plotter.savefig('/login2/jofa/first/fig1.png')



            # ??? get boundary REDUNDANt? matlab does it
            #tmp = np.asarray(bdy_msk[930,1312:1320])

        ##### PLOT unfinished
        '''
        from pyproj import Proj
        from pylab import plot, show

        zone = int(round((bdy_msk[:,1].mean() +180) / 6.))

        p = Proj(proj='utm', zone=zone, ellps='WGS84')
        xs,ys = p(bdy_msk[:,:], bdy_msk[:,:])

        lines = plot(xs, ys, 'r.-')
        show()

        #wolf = plot(p(bdy_msk[:,:], 'r.-'))
        #show()
        '''
        #import matplotlib
        #matplotlib.rcParams['backend'] = 'Qt4Agg'

        #import matplotlib.pyplot as plotter

        #plotter.plot(bdy_msk[0:200,0:200], 'r+')

"""
