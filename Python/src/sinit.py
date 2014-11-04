# first init

from netCDF4 import Dataset
import numpy
import copy
import pickle

use_file = False


try:
	f = open('F:/NEMO_bdy_tools/first/grid.pkl', 'r')
	bdy_msk = pickle.load(f)
	use_file = True

except IOError:


	nc = Dataset('F:/NEMO_bdy_tools/bdy_matlab/grid_C/NNA_R12_bathy_meter_bench.nc', 'r', format='NetCDF3_CLASSIC')

	bdy_msk = numpy.asarray(nc.variables['Bathymetry'][:,:])

# since its opened up for reading only, this is redundant
# copy data
#bdy_msk = copy.deepcopy(bathy)


#def debug():
#	print 'DEBUG: bdy_msk 90 0-150' 
#	print bdy_msk[90,0:150]
#	print '\n ' 
#
#debug()

# Go thru values and make x > 0 == 1
	xl = range(len(bdy_msk))
	yl = range(len(bdy_msk[0]))

# Find quicker way later
	for x in xl:
		for y in yl:
			if bdy_msk[x][y] > 0:
				bdy_msk[x][y] = 1

#debug()

#create boundary  1s to -1s	
	for x in range(1312,1320):
		if bdy_msk[930][x] == 1:
			bdy_msk[930][x] = -1

	for x in range(len(bdy_msk[0])):
		if bdy_msk[0][x] == 1:
			bdy_msk[0][x] = -1



#Now lets pickle the data to save reprocessing next time
	f = open('F:/NEMO_bdy_tools/first/grid.pkl', 'w')
	pickle.dump(bdy_msk, f)



# get boundary REDUNDANt? matlab does it
#tmp = numpy.asarray(bdy_msk[930,1312:1320])

#print 'DEBUG: tmp \n ', tmp, '\n'

# set boundary 1s to -1s
#for x in range(7):
#	if tmp[0][x] == 1:
#		tmp[0][x] = -1
#
#print 'DEBUG: tmp \n', tmp, '\n'






##### PLOT
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

import matplotlib.pyplot as plotter

#plotter.plot(bdy_msk[0:200,0:200], 'r+')
plotter.pcolor(bdy_msk[0:930,0:1600])



plotter.show()
plotter.savefig('F:/NEMO_bdy_tools/first/fig1svr.png')






if use_file == False:
	nc.close()
	
