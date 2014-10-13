'''
Created on 3 Oct 2014

@author: Mr. Srikanth Nagella
Netcdf writer for the bdy output
'''
from netCDF4 import Dataset
import numpy as np
def WriteDataToFile(filename, variableName, data):
    ncid = Dataset(filename, 'a', clobber=False, format='NETCDF4')
    count = data.shape
    
    threeDimensionVariables = ['votemper', 'vosaline', 'N1p', 'N3n', 'N5s']
    twoDimensionVariables = ['sossheig','vobtcrtx', 'vobtcrty', 'iicethic', 'ileadfra', 'isnowthi']
    
    if variableName in threeDimensionVariables:
        if len(count) == 3:
            count += (1L,)
        ncid.variables[variableName][:,:,:,:] = np.reshape(data, count)[:,:,:,:]
    elif variableName in twoDimensionVariables:
        if len(count) == 2:
            count += (1L,)
        elif len(count) == 1:
            count += (1L,1L,)
        ncid.variables[variableName][:,:,:] = np.reshape(data, count)[:,:,:]
    elif variableName == 'time_counter':
        ncid.variables[variableName][:] = data[:]
    else:
        if len(count) == 1:
            ncid.variables[variableName][:] = data[:]
        elif len(count) == 2:
            ncid.variables[variableName][:,:] = data[:,:]
        elif len(count) == 3:
            ncid.variables[variableName][:,:,:] = data[:,:,:]
                                         
        
    ncid.close()
        
        
        