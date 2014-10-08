'''
Created on 3 Oct 2014

@author: Mr. Srikanth Nagella
Netcdf writer for the bdy output
'''
from netCDF4 import Dataset
import numpy as np
def WriteDataToFile(filename, variableName, data):
    ncid = Dataset(filename, 'w', clobber=True, format='NETCDF4')
    count = data.shape
    
    threeDimensionVariables = ['votemper', 'vosaline', 'N1p', 'N3n', 'N5s']
    twoDimensionVariables = ['sossheig','vobtcrtx', 'vobtcrty', 'iicethic', 'ileadfra', 'isnowthi']
    
    if variableName in threeDimensionVariables:
        if count.size == 3:
            count = [count, 1]
        ncid.variables[variableName] = np.reshape(data, count)
    elif variableName in twoDimensionVariables:
        if count.size == 2:
            count = [count, 1]
        ncid.variables[variableName] = np.reshape(data, count)
    elif variableName == 'time_counter':
        ncid.variables[variableName] = data
    else:
        ncid.variables[variableName] = data
        
        
        