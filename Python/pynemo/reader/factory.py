'''
This is generic file loader factory.

@author: Mr. Srikanth Nagella
'''

#Global Imports
import os
#Local Imports
from netcdf import GetFile as netcdf_get_file
from netcdf import GetRepository as netcdf_repository
from ncml import GetFile as ncml_get_file 
from ncml import GetRepository as ncml_repository

    
    
import os
def GetRepository(src_dir, grid, t_adjust, reader_type="netcdf"):
    """ Generic method to return the repository either netcdf or ncml based on some
    logic, now passing as reader_type to choose """
    if reader_type == "netcdf":
        return netcdf_repository(src_dir, grid, t_adjust)
    elif reader_type == "ncml":
        return ncml_repository(src_dir, grid, t_adjust)
      
def GetFile(file_path, reader_type="netcdf"):
    """ Generic method to return the file object either netcdf or ncml based on some
    logic, now passing as reader_type to choose"""
    if reader_type == "netcdf":
        return netcdf_get_file(file_path)
    elif reader_type == "ncml":
        return ncml_get_file(file_path)
