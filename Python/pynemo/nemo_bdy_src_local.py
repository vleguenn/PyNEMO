'''
This is a interface for the getting the files.

@author: Mr. Srikanth Nagella
'''
# pylint: disable=E1103
# pylint: disable=no-name-in-module
#External imports
from os import listdir
from netCDF4 import Dataset
from netCDF4 import netcdftime
import logging
#Local Imports

class SourceData(object):
    logger = logging.getLogger(__name__)
    def __init__(self):
        """ This class is the source data that holds the dataset information """
        self.file_name = None
        self.units = None
        self.calendar = None
        self.date = None
        self.seconds = None

    def get_data(self, variable):
        """ Returns the requested data corresponding to the variable """
        try:
            dataset = Dataset(self.file_name, 'r')
            return dataset[variable]
        except KeyError:
            self.logger.error('Cannot find the requested variable '+variable)
        except (IOError, RuntimeError):
            self.logger.error('Cannot open the file '+self.file_name)
        return None
    
    def get_meta_data(self, variable):
        """ Returns a dictionary with meta data information correspoinding to the variable """
        source_dic = {}
        try:
            dataset = Dataset(self.file_name, 'r')
            source_att = dataset.variables[variable].ncattrs()
            varid = dataset.variables[variable]
            source_dic['sf'] = 1
            source_dic['os'] = 0
            for i in range(len(source_att)):
                if source_att[i] == 'missing_value':
                    source_dic['mv'] = varid.missing_value[:]
                elif source_att[i] == 'scale_factor':
                    source_dic['sf'] = varid.scale_factor[:]
                elif source_att[i] == 'add_offset':
                    source_dic['os'] = varid.add_offset[:]
                elif source_att[i] == '_FillValue':
                    source_dic['fv'] = varid._FillValue #[:]
            dataset.close()
            return source_dic            
        except KeyError:
            self.logger.error('Cannot find the requested variable '+variable)
        except (IOError, RuntimeError):
            self.logger.error('Cannot open the file '+self.file_name)
        return None        
            
           
class LocalRepository(object):
    '''
    This manages the files in the local directory
    '''

    def __init__(self, directory, grid_type_list, time_adjust):
        '''
        This takes in directory path as input and returns the required information to the
        bdy.
        Keyword arguments:
        directory -- The directory in which to look for the files
        grid_type_list -- The list of grid types to collect the information
        time_adjust -- amount of time to be adjusted to the time read from file.
        '''
        self.directory = directory
        self.day_interval = 1
        self.hr_interval = 0
        self.grid_source_data = {}
        for grid_type in grid_type_list:
            self.grid_source_data[grid_type] = self._get_source_timedata(grid_type, time_adjust)
        if grid_type_list is not None and len(self.grid_source_data) != 0:
            self._calculate_time_interval()

    def _get_dir_list(self, grid):
        """
        This method scans the directory for a input grid related netcdf files. i.e ending with
        the grid name.
        Keyword arguments:
        grid -- grid name eg. 't','v','u','i'
        """
        fend = '%s.nc' %grid.upper()
        dir_list = listdir(self.directory)
        for i in range(len(dir_list)):
            if dir_list[i][-4:] != fend:
                dir_list[i] = ''
            else:
                dir_list[i] = self.directory + dir_list[i]

        dir_list.sort()
        return filter(None, dir_list)

    def _delta_time_interval(self, time1, time2):
        """ Get the difference between the two times in days and hours"""
        timedif = time2 - time1
        days = timedif / (60 * 60 * 24)
        hrs = timedif % (60 * 60 * 24)
        hrs = hrs / (60 * 60)
        return days, hrs
    
    def _get_source_timedata(self, grid, t_adjust):
        """ Get the source time data information. builds up sourcedata objects of a given grid """
        dir_list = self._get_dir_list(grid)
        src_time = []
        for filename in dir_list:            
            nc = Dataset(filename, 'r')
            varid = nc.variables['time_counter']
            src_data = SourceData()
            src_data.file_name = filename
            src_data.units = varid.units
            src_data.calendar = varid.calendar
            src_data.seconds = varid[0]+t_adjust
            src_data.date = netcdftime.num2date(varid[0]+t_adjust, src_data.units,
                                                src_data.calendar)
            src_time.append(src_data)
        src_time = sorted(src_time, key=lambda x: x.seconds) # sort the list based on time
        return src_time

    def _calculate_time_interval(self):
        """ This method will calculate the time interval of the each grid. If all the grids
        get the same interval then it sets it to the days and hours otherwise it throws an
        error"""
        days = set()
        hrs = set()
        for grid_type in self.grid_source_data.keys():
            day, hr = self._delta_time_interval(self.grid_source_data[grid_type][0].seconds,
                                                self.grid_source_data[grid_type][1].seconds)
            days.add(day)
            hrs.add(hr)
        if len(days) != 1 or len(hrs) != 1:
            raise Exception('All the Grid time interval is not same')
        self.day_interval = list(days)[0]
        self.hr_interval = list(hrs)[0] 


from thredds_crawler.crawl import Crawl
class OpenDAPRepository(LocalRepository):
    def _get_dir_list(self, grid):
        """
        This method scans the opendap server directory for a input grid related netcdf files. i.e ending with
        the grid name.
        Keyword arguments:
        grid -- grid name eg. 't','v','u','i'
        """
        fend = '%s.nc' %grid.upper()
        c = Crawl(self.directory,select=[".*"+fend])
        urls = [s.get("url") for d in c.datasets for s in d.services if s.get("service").lower() == "opendap"]
        urls.sort()
        return urls
    
import os
def GetRepository(src_dir, grid, t_adjust):
    """ This method looks at teh src_dir and decides whether its a local directory or
    a OpenDAP repository and creates appropriate object"""
    if os.path.isdir(src_dir):
        return LocalRepository(src_dir, grid, t_adjust)
    elif url_exists(src_dir):
        return OpenDAPRepository(src_dir, grid, t_adjust)
    else:
        print 'Invalid Data Location'
        return None
        
        
import urllib2
def url_exists(location):
    request = urllib2.Request(location)
    request.get_method = lambda : 'HEAD'
    try:
        response = urllib2.urlopen(request)
        return True
    except urllib2.HTTPError:
        return False