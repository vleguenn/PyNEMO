'''
NcML reading implementation using pyjnius
@author: Mr. Srikanth Nagella
'''
# pylint: disable=E1103
# pylint: disable=no-name-in-module
#Loading of NCML jar file
import os
import string
import logging
import numpy as np
import jnius_config
ncmlpath, file_name = os.path.split(__file__)
ncmlpath = os.path.join(ncmlpath, "jars", "netcdfAll-4.6.jar") 
jnius_config.set_classpath('.',ncmlpath)
try:
    if os.environ['http_proxy'] is not None:
        #split the proxy name and port
        proxylist = string.split(os.environ['http_proxy'],':')
        proxy_host = proxylist[0]
        proxy_port = proxylist[1]        
        jnius_config.add_options('-Dhttp.proxyHost='+proxy_host,'-Dhttp.proxyPort='+proxy_port)
except:
    print "Didn't find a proxy environment variable"
NetcdfDataset = None
NcMLReader = None
Section = None
try:
    from jnius import autoclass
    def init_jnius():
        global NetcdfDataset
        global NcMLReader
        global Section
        NetcdfDataset = autoclass('ucar.nc2.dataset.NetcdfDataset')
        NcMLReader = autoclass('ucar.nc2.ncml.NcMLReader')
        Section = autoclass('ucar.ma2.Section')
    init_jnius()
except ImportError:
    print 'Warning: Please make sure pyjnius is installed and jvm.dll/libjvm.so/libjvm.dylib is in the path'
    
class Data(object):

    def __init__(self, filename, variable):
        self.logger = logging.getLogger(__name__)
        self.variable = variable
        self.file_name = filename
        self.dimensions = self._get_dimensions()
        
    def __str__(self):
        return "PyNEMO Data Object from file: %s and variable %s" % self.file_name, self.variable
    
    def __len__(self):
        """ Returns the length of the variable """
        try:
            dataset = NetcdfDataset.openFile(self.file_name, None)
            dvar = dataset.findVariable(self.variable)
            if dvar is None:
                dataset.close()
                raise KeyError()
            val  = dvar.getDimension(0).getLength()
            dataset.close()
            return val
        except KeyError:
            self.logger.error('Cannot find the requested variable '+self.variable)
        except (IOError, RuntimeError):
            self.logger.error('Cannot open the file '+self.file_name)
        return None
            
    def __getitem__(self, val):
        """ Returns the data requested """
        if type(val) != tuple:
            val = (val,)
        try:
            dataset = NetcdfDataset.openFile(self.file_name, None)
            dvar = dataset.findVariable(self.variable)
            if dvar is None:
                dataset.close()
                raise KeyError()
            dims = dvar.getShape()
            # get the requested slice and extract that information from jnius
            # Check if the request data is with in the dataset dimensions
            start = [0]*len(dims)
            stop = dims
            stride = [1]*len(dims)
            new_dims = tuple()
            np_input = False
            for idx in range(0,len(dims)):
                try:
                    if val[idx].start is not None:
                        start[idx] = val[idx].start
                    if val[idx].step is not None:
                        stride[idx] = val[idx].step
                    if val[idx].stop is not None:
                        if val[idx].stop == -1:
                            val[idx] = stop[idx] - 1
                        elif val[idx].stop > stop[idx]:
                            val[idx].stop = stop[idx]
                        stop[idx] = (val[idx].stop - start[idx])//stride[idx]
                        if (val[idx].stop - start[idx])%stride[idx] != 0:
                            stop[idx] = stop[idx] + 1
                    new_dims = new_dims+(stop[idx],)
                except IndexError:
                    pass
                except AttributeError:
                    if isinstance(val[idx],int):
                        start[idx] = val[idx]
                        stop[idx] = 1
                    elif isinstance(val[idx],np.ndarray):
                        new_dims = new_dims+(val[idx].shape)
                        np_input = True
            # Create a section object that represents the requested slice 
            start = [int(i) for i in start]
            stop = [int(i) for i in stop]
            stride = [int(i) for i in stride]
            selected_section = Section(start,stop,stride)
            data_array = dvar.read(selected_section)
            retval = data_array.copyToNDJavaArray() 
            #TODO: copy it into numpy instead of Java array and then convert to numpy
            # convert to numpy array
            retval = np.asarray(retval)
            self.logger.info(retval.shape)
            if np_input: #if an array is passed as selection
                ret_dim_list = ()
                for idx in range(0,len(dims)):
                    if isinstance(val[idx], np.ndarray): 
                        ret_dim_list2 = ret_dim_list+(val[idx],)
                        # can't do all the reductions at once due to Index Error: shape mismatch
                        retval = retval[ret_dim_list2]  
                    ret_dim_list = ret_dim_list+(slice(None,None,None),)
                self.logger.info(ret_dim_list)                        
                self.logger.info(retval.shape)
                self.logger.info(ret_dim_list)
            # reshape to reflect the request
            retval = np.reshape(retval, new_dims)
            dataset.close()
            return retval
        except KeyError:
            self.logger.error('Cannot find the requested variable '+self.variable)
        except (IOError, RuntimeError):
            self.logger.error('Cannot open the file '+self.file_name)
        return None
         
    def _get_dimensions(self):
        """ Returns the dimensions of the variables """
        try:
            dataset = NetcdfDataset.openFile(self.file_name, None)
            dvar = dataset.findVariable(self.variable)
            if dvar is None:
                dataset.close()
                raise KeyError()
            retval = tuple(dvar.getDimensionsString().split(' '))
            dataset.close()
            return retval
        except KeyError:
            self.logger.error('Cannot find the requested variable '+self.variable)
        except (IOError, RuntimeError):
            self.logger.error('Cannot open the file '+self.file_name)
        return None     
    
    def get_attribute_value(self, attr_name):
        """ Returns the attribute value of the variable """
        try:
            dataset = NetcdfDataset.openFile(self.file_name, None)
            dvar = dataset.findVariable(self.variable)
            if dvar is None:
                dataset.close()
                raise KeyError()
            attr = dvar.findAttributeIgnoreCase(attr_name)
            if attr is not None:
                retval = attr.getValue(0)            
            dataset.close()
            return retval
        except KeyError:
            self.logger.error('Cannot find the requested variable '+self.variable)
        except (IOError, RuntimeError):
            self.logger.error('Cannot open the file '+self.file_name)
        return None  


class SourceData(object):
    logger = logging.getLogger(__name__)
    def __init__(self):
        """ This class is the source data that holds the dataset information """
        self.file_name = None
        self.units = None
        self.calendar = None
        self.date = None
        self.seconds = None

    def __getitem__(self, val):
        """ Returns the data requested """
        try:
            return Data(self.file_name, val)
        except KeyError:
            self.logger.error('Cannot find the requested variable '+self.variable)
        except (IOError, RuntimeError):
            self.logger.error('Cannot open the file '+self.file_name)
        return None
       
    def get_meta_data(self, variable, source_dic):
        """ Returns a dictionary with meta data information correspoinding to the variable """
        #source_dic = {}
        try:
            dataset = NetcdfDataset.openFile(self.file_name, None)
            dvar = dataset.findVariable(self.variable)
            if dvar is None:
                dataset.close()
                raise KeyError()
            source_dic['sf'] = 1
            source_dic['os'] = 0
            mv_attr = dvar.findAttributeIgnoreCase('missing_value')
            if mv_attr is not None:
                source_dic['mv'] = mv_attr.getValues().copyToNDJavaArray()
            sf_attr = dvar.findAttributeIgnoreCase('scale_factor')
            if sf_attr is not None:
                source_dic['sf'] = sf_attr.getValues().copyToNDJavaArray()
            os_attr = dvar.findAttributeIgnoreCase('add_offset')
            if os_attr is not None:
                source_dic['os'] = os_attr.getValues().copyToNDJavaArray()
            fv_attr = dvar.findAttributeIgnoreCase('_FillValue')
            if fv_attr is not None:
                source_dic['mv'] = fv_attr.getValue(0)
            dataset.close()
            return source_dic            
        except KeyError:
            self.logger.error('Cannot find the requested variable '+variable)
        except (IOError, RuntimeError):
            self.logger.error('Cannot open the file '+self.file_name)
        return None        
            
    def close(self):
        """ This is not yet implemented. TODO: keep the netcdf file open until its expicitly 
        closed """
        pass   
    
    
from os import listdir
from netCDF4 import netcdftime
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
            varid = Data(filename, 'time_counter')
            src_data = SourceData()
            src_data.file_name = filename
            src_data.units = varid.get_attribute_value('units')
            src_data.calendar = varid.get_attribute_value('calendar')
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
        
      
def GetFile(file_path):
    """ This method looks at the file_path and decides on appropriate file reader """
    if os.path.isfile(file_path) and file_path.endswith(".nc"):
        data = SourceData()
        data.file_name = file_path
        return data
    elif url_exists(file_path) and file_path.endswith(".nc"):
        data = SourceData()
        data.file_name = file_path
        return data
    else:
        raise NotImplementedError, "The file reader for the requested file is not supported"
     
import urllib2
def url_exists(location):
    request = urllib2.Request(location)
    request.get_method = lambda : 'HEAD'
    try:
        response = urllib2.urlopen(request)
        return True
    except urllib2.HTTPError:
        return False