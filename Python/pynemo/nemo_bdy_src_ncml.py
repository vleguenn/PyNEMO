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
from jnius import autoclass
def init_jnius():
    global NetcdfDataset
    global NcMLReader
    global Section
    NetcdfDataset = autoclass('ucar.nc2.dataset.NetcdfDataset')
    NcMLReader = autoclass('ucar.nc2.ncml.NcMLReader')
    Section = autoclass('ucar.ma2.Section')

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
            for idx in range(0,len(dims)):
                if isinstance(val[idx],int):
                    start[idx] = val[idx]
                    stop[idx] = 1
                    continue
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
 #               except AttributeError:
 #                   start[idx] = val[idx]
            # Create a section object that represents the requested slice 
            start = [int(i) for i in start]
            stop = [int(i) for i in stop]
            stride = [int(i) for i in stride]
            self.logger.info(start)
            self.logger.info(stop)
            self.logger.info(stride)
            self.logger.info(new_dims)
            selected_section = Section(start,stop,stride)
            data_array = dvar.read(selected_section)
            retval = data_array.copyToNDJavaArray() 
            #TODO: copy it into numpy instead of Java array and then convert to numpy
            # convert to numpy array
            retval = np.asarray(retval)
            self.logger.info(retval.shape)
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