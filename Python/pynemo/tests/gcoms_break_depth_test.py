'''
This is unit test for gcoms_break_depth

@author: Mr. Srikanth Nagella
'''
import unittest
from pynemo.utils import gcoms_break_depth
from netCDF4 import Dataset
import numpy as np
class Test(unittest.TestCase):


    def setUp(self):
        self.nc = Dataset('data/gebco_1.nc')
        self.bathy = self.nc.variables['topo']
        r = 18
        self.bathy = self.bathy[2:10801:6, 2:21600:6]
        self.bathy = self.bathy[991-r:1295+r,1556-r:1801+r]
        self.bathy[self.bathy>=0]=0
        self.bathy = -1*self.bathy
        self.bathy[np.isnan(self.bathy)] = -1
        
        #gcoms_break_depth.gcoms_boundary_masks(self.bathy, -1, 0)
        
        print self.bathy.shape
        print self.bathy
        pass


    def tearDown(self):
        self.nc.close()
        pass


    def testGcomsBreakDepth(self):
        gcoms_break_depth.gcoms_break_depth(self.bathy)
        pass


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()