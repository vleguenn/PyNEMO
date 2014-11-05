'''
Unit test for nemo_bdy_gen_c. Boundary class

@author: Mr. Srikanth Nagella
'''
import unittest
from pynemo.nemo_bdy_gen_c import *
from pynemo.nemo_bdy_msk_c import *
class Test(unittest.TestCase):

    settings={}
    def setUp(self):
        self.settings['rimwidth'] = 9
        self.Mask = Mask('data/grid_C/NNA_R12_bathy_meter_bench.nc', med=1, blk=1, hud=1, bal=1, v2=1, custom_areas=None)
        

    def tearDown(self):
        self.settings={}


    def testInvalidGridType(self):
       self.assertRaises(ValueError,Boundary,self.Mask.bdy_msk,self.settings,'x')
       
    def testGridTypeT(self):
        Grid_T = Boundary(self.Mask.bdy_msk,self.settings,'t')
        self.assertEqual(Grid_T.grid_type,'t','Grid type is not T')
        


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()