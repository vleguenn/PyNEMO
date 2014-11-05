'''
Unit test to nemo_bdy_msk_c

@author: Srikanth Nagella
'''
import unittest
from pynemo.nemo_bdy_msk_c import * 

class Test(unittest.TestCase):


    def testMaskData(self):
        mask = Mask('data/grid_C/NNA_R12_bathy_meter_bench.nc', med=1, blk=1, hud=1, bal=1, v2=1, custom_areas=None)
        print mask.bdy_msk.shape
        self.assertEqual(mask.bdy_msk.shape,(401L,351L),'Mask reading failed')


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()