'''
Unit tests for Ncml data reading module using pyjnius

@author: Mr. Srikanth Nagella
'''
import unittest

from pynemo.nemo_bdy_src_ncml import init_jnius
from pynemo.nemo_bdy_src_ncml import Data

import os
class Test(unittest.TestCase):


    def testDataInit(self):
        init_jnius()
        testpath, file_name = os.path.split(__file__)
        testfile = os.path.join(testpath, "test.ncml")
        dataset = Data(testfile,"votemper")
        data = dataset[0,0,0,0:10]
        self.assertEquals(data.shape[0],10,"Extracted requested dimension of data")
        self.assertAlmostEquals(data[0],18.945175,6)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()