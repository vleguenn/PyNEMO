'''
Testing the nemo_bdy_src_local

@author: Mr. Srikanth Nagella
'''
import unittest

from pynemo.nemo_bdy_src_local import LocalRepository
from pynemo.nemo_bdy_src_local import OpenDAPRepository
from netcdftime import utime
class Test(unittest.TestCase):

    def testGetTimeInterval(self):
        repo = LocalRepository('../../data/srcdata_low_res_C',[],0)
        t0 = 6.88392E8
        t1 = 6.88824E8
        days, hrs = repo._delta_time_interval(t0, t1)
        self.assertEquals(days, 5, 'return a valid time difference in days')
        self.assertEquals(hrs, 0, 'return a valid time difference in hours')

    def testGetDirList(self):
        repo = LocalRepository('data/srcdata_low_res_C/',['t'],0)
        file_list = repo._get_dir_list('t')
        print file_list
        self.assertGreater(len(file_list),0,'No files found in the directory')    
    
    def testGetMetaData(self):
        repo = LocalRepository('data/srcdata_low_res_C/',['t'],0)
        self.assertGreater(repo.grid_source_data['t'][0].get_meta_data('votemper')['fv'], 9.969201E36, 'Fill Value is not read properly')
        
    def testRemoteGetDirList(self):
        repo = OpenDAPRepository('http://esurgeod.noc.soton.ac.uk:8080/thredds/catalog/PyNEMO/data/catalog.xml',['t'],0)
        file_list = repo._get_dir_list('t')
        self.assertGreater(len(file_list),0,'No files found in the directory')
        
    def testRemoteGetTimeInterval(self):
        repo = OpenDAPRepository('http://esurgeod.noc.soton.ac.uk:8080/thredds/catalog/PyNEMO/data/catalog.xml',['t'],0)
        self.assertEquals(repo.day_interval,5,'Not a 5 day interval repository')
                
    def testInit(self):
        repo = LocalRepository('data/srcdata_low_res_C/',['t'],0)
        self.assertEquals(repo.day_interval, 5, 'return day interval is not valid')
        self.assertEquals(repo.hr_interval, 0, 'return hour interval is not valid')
        self.assertEquals(repo.grid_source_data['t'][0].file_name, 'data/srcdata_low_res_C/ORCA025-N206_19791101d05T.nc',
                          'first t grid file name in the directory' )
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()