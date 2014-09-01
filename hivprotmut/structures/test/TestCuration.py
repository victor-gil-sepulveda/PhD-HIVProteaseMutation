'''
Created on 1/9/2014

@author: victor
'''
import unittest

import hivprotmut.structures.test.data as test_data
import os
from hivprotmut.tools import get_pdb_from_remote_or_db
from hivprotmut.structures.pdbcuration import curate_struct,\
    process_water_structures

class TestCuration(unittest.TestCase):

    # Normal, regression
    
    # Special "difficult" cases
    # 1w5y has more than 2 waters
    
    def test_NoMoreThan2Waters(self):
        pdb, path = get_pdb_from_remote_or_db("1w5y", "all", test_data.__path__[0])
        alignment = {"pdb":{}}
        os.remove(path)
        
        waters =  process_water_structures(pdb, ["A", "B"])
        expected_lengths = {
                            "2036:B":1,
                            "2035:B":1
                            }
        for wat_id in expected_lengths:
            self.assertEqual(expected_lengths[wat_id], len(waters[wat_id]))


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()