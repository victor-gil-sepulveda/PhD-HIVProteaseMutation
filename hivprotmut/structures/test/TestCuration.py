'''
Created on 1/9/2014

@author: victor
'''
import unittest

import hivprotmut.structures.test.data as test_data
import os
from hivprotmut.tools import get_pdb_from_remote_or_db
from hivprotmut.structures.pdbcuration import curate_struct,\
    process_water_structures, CurationSelections

class TestCuration(unittest.TestCase):

    # Normal, regression
    
    # Special "difficult" cases
    # 1w5y has more than 2 waters
    # 1izi has 2 ions
    
    def test_NoMoreThan2Waters(self):
        """
        Any structure has at most 2 waters.
        """
        pdb, path = get_pdb_from_remote_or_db("1w5y", "all", test_data.__path__[0])
        os.remove(path)
        
        waters =  process_water_structures(pdb, ["A", "B"])
        expected_lengths = {
                            "2036:B":1,
                            "2035:B":1
                            }
        for wat_id in expected_lengths:
            self.assertEqual(expected_lengths[wat_id], len(waters[wat_id]))
            
    def test_NoIons(self):
        """
        No ions must be present in the structure.
        """
        pdb, path = get_pdb_from_remote_or_db("1izi", "all", test_data.__path__[0])
        os.remove(path)
        self.assertItemsEqual(['Q50'], set(pdb.select(CurationSelections.LIGAND_SELECTION).getResnames()))
        
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()