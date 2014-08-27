'''
Created on 27/8/2014

@author: victor
'''
import unittest
import cStringIO
import hivprotmut.test.data as test_data
from hivprotmut.tools import get_all_ids_from_file_handler

class TestTools(unittest.TestCase):


    def test_get_pdb_list(self):
        expected_ids = ["2wl0", "1sgu", "1izi", "1d4s", "2fdd", "3s43", "3nu3", "3em6"]
        database = cStringIO.StringIO(test_data.id_file_contents)
        self.assertItemsEqual(expected_ids, get_all_ids_from_file_handler(database))

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_get_pdb_listName']
    unittest.main()