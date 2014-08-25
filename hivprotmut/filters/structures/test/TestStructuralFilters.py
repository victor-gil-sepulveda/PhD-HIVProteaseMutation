'''
Created on 25/8/2014

@author: victor
'''
import unittest
from hivprotmut.filters.structures.filters import NoResiduesNamed,\
    StructureAlignmentFilter

class pdb_mock(object):
    def __init__(self, residues):
        self.residues = residues
        
    def getResnames(self):
        return self.residues

class TestStructuralFilters(unittest.TestCase):

    def test_NoResiduesNamed(self):
        pdb = pdb_mock(["ALA","ILE","ULE", "ELE", "OLE!!"])
        self.assertEqual(True, NoResiduesNamed.is_filtered(pdb, ["ELE"]))
        self.assertEqual(False, NoResiduesNamed.is_filtered(pdb, ["GLU"]))
        self.assertEqual(True, NoResiduesNamed.is_filtered(pdb, ["OLE","OLE!!"]))
        
        myfilter =  StructureAlignmentFilter()
        myfilter.add_filter(NoResiduesNamed, ["ELE"])
        self.assertEqual(True, myfilter.is_filtered(pdb))
        
        myfilter =  StructureAlignmentFilter()
        myfilter.add_filter(NoResiduesNamed, ["GLU"])
        self.assertEqual(False, myfilter.is_filtered(pdb))
        
        myfilter =  StructureAlignmentFilter()
        myfilter.add_filter(NoResiduesNamed, ["OLE","OLE!!"])
        self.assertEqual(True, myfilter.is_filtered(pdb))
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_NoResiduesNamed']
    unittest.main()