'''
Created on 14/8/2014

@author: victor
'''
import unittest
from hivprotmut.mutation.pdbPseudoMutation import PdbPseudoMutation
import cStringIO
import hivprotmut.mutation.test.data as test_data

class TestMutation(unittest.TestCase):


    def test_get_mutated_residues(self):
        master_sequence =  "PQVTLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF"
        mutated_sequence = "PQVTLWQGPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKLKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRPLLTQIGCTLNF"
        #                          X                                   X                                           X
        
        expected_mutated = {8:"GLY", 44:"LEU", 88:"PRO"}
        self.assertDictEqual(expected_mutated, PdbPseudoMutation.get_mutated_residues(master_sequence, mutated_sequence))                           
    
    def test_process_pdb(self):
        input_pdb_handler = cStringIO.StringIO(test_data.small_model)
        output_pdb_handler = cStringIO.StringIO()
        
        mutations = {4:"ALA"}
        parameters = {}
        PdbPseudoMutation.process_pdb_handlers(input_pdb_handler, 
                                               output_pdb_handler, 
                                               mutations, 
                                               parameters)
        self.assertEqual(test_data.mutated_small_model, output_pdb_handler.getvalue())
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()