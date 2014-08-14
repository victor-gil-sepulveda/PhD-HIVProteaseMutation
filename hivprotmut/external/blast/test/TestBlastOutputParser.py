'''
Created on 14/8/2014

@author: victor
'''
import unittest
from hivprotmut.external.blast.blastOutputParser import BlastOutputParser
import  hivprotmut.external.blast.test.data as test_data
import os

class Test(unittest.TestCase):
    
    def test_parse(self):
        """
        Tests if it can read the xml alignments.
        """
        expected_alignments = [
                               {
                                'score': 436.0, 
                                'query_seq': 'PQVTLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF', 
                                'midline__': 'PQ+TLW+RPLVTI+IGGQLKEALLDTGADDTV+EE +LPG WKPK IGGI GFIKVRQYDQI +EI GHKAIGTVLVGPTPVNIIGRNLLTQIG TLNF', 
                                'hit_seq__': 'PQITLWKRPLVTIRIGGQLKEALLDTGADDTVIEEXNLPGXWKPKXIGGIXGFIKVRQYDQIPVEIXGHKAIGTVLVGPTPVNIIGRNLLTQIGXTLNF', 
                                'gaps': 0, 
                                'hit_chain': 'A', 
                                'pdb': {
                                        'id': '3HZC'
                                        }
                                }
                               ]
        bo = BlastOutputParser(os.path.join(test_data.__path__[0], "two_hits_example"), True)
        
        self.assertEqual(len(bo.alignments), 1)
        self.assertDictEqual(expected_alignments[0], bo.alignments[0])
        
    def test_two_hits(self):
        """
        If an alignment has two hits, we only store the first.
        """
        expected_alignments = [{
          'score': 436.0, 
          'gaps': 0, 
          'query_seq': 'PQVTLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF', 
          'midline__': 'PQ+TLW+RPLVTI+IGGQLKEALLDTGADDTV+EE +LPG WKPK IGGI GFIKVRQYDQI +EI GHKAIGTVLVGPTPVNIIGRNLLTQIG TLNF', 
          'hit_seq__': 'PQITLWKRPLVTIRIGGQLKEALLDTGADDTVIEEXNLPGXWKPKXIGGIXGFIKVRQYDQIPVEIXGHKAIGTVLVGPTPVNIIGRNLLTQIGXTLNF', 
          'hit_chain': 'A', 
          'pdb': {
                  'id': '3HZC'
          }
        }, 
         
         {
          'score': 407.0, 
          'gaps': 0, 
          'query_seq': 'PQVTLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF', 
          'midline__': 'PQ+TLW+RPLVTI+IGGQLKEALLDTGADDTV+EE +LPG WKPK IGGIGGFIKVRQYDQI +EI GHKAIGTVLVGPTPVNIIGRNLLTQIG TLNF', 
          'hit_seq__': 'PQITLWKRPLVTIRIGGQLKEALLDTGADDTVIEEXNLPGXWKPKXIGGIGGFIKVRQYDQIPVEIXGHKAIGTVLVGPTPVNIIGRNLLTQIGXTLNF', 
          'hit_chain': 'A', 
          'pdb': {
                  'id': '3HZC'
          }
        }]

        bo = BlastOutputParser(os.path.join(test_data.__path__[0], "two_hits_example"), False)
        
        self.assertEqual(len(bo.alignments), 2)
        self.assertDictEqual(expected_alignments[0], bo.alignments[0])
        self.assertDictEqual(expected_alignments[1], bo.alignments[1])

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_two_hits']
    unittest.main()