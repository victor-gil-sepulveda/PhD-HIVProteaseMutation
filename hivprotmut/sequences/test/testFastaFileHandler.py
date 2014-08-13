'''
Created on 13/8/2014

@author: victor
'''
import unittest
import hivprotmut.sequences.test.data as data
import os
from hivprotmut.sequences.fastaFile import FastaFile, FastaFileHandler
import cStringIO

class TestFasta(unittest.TestCase):
    
    def test_open_and_close(self):
        fasta_handler = FastaFile.open(os.path.join(data.__path__[0], "HIV.fasta"), "r")
        
        expected_contents = open(os.path.join(data.__path__[0], "HIV.fasta"),"r").read()
        fasta_contents = fasta_handler.handler.read()
        fasta_handler.close()
        
        self.assertEqual(expected_contents, fasta_contents)
        self.assertRaises(IOError, fasta_handler.close)

    def test_write(self):
        sequence = open(os.path.join(data.__path__[0], "HIV.seq"), "r").read()
        expected_contents = open(os.path.join(data.__path__[0], "HIV.fasta"),"r").read()
        handler = cStringIO.StringIO()
        fasta_handler = FastaFileHandler(handler)
        fasta_handler.write("HIV", sequence)
        self.assertEqual(expected_contents, handler.getvalue())
        fasta_handler.close()

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()