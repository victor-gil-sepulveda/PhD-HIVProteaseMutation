'''
Created on 25/8/2014

@author: victor
'''
import unittest
from hivprotmut.filters.sequences.filters import ExactlyThisLengthFilter,\
    NoGapsFilter, SequenceAlignmentFilter


class TestSequenceFilters(unittest.TestCase):

    def test_nogap_filter(self):
        sequence = "en un lugar de la mancha de cuyo nombre ... "
        alignments = [
                      { 
                       "query_seq":   sequence,
                       "gaps" : 10                 
                      },
                      { 
                       "query_seq":   sequence[0:2],
                       "gaps": 0             
                      },
                      { 
                       "query_seq":   sequence[0:30],
                       "gaps":   7          
                      }]
        self.assertItemsEqual( NoGapsFilter.filter(alignments),
                [{'query_seq': 'en', 'gaps': 0}])
    
    def test_seq_len_filter(self):
        sequence = "en un lugar de la mancha de cuyo nombre ... "
        alignments = [
                      { 
                       "query_seq":   sequence                 
                      },
                      { 
                       "query_seq":   sequence[0:30]                 
                      },
                      { 
                       "query_seq":   sequence[0:20]                 
                      }]
        self.assertItemsEqual( ExactlyThisLengthFilter.filter(alignments, 44),
                              [{'query_seq': 'en un lugar de la mancha de cuyo nombre ... '}])
        self.assertItemsEqual( ExactlyThisLengthFilter.filter(alignments, 30),
                               [{'query_seq': 'en un lugar de la mancha de cu'}])
        self.assertItemsEqual( ExactlyThisLengthFilter.filter(alignments, 20),
                               [{'query_seq': 'en un lugar de la ma'}])
    
    def test_more_than_one_filter(self):
        sequence = "en un lugar de la mancha de cuyo nombre ... "
        sequence2 = "EnUnLugarDelamanchadecuyonombre...__________"
        alignments = [
                      { 
                       "query_seq":   sequence,
                       "gaps" : 10                 
                      },
                      { 
                       "query_seq":   sequence2,
                       "gaps" : 0                 
                      },
                      { 
                       "query_seq":   sequence[0:2],
                       "gaps": 0             
                      },
                      { 
                       "query_seq":   sequence2[0:2],
                       "gaps": 3             
                      },
                      { 
                       "query_seq":   sequence[0:30],
                       "gaps":   7          
                      }]
        al_filter = SequenceAlignmentFilter()
        al_filter.add_filter(NoGapsFilter)
        al_filter.add_filter(ExactlyThisLengthFilter, 44)
        filtered_alignments = al_filter.filter(alignments)
        self.assertItemsEqual(  filtered_alignments,
                                [{'query_seq': 'EnUnLugarDelamanchadecuyonombre...__________', 'gaps': 0}])

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()