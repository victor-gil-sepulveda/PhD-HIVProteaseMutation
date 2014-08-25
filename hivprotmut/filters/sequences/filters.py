"""
created on 14/8/2014

@author: victor
"""

class ExactlyThisLengthFilter(object):
    
    def __init__(self):
        pass
    
    @classmethod
    def filter(cls, alignments, seq_len):
        return [alignment_info for alignment_info in alignments if len(alignment_info["query_seq"]) == seq_len]
  
class NoGapsFilter(object):
    
    def __init__(self):
        pass
    
    @classmethod
    def filter(cls, alignments):
        return [alignment_info for alignment_info in alignments if alignment_info["gaps"] == 0]

class SequenceAlignmentFilter(object):
    """
    Alignment filters only use alignment information (e.g. sequence)
    """
    def __init__(self):
        self.filters = []
    
    def add_filter(self, this_filter, params = None):
        self.filters.append((this_filter, params))
    
    def filter(self, alignments):
        filtered_alignments = alignments
        for filter_class, params  in self.filters:
            if params is None:
                filtered_alignments = filter_class.filter(filtered_alignments)
            else:
                filtered_alignments = filter_class.filter(filtered_alignments, params)
        return filtered_alignments
