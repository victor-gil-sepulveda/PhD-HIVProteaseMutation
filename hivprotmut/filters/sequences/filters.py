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
    
class OnlyOneAlignmentPerStructure(object):
    
    def __init__(self):
        pass
    
    @classmethod
    def filter(cls, alignments):
        """
        Brute force!
        """
        already_seen_ids = []
        filtered = []
        for alignment_info in alignments:
            pdb_id = alignment_info["pdb"]["id"]
            if not pdb_id in already_seen_ids:
                already_seen_ids.append(pdb_id)
                filtered.append(alignment_info)
        return filtered
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
