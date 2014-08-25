

class NoResiduesNamed(object):
    def __init__(self):
        pass
    
    @classmethod
    def is_filtered(cls, pdb, residue_names):
        """
        pdb is a prody pdb.
        """
        for residue_name in residue_names:
            if residue_name in pdb.getResnames():
                return True
        return False
        
class StructureAlignmentFilter(object):
    """
    Alignment filters use pdb information and 
    """
    def __init__(self):
        self.filters = []
    
    def add_filter(self, this_filter, params = None):
        self.filters.append((this_filter, params))
    
    def is_filtered(self, pdb):
        filtered = True
        
        for filter_class, params  in self.filters:
            if params is None:
                filtered = filtered and filter_class.is_filtered(pdb)
            else:
                filtered = filtered and filter_class.is_filtered(pdb, params)
        
        return filtered
