import prody

class NumChainsIs(object):
    def __init__(self):
        pass
    
    @classmethod
    def is_filtered(cls, pdb, num_chains):
        """
        Checks if the structure has at least n proteic chains.
        """
        hw = prody.HierView(pdb.select("protein"))
        return hw.numChains() != num_chains

class NumChainsIsAtLeast(object):
    def __init__(self):
        pass
    
    @classmethod
    def is_filtered(cls, pdb, min_num_chains):
        """
        Checks if the structure has at least n proteic chains.
        """
        hw = prody.HierView(pdb.select("protein"))
        return hw.numChains() < min_num_chains


class EqualChainSequences(object):
    def __init__(self):
        pass
    
    @classmethod
    def must_be_filtered(cls, pdb):
        """
        Checks if the structure has at least n proteic chains.
        """
        hw = prody.HierView(pdb.select("protein"))
        
        return len(set([chain.getSequence() for chain in hw.iterChains()])) != 1

class NoResiduesNamed(object):
    def __init__(self):
        pass
    
    @classmethod
    def must_be_filtered(cls, pdb, residue_names):
        """
        pdb is a prody pdb.
        """
        for residue_name in residue_names:
            if residue_name in pdb.getResnames():
                return True
        return False
    
class CrystalHasLigandAndWater(object):
    
    
    def __init__(self):
        pass
    
    @classmethod
    def must_be_filtered(cls, pdb, min_ligand_atoms):
        """
        The structure must have ligand and waters. We expect that a ligand has at least
        MIN_LIGAND_ATOMS, but is impossible to know wether is the expected ligand or 
        just a bunch of ions.
        """
        ligand = pdb.select("hetero and not water")
        waters = pdb.select("water")
        return ligand is None or ligand.numAtoms() < min_ligand_atoms or waters is None 
        
class StructureAlignmentFilter(object):
    """
    Alignment filters use pdb information and 
    """
    def __init__(self):
        self.filters = []
    
    def add_filter(self, this_filter, params = None):
        self.filters.append((this_filter, params))
    
    def must_be_filtered(self, pdb):
        """
        Checks if a pdb must be filtered or not.
        """
        filtered = True
        failed = []
        
        for filter_class, params  in self.filters:
            if params is None:
                not_ok = filter_class.must_be_filtered(pdb)
            else:
                not_ok = filter_class.must_be_filtered(pdb, params)
            
            if not_ok:
                failed.append(filter_class.__name__)
            print "Filter: ",filter_class.__name__, not_ok
            filtered = filtered and not_ok
        return failed
