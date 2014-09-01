"""
Created on 25/8/2014

@author: victor
"""

import prody
import numpy

class CurationSelections():
    LIGAND_SELECTION = "hetero not water not ion"
    PROTEIN_CHAIN_TEMPLATE = "protein chain %s"
    
    def __init__(self):
        pass

def choose_main_chains(initial_pdb):
    """
    We can have complexes attached to the chain or even duplicated chains
    that cover the same space (ex. in the same model, A and B are one structure
    and C and B form a duplicated protein). We only have to leave two of that 
    main chains, and that's what this function does :) .
    
    :param initial_pdb: The pdb (prody structure) we want to extract the chains.
    
    :return: An array containing the chain ids of the main chains.
    """
    hw = prody.HierView(initial_pdb.select("protein"))
    chain_lengths = []
    for chain in hw.iterChains():
        chain_lengths.append((len(chain.getSequence()), chain.getChid())) 
    
    leave_chains = sorted(chain_lengths)[-2:]
    leave_chains = [chain_id for _, chain_id in leave_chains]
    return leave_chains

def process_water_structures(initial_pdb, main_chains ):
    """
    Detects the waters we have to keep (important for the simulation) and returns 
    a structure holding them.
    Important waters are the ones closer to Template residue 50(Ile), the aa is not 
    but it is not guaranteed to be conserved, which means we have to rely into the 
    residue number to choose it, and take any offset into account if needed.
    
    :param initial_pdb: The pdb (prody structure) we want to extract the chains.
    
    :return: A dictionary indexed by the water id (res. num. + chain id) holding the prody pdb
    structure of that water.
    """
    hw = prody.HierView(initial_pdb.select("protein"))
    water_structs = {}
    for chain in hw.iterChains():
        if chain.getChid() in main_chains:
            # We cannot do a direct selection, instead we iterate
            for i, residue in enumerate(chain.iterResidues()):
                if i == 50: # 50th residue
                    break
            
            residue_com = prody.calcCenter(residue)
            
            # Identify closer water 
            waters = initial_pdb.select("name O and water")
            if waters is not None:
                distances = numpy.sqrt(((residue_com - waters.getCoords())**2).sum(axis=1))
                
                min_dist = numpy.min(distances)
                min_dist_index = numpy.where(distances == min_dist)
                water_resnum = waters.getResnums()[min_dist_index]
                water_chid = waters.getChids()[min_dist_index][0]
                water_id = "%d:%s"%(water_resnum, water_chid)
                # We use a dict in order to get rid of repeats
                selection_string = "resnum %d and chain %s"%(water_resnum,
                                                             water_chid)
                water_structs[water_id] = initial_pdb.water.select(selection_string).copy()
    return water_structs

def curate_struct(initial_pdb, main_chains, pdb_alignment, parameters):
    """
    Returns the "curated" pdb. A curated pdb has potentially 2 waters around residue 
    50 of each chain, a ligand and two main (symmetric) chains; everything else must be 
    deleted. This function will work even in the case that the 2 later are not present, 
    which can happen when processing any of the "mandatory" structures (those can pass 
    the filters automatically).
    
    :param initial_pdb: The prody pdb structure we want to extract the chains.
    
    :return: The "curated" pdb.
    """
    # Get chain info (without ligand or waters)
    hw = prody.HierView(initial_pdb.select("protein"))
    pdb_alignment["pdb"]["num_chains"] = hw.numChains()

    # Pick main chains
    prot_struct = initial_pdb.select(CurationSelections.PROTEIN_CHAIN_TEMPLATE%(" ".join(main_chains))).copy()
    
    # Add the ligand (if found), must be part of other chains (not main_chains)
    ligand_struct = initial_pdb.select(CurationSelections.LIGAND_SELECTION)
    if ligand_struct is not None and ligand_struct.numAtoms() >= parameters["min_ligand_atoms"]:
        tmp_struct = prot_struct + ligand_struct.copy()
    else:
        tmp_struct = prot_struct
    
    # Add "important" waters, if found
    water_structs = process_water_structures(initial_pdb, main_chains)
    pdb_alignment["pdb"]["waters"] = water_structs.keys() # Keep track of added waters in the alignment file
    for water_id in water_structs:
        tmp_struct = tmp_struct + water_structs[water_id]
    
    return tmp_struct
