"""
Created on 25/8/2014

@author: victor
"""

import prody
import numpy

def curate_struct(initial_pdb, pdb_alignment):
    """
    Returns the "curated" pdb.
    """
    # Get chain info (without ligand or waters)
    hw = prody.HierView(initial_pdb.select("protein"))
     
    pdb_alignment["pdb"]["num_chains"] = hw.numChains()
    chain_lengths = []
    for chain in hw.iterChains():
        pdb_alignment["pdb"]["chain_"+chain.getChid()] = chain.getSequence()
        chain_lengths.append((len(chain.getSequence()), chain.getChid())) 
    
    # Sometimes we have complexes attached to the chain, sometimes we have
    # more than 2 chain description even if they cover the same space
    # We only have to leave two of that chains.
    leave_chains = sorted(chain_lengths)[0:2]
    leave_chains = [chain_id for _, chain_id in leave_chains]
    
    # Detect important waters and keep them
    water_structs = {}
    for chain in hw.iterChains():
        # Template residue is Ile 50, but it is not conserved
        # we have to index it by number of residue, taking into account
        # the offset if any
        # we need residue 50 and not residue with resnum 50 (the 49th)
        i = 0
        for residue in chain.iterResidues():
            if i == 49:
                break
            i = i+1
        
        residue_com = prody.calcCenter(residue)
        
        # Identify closer water 
        waters = initial_pdb.select("name O and water").copy()
        distances = numpy.sqrt(((residue_com - waters.getCoords())**2).sum(axis=1))
        
        min_dist = numpy.min(distances)
        min_dist_index = numpy.where(distances == min_dist)
        water_id = "%d:%s"%(waters.getResnums()[min_dist_index], waters.getChids()[min_dist_index][0])
        # we use a dict in order to get rid of repeats
        water_structs[water_id] = initial_pdb.select("resnum %d"%waters.getResnums()[min_dist_index]).copy()
    
    prot_struct = initial_pdb.select("protein chain %s"%(" ".join(leave_chains))).copy()
    ligand_struct = initial_pdb.select("hetero not water").copy()
    tmp_struct = prot_struct + ligand_struct
    pdb_alignment["pdb"]["waters"] = []
    for water_id in water_structs:
        # Keep track of added waters
        pdb_alignment["pdb"]["waters"].append(water_id)
        tmp_struct = tmp_struct + water_structs[water_id]
    
    return tmp_struct
