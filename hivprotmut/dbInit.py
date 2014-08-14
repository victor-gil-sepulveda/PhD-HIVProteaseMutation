"""
Created on 12/8/2014

@author: victor
"""
import os
import hivprotmut.tools as tools
import prody
import numpy
from hivprotmut.external.blast.blastpCommands import BlastpCommands
import sys

###############################
### DB CONFIG
PDB_TMP_DATABASE_FOLDER = "tmp_db"
PDB_DATABASE_FOLDER = "db"
###############################

###############################
### PREPROCESSING CONFIG
LOAD_SELECTION_STRING = "all"
MAX_ALLOWED_BEGINNING_GAPS = 3
MAX_ALLOWED_ENDING_GAPS = 3
###############################


if __name__ == '__main__':
    
    parameters = tools.remove_comments(open(sys.argv[1]).read())
    
    tools.create_folder(PDB_TMP_DATABASE_FOLDER)
    tools.create_folder(PDB_DATABASE_FOLDER)
    
    alignments = BlastpCommands.find_closest_sequences("HIV.fasta", 
                                                       parameters["blastp"])
    
    print "Found %d alignments"%(len(alignments))
 
    # Get the ids of pdbs without gap (backbones must be equal)
    filtered_alignments = [alignment_info for alignment_info in alignments[0:5] if alignment_info["gaps"] == 0]
    # IMPROVEMENT : LEAVE STRUCTURES WHERE THE GAPS ARE AT THE BEGGINING OR THE END TO SOME EXTENT
    # we need to store the offset
    print "We have filtered %d structures because the alignment had gaps."%(len(alignments) - len(filtered_alignments))
     
    # Get the pdbs
    for alignment in filtered_alignments:
        pdb, pdb_path = tools.get_pdb(alignment["pdb"]["id"], LOAD_SELECTION_STRING)
        tmp_path = os.path.join(PDB_TMP_DATABASE_FOLDER, os.path.basename(pdb_path))
        os.remove(pdb_path)
         
        # Get chain info (without ligand or waters)
        hw = prody.HierView(pdb.select("protein"))
         
        alignment["pdb"]["num_chains"] = hw.numChains()
        for chain in hw.iterChains():
            alignment["pdb"]["chain_"+chain.getChid()] = chain.getSequence()
         
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
            waters = pdb.select("name O and water").copy()
            distances = numpy.sqrt(((residue_com - waters.getCoords())**2).sum(axis=1))
             
            min_dist = numpy.min(distances)
            min_dist_index = numpy.where(distances == min_dist)
            water_id = "%d:%s"%(waters.getResnums()[min_dist_index], waters.getChids()[min_dist_index][0])
            # we use a dict in order to get rd of repeats
            water_structs[water_id] = pdb.select("resnum %d"%waters.getResnums()[min_dist_index]).copy()
         
        prot_struct = pdb.select("protein").copy()
        ligand_struct = pdb.select("hetero not water").copy()
        tmp_struct = prot_struct + ligand_struct
        alignment["pdb"]["waters"] = []
        for water_id in water_structs:
            # Keep track of added waters
            alignment["pdb"]["waters"].append(water_id)
            tmp_struct = tmp_struct + water_structs[water_id]
         
        prody.writePDB(tmp_path+".prot_lig_water", tmp_struct)
         
    # TODO: PROCESS ALT LOCS? I NEED AN EXAMPLE
    # TODO: STORE AA NAME IN POS 50
    # TODO: FIND GAP EXAMPLE. 
    # TODO: ADD TOLERANCE TO BEGINNING AND ENDING GAPS
    
    # Keep one water or two?
    # Ile 50 no siempre se conserva!
    # Mirar dist de centro de masas? o mejor de un atomo en concreto?
    
    BlastpCommands.create_database_from_alignments(filtered_alignments,
                                                  parameters["blast_database_creation"])
