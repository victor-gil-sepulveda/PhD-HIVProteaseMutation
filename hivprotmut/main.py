"""
Created on 12/8/2014

@author: victor
"""
import os
from hivprotmut.BlastOutputParser import BlastOutputParser
import hivprotmut.tools as tools
import prody
import math
import numpy

###############################
### GLOBAL CONFIG
MAX_SEQUENCES = 50
###############################

###############################
### BLAST CONFIG
BLASTP_EXEC = "blastp"
BLASTP_ARGS =  " -query %s -out %s -remote -db %s -max_target_seqs %d -outfmt 5"
BLAST_DB = "pdbaa"
###############################

###############################
### DB CONFIG
PDB_TMP_DATABASE_FOLDER = "tmp_db"
PDB_DATABASE_FOLDER = "db"
###############################

###############################
### PREPROCESSING CONFIG
LOAD_SELECTION_STRING = "all"
WATER_SEARCH_CUTOFF = 7 # 7 A
###############################


def find_closest_sequence_pdbs(fasta_file, output_file):
    """
    """
    os.system(BLASTP_EXEC + BLASTP_ARGS%(fasta_file, output_file, BLAST_DB, MAX_SEQUENCES))

if __name__ == '__main__':
    tools.create_folder(PDB_TMP_DATABASE_FOLDER)
    tools.create_folder(PDB_DATABASE_FOLDER)
    
    # Find close sequences 
#     find_closest_sequence_pdbs("hivprotmut/test/HIV.fasta", "sequences.xml")
    
    # Parse results
    bot =  BlastOutputParser("sequences.xml")
    bot.save("alignments.json")
    print "Found %d alignments"%(len(bot.alignments))

    # Get the ids of pdbs without gap (backbones must be equal)
    filtered_alignments = [alignment_info for alignment_info in bot.alignments if alignment_info["gaps"] == 0]
    # IMPROVEMENT : LEAVE STRUCTURES WHERE THE GAPS ARE AT THE BEGGINING OR THE END TO SOME EXTENT
    # we need to store the offset
    print "We have filtered %d structures because the alignment had gaps."%(len(bot.alignments) - len(filtered_alignments))
    
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
        
        water_structs = []
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
            water_structs.append(pdb.select("resnum %d"%waters.getResnums()[min_dist_index]).copy())
        
        prot_struct = pdb.select("protein").copy()
        ligand_struct = pdb.select("hetero not water").copy()
        tmp_struct = prot_struct + ligand_struct
        for water in water_structs:
            tmp_struct = tmp_struct + water
        
        prody.writePDB(tmp_path+".prot_lig_water", tmp_struct)
#             # keep one water or two?
# Ile 50 no siempre se conserva
# Mirar dist de centro de masas?

