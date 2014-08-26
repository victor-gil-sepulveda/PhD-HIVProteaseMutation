"""
Created on 12/8/2014

@author: victor
"""
import os
import hivprotmut.tools as tools
import prody
from hivprotmut.external.blast.blastpCommands import BlastpCommands
import sys
import json
from hivprotmut.filters.sequences.filters import SequenceAlignmentFilter, NoGapsFilter,\
    ExactlyThisLengthFilter
from hivprotmut.filters.structures.filters import StructureAlignmentFilter,\
    NoResiduesNamed
from hivprotmut.structures.pdbcuration import curate_struct
from hivprotmut.tools import save_json
from datetime import datetime

if __name__ == '__main__':
    
    parameters = json.loads(tools.remove_comments(open(sys.argv[1]).read()))

    log = open(parameters["global"]["log_name"],"w")
    log.write(datetime.now().strftime("Started on %A, %d. %B %Y %I:%M%p\n"))
    
    tools.create_folder(parameters["global"]["structure_database"])
    

    ##########
    # Find alignments
    ##########
    alignments = BlastpCommands.find_closest_sequences("HIV.fasta", 
                                                       parameters["blastp"])
    log.write("Found %d alignments\n"%(len(alignments)))
    
    ##########
    # Filter them
    ##########
    # Get the ids of pdbs without gap (backbones must be equal)
    al_filter = SequenceAlignmentFilter()
    al_filter.add_filter(NoGapsFilter)
    al_filter.add_filter(ExactlyThisLengthFilter, 99)
    filtered_alignments = al_filter.filter(alignments)
    # IMPROVEMENT : LEAVE STRUCTURES WHERE THE GAPS ARE AT THE BEGGINING OR THE END TO SOME EXTENT
    # we need to store the offset
    log.write("We have filtered %d structures because of their sequence features.\n"%(len(alignments) - len(filtered_alignments)))
    

    ##########
    # Curate structures
    ##########    
    structure_filter = StructureAlignmentFilter()
    structure_filter.add_filter(NoResiduesNamed, 
                                parameters["pdb_preparation"]["forbidden_residues"])
    for alignment in filtered_alignments:
        pdb, pdb_path = tools.get_pdb(alignment["pdb"]["id"], 
                                      parameters["pdb_preparation"]["load_selection"])
        if not structure_filter.is_filtered(pdb):
            alignment["rejected"] = False
            tmp_path = os.path.join(parameters["global"]["structure_database"], 
                                    os.path.basename(pdb_path))
            os.remove(pdb_path)
            curated_pdb = curate_struct(pdb, alignment)
            prody.writePDB(tmp_path+".prot_lig_water", curated_pdb)
        else:
            os.remove(pdb_path)
            alignment["rejected"] = True
            log.write("PDB: %s not processed because it has a residue name in %s.\n"%(alignment["pdb"]["id"],
                                                                                str(parameters["pdb_preparation"]["forbidden_residues"])))
    
    
    BlastpCommands.create_database_from_alignments(filtered_alignments,
                                                  parameters["blast_database_creation"])
    save_json(filtered_alignments, 
              parameters["blastp"]["alignments_file"])
    
    log.write(datetime.now().strftime("Finished on %A, %d. %B %Y %I:%M%p\n"))
    log.close()
    
    # TODO: PROCESS ALT LOCS? I NEED AN EXAMPLE
    # TODO: STORE AA NAME IN POS 50
    # TODO: FIND GAP EXAMPLE. 
    # TODO: ADD TOLERANCE TO BEGINNING AND ENDING GAPS
    
    # Keep one water or two?
    # Ile 50 no siempre se conserva!
    # Mirar dist de centro de masas? o mejor de un atomo en concreto?
    ###############################
    ### PREPROCESSING CONFIG
    # MAX_ALLOWED_BEGINNING_GAPS = 3
    # MAX_ALLOWED_ENDING_GAPS = 3
    ###############################
    
    