"""
Created on 12/8/2014

@author: victor
"""
import os
import prody
import sys
import json
from datetime import datetime
import hivprotmut.tools as tools
from hivprotmut.sequences.fastaFile import FastaFile
from hivprotmut.structures.pdbcuration import curate_struct
from hivprotmut.external.blast.blastpCommands import BlastpCommands
from hivprotmut.filters.sequences.filters import SequenceAlignmentFilter, NoGapsFilter,\
    ExactlyThisLengthFilter, OnlyOneAlignmentPerStructure
from hivprotmut.filters.structures.filters import StructureAlignmentFilter,\
    NoResiduesNamed, CrystalHasLigandAndWater, EqualChainSequences,\
    NumChainsIs, NumChainsIsAtLeast

if __name__ == '__main__':
    
    parameters = json.loads(tools.remove_comments(open(sys.argv[1]).read()))

    log = open(parameters["global"]["log_name"],"w")
    log.write(datetime.now().strftime("Started on %A, %d. %B %Y %I:%M%p\n"))
    
    tools.create_folder(parameters["global"]["curated_structure_database"])
    
    ##########
    # Find alignments
    ##########
    # Prepare fasta file with input sequence
    input_fasta = FastaFile.open(parameters["query"]["fasta_file"])
    input_fasta.write("QUERY", parameters["query"]["sequence"])
    input_fasta.close()
    # Use blast
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
    al_filter.add_filter(OnlyOneAlignmentPerStructure)
    filtered_alignments = al_filter.filter(alignments)
    # IMPROVEMENT : LEAVE STRUCTURES WHERE THE GAPS ARE AT THE BEGGINING OR THE END TO SOME EXTENT
    # we need to store the offset
    log.write("We have filtered %d structures because of their sequence features.\n"%(len(alignments) - len(filtered_alignments)))
    
    ##########
    # Curate structures
    ##########
    structure_filter = StructureAlignmentFilter()
    
    structure_filter.add_filter(CrystalHasLigandAndWater)
    #structure_filter.add_filter(NumChainsIs, 2)
    structure_filter.add_filter(NumChainsIsAtLeast, 2)
    structure_filter.add_filter(EqualChainSequences)
    structure_filter.add_filter(NoResiduesNamed, 
                                parameters["pdb_preparation"]["forbidden_residues"])
    
    structure_db_path = parameters["pdb_preparation"]["structure_database_path"] if "structure_database_path" in parameters["pdb_preparation"] else ""
    log.write("While processing the structures:\n")
    for alignment in filtered_alignments:
        pdb, pdb_path = tools.get_pdb_from_remote_or_db(alignment["pdb"]["id"], 
                                      parameters["pdb_preparation"]["load_selection"],
                                      structure_db_path)
        
        failed_filters = structure_filter.is_filtered(pdb)
        if len(failed_filters) == 0:
            alignment["rejected"] = False
            tmp_path = os.path.join(parameters["global"]["curated_structure_database"], 
                                    os.path.basename(pdb_path))
            
            curated_pdb = curate_struct(pdb, alignment)
            prody.writePDB(tmp_path, curated_pdb)
        else:
            alignment["rejected"] = True
            log.write("\t - PDB %s is not going to be stored because it has failed this filters: %s.\n"%(alignment["pdb"]["id"],
                                                                                str(failed_filters)))
        
        if structure_db_path is not None:# We do not want to delete our database!
            os.remove(pdb_path)
    
    # Blast DB creation
    BlastpCommands.create_database_from_alignments(filtered_alignments,
                                                  parameters["blast_database_creation"])
    tools.save_json(filtered_alignments, 
              parameters["blastp"]["alignments_file"])
    
    log.write(datetime.now().strftime("Finished on %A, %d. %B %Y %I:%M%p\n"))
    log.close()
