"""
Created on 14/8/2014

@author: victor
"""
import os


class PWCommands(object):
        
    PW_EXEC = "%s %s %s %s"
    
    def __init__(self):
        pass
    
    @classmethod
    def apply_protein_wizard(cls, input_pdb, parameters):
        os.system(PWCommands.PW_EXEC%(parameters["exec"], 
                                      parameters["options"],
                                      input_pdb,
                                      parameters["output_file"]))
