"""
Created on 14/8/2014

@author: victor
"""
import os


class PlopCommands(object):
        
    PLOP_CONTROL_TEMPLATE = """load pdb %s het yes ions yes wat yes 
    write pdb %s
    """
    
    PLOP_EXECUTION = "%s %s"
    
    def __init__(self):
        pass
    
    @classmethod
    def reconstruct_pdb(cls, input_pdb, output_pdb, parameters):
        script = PlopCommands.PLOP_CONTROL_TEMPLATE%(input_pdb, output_pdb)
        open(parameters["control_file"],"w").write(script)
        os.system(PlopCommands.PLOP_EXECUTION%(parameters["exec"], parameters["control_file"]))
    
