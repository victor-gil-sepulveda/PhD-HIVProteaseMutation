"""
Created on 14/8/2014

@author: victor
"""
import os


class PlopCommands(object):
        
    PLOP_CONTROL_TEMPLATE = """data %s
    load pdb %s het %s ions yes wat yes opt %s
    write pdb %s
    """
    
    PLOP_EXECUTION = "%s %s"
    
    def __init__(self):
        pass
    
    @classmethod
    def reconstruct_pdb(cls, input_pdb, parameters):
        script = PlopCommands.PLOP_CONTROL_TEMPLATE%(parameters["data_path"], 
                                                     input_pdb, 
                                                     parameters["options"]["hetero"],
                                                     parameters["options"]["optimize"],
                                                     parameters["output_file"])
        open(parameters["control_file"],"w").write(script)
        os.system(PlopCommands.PLOP_EXECUTION%(parameters["exec"], parameters["control_file"]))
    
