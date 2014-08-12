"""
Created on 12/8/2014

@author: victor
"""

import xml.etree.ElementTree
import json

class BlastOutputParser(object):


    def __init__(self, blast_output_file):
        tree =  xml.etree.ElementTree.parse(blast_output_file)
        self.alignments = self.parse_hits(tree)
        
    def parse_hits(self, tree):
        root = tree.getroot()
        pdbs = []
        for hit in root.iter("Hit"):
            acc_id =  hit.find("Hit_accession").text
            hit_info ={
                       "pdb": {"id": acc_id[:4]}, 
                       "hit_chain": acc_id[5:6],
                       "query_seq": [seq.text for seq in hit.iter("Hsp_qseq")][0],
                       "hit_seq__": [seq.text for seq in hit.iter("Hsp_hseq")][0], 
                       "score": float([seq.text for seq in hit.iter("Hsp_score")][0]),
                       "gaps": int([seq.text for seq in hit.iter("Hsp_gaps")][0])
                       }
            pdbs.append(hit_info)
        return pdbs
    
    def save(self, file_name):
        open(file_name, "w").write(json.dumps(self.alignments, 
                                              sort_keys=False, 
                                              indent=4, 
                                              separators=(',', ': ')))
        