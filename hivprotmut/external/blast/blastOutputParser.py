"""
Created on 12/8/2014

@author: victor
"""

import xml.etree.ElementTree

class BlastOutputParser(object):


    def __init__(self, blast_output_file, only_one_alignment_per_hit):
        tree =  xml.etree.ElementTree.parse(blast_output_file)
        self.alignments = self.parse_hits(tree, only_one_alignment_per_hit)
        
    def parse_hits(self, tree, only_one_alignment_per_hit):
        root = tree.getroot()
        alignments = []
        for hit in root.iter("Hit"):
            acc_id =  hit.find("Hit_accession").text
            for hsp in hit.iter("Hsp"):
                hit_info ={
                           "pdb": {"id": acc_id[:4]}, 
                           "hit_chain": acc_id[5:6],
                           "query_seq": hsp.find("Hsp_qseq").text,
                           "hit_seq__": hsp.find("Hsp_hseq").text, 
                           "midline__": hsp.find("Hsp_midline").text,
                           "score": float(hsp.find("Hsp_score").text),
                           "gaps": int(hsp.find("Hsp_gaps").text)
                           }
                
                alignments.append(hit_info)
                
                if only_one_alignment_per_hit:
                    break # We only store one hsp
        return alignments
        
        