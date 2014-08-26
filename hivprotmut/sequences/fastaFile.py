"""
Created on 13/8/2014

@author: victor
"""

def sequence_iterator(sequence, max_chars):
    """
    Gets chunks of max_char from a sequence.
    """
    beginning = 0
    ending = beginning + max_chars
    while not beginning >= len(sequence):
        yield sequence[beginning:ending]
        beginning = ending
        ending = min(beginning + max_chars, len(sequence))


class FastaFile(object):
    
    def __init__(self):
        pass
    
    @classmethod
    def open(cls, filename, mode = "w"): 
        handler = open(filename, mode)
        return FastaFileHandler(handler)
        
class FastaFileHandler(object):
    MAX_COLUMNS = 80 
    
    def __init__(self, handler):
        self.handler = handler
    
    def close(self):
        if self.handler is not None:
            self.handler.close()
            self.handler = None
        else:
            raise IOError((0,"Handler already close."))
        
    def write(self, item_id, sequence):
        if self.handler is not None:
            # Header
            self.handler.write(">%s\n"%(item_id))
            # Sequence
            for sequnce_chunk in sequence_iterator(sequence, FastaFileHandler.MAX_COLUMNS):
                self.handler.write("%s\n"%sequnce_chunk)
        else:
            raise IOError((0,"Handler already close."))
        
        
        