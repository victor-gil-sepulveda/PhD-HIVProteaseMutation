"""
Created on 12/8/2014

@author: victor
"""

import prody
import os
import errno
import json


def get_pdb(pdb_id, selection):
    """
    Downloads a pdb from the Protein Data Bank (if necessary) and removes all models so that it only has one
    model.

    :param pdb_id: A 4 letter pdb id

    :return: The downloaded pdb data structure.
    """
    # Download pdb
    path = prody.fetchPDB(pdb_id, compressed=False)
    
    # Get pdb data structure
    pdb = prody.parsePDB(path)
    pdb = pdb.select(selection).copy()
    number_of_models = pdb.numCoordsets()
    
    # Delete all coordsets but coordset 0
    [pdb.delCoordset(1) for _ in range(1, number_of_models)]
    return pdb, path

def create_folder(directory_path, ensure_writability=False):
    """
    Creates a directory (with subdirs) if it doesn't exist.

    :param directory_path: the path of the directory and subdirectories to be created.

    :return: True if it was possible to create the folder, False otherwise.
    """
    if ensure_writability:
        if not os.access(os.path.dirname(directory_path), os.W_OK):
            return False
    try:
        os.makedirs(directory_path)
        return True
    except OSError, e:
        if e.errno != errno.EEXIST:
            raise
    return False

def save_json(items, path):
    open(path, "w").write(json.dumps(items, 
                                              sort_keys=False, 
                                              indent=4, 
                                              separators=(',', ': ')))
