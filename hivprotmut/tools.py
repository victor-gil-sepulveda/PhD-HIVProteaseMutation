"""
Created on 12/8/2014

@author: victor
"""

import prody
import os
import errno
import json
import re


def get_pdb(pdb_id, selection):
    """
    Downloads a pdb from the Protein Data Bank (if necessary) and removes all models so that it only has one
    model.

    :param pdb_id: A 4 letter pdb id

    :return: The downloaded pdb prody data structure and the path to the downloaded file.
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

def get_pdb_from_remote_or_db(pdb_id, selection, source_folder = ""):
    """
    If called without "source_folder" the behaviour of this function is actually that
    of 'get_pdb'. If "source_folder" is defined (a user-defined folder) Prody will look
    the availability of the pdb there before trying to download it from rcbs pdb (slower).  
    """
    prody.pathPDBFolder(source_folder)
    return get_pdb(pdb_id, selection)

def get_ids_from_file_handler(database_file_handler):
    """
    Reads the pdb ids inside a file handler.
    
    :param database_file_path: The handler of the (opened) file holding the ids.
    
    :return: A list with all read ids (with no repeats).
    """
    all_ids = remove_comments(database_file_handler.read())
    all_ids = all_ids.split()
    all_ids = filter(lambda x: len(x) == 4, all_ids)
    all_ids = set([one_id.lower() for one_id in all_ids])
    return list(all_ids)

def get_ids_from_file(database_file_path):
    """
    Loads a file with PDB ids, filters and returns them.

    :param database_file_path: The path (and filename) of the text file.
    
    :return:  An array containing the read ids (with no repeats).
    """
    database = open(database_file_path, "r")
    all_ids = get_ids_from_file_handler(database)
    database.close()
    return all_ids

def get_ids_from_file_list(file_list):
    """
    Gets the pdb ids of different files and returns them without repeated elements.
    
    :param file_list: An array of file paths.
    
    :return: Non repeated ids in the input files.
    """
    all_ids = []
    for file_path in file_list:
        all_ids.extend(get_ids_from_file(file_path))
    return list(set(all_ids))

def get_all_ids_from_folder(folder):
    """
    Returns the pdb ids present in a folder (if the original file has a canonical name).
    """
    contents = os.listdir(folder)
    return dict([(x[0:4].lower(), x) for x in filter(lambda x: ".pdb" in x, contents)])

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
    """
    Saves a dictionary in JSON format.
    """
    open(path, "w").write(json.dumps(items, sort_keys=False, indent=4, separators=(',', ': ')))
    
def remove_comments(string):
    """
    Removes /**/ and // comments from a string (used with the control script).
    From http://stackoverflow.com/questions/2319019/using-regex-to-remove-comments-from-source-files
    """
    string = re.sub(re.compile("/\*.*?\*/",re.DOTALL ) ,"" ,string) 
    string = re.sub(re.compile("//.*?\n" ) ,"\n" ,string)
    return string


aa_dic_standard = {'asp': 'D', 'glu': 'E', 'lys': 'K', 'his': 'H', 'arg': 'R',
                 'gln': 'Q', 'asn': 'N', 'ser': 'S', 'asx': 'B', 'glx': 'Z',
                 'phe': 'F', 'trp': 'W', 'tyr': 'Y', 'gly': 'G', 'ala': 'A',
                 'ile': 'I', 'leu': 'L', 'cys': 'C', 'met': 'M', 'thr': 'T',
                 'val': 'V', 'pro': 'P', 'cyx': 'C', 'hid': 'H', 'hie': 'H',
                 'hip': 'H', 'unk': 'X', 'ace': 'X', 'nme': 'X'}

def get_protein_sequence(pdb):
    """
    Generates the 1 letter per residue sequence for a protein. Uses a dictionary that maps the 3 letter naming with the 1 letter naming convention
    Source:
        - Biskit (http://biskit.pasteur.fr/)

    @param pdb: A prody pdb data structure.

    @return: A string with the sequence of this protein.
    """
    # One-liner just for the sake of the challenge
    return "".join([aa_dic_standard[resname] if resname in aa_dic_standard else "X" for resname in
                    [residue.getResname().lower() for residue in prody.HierView(pdb).iterResidues()]])

def one_letter_residue_to_three(char_resid):
    """
    Converts a 1-letter residue id string into a 3-letter residue id list.
    http://stackoverflow.com/questions/483666/python-reverse-inverse-a-mapping
    """
    inv_map = {v:k for k, v in aa_dic_standard.items()}
    three_char_res_id = "UNK"
    try:
        three_char_res_id =  inv_map[char_resid]
    except Exception:
        pass
    return three_char_res_id
