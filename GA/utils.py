from rdkit.Chem import AllChem
from rdkit.Chem import MolFromSmiles
from rdkit.Chem import MolToSmiles
import rdkit
from rdkit import Chem
import logging

def set_logger(path):
    logger = logging.getLogger(__name__)
    streamHandler = logging.StreamHandler()
    logger.addHandler(streamHandler)

    fileHandler = logging.FileHandler(path)
    logger.addHandler(fileHandler)
    logger.setLevel(level=logging.DEBUG)

    return logger

def read_dataset(filename):
    '''Return a list of smiles contained in file filename

    Parameters:
    filename (string) : Name of file containg smiles seperated by '\n'

    Returns
    content  (list)   : list of smile string in file filename
    '''
    with open(filename) as f:
        content = f.readlines()
    content = [x.strip() for x in content]

    return content

def get_chunks(arr, num_processors, ratio):
    """
    Get chunks based on a list
    """
    chunks = []  # Collect arrays that will be sent to different processorr
    counter = int(ratio)
    for i in range(num_processors):
        if i == 0:
            chunks.append(arr[0:counter])
        if i != 0 and i<num_processors-1:
            chunks.append(arr[counter-int(ratio): counter])
        if i == num_processors-1:
            chunks.append(arr[counter-int(ratio): ])
        counter += int(ratio)
    return chunks

def sanitize_smiles(smi):
    '''Return a canonical smile representation of smi

    Parameters:
    smi (string) : smile string to be canonicalized

    Returns:
    mol (rdkit.Chem.rdchem.Mol) : RdKit mol object                          (None if invalid smile string smi)
    smi_canon (string)          : Canonicalized smile representation of smi (None if invalid smile string smi)
    conversion_successful (bool): True/False to indicate if conversion was  successful
    '''
    try:
        mol = MolFromSmiles(smi, sanitize=True)
        smi_canon = MolToSmiles(mol, isomericSmiles=False, canonical=True)
        return (mol, smi_canon, True)
    except:
        return (None, None, False)

def sampling_smi(target, population_size):
    mol = Chem.MolFromSmiles(target)
    if mol == None:
        raise Exception('Invalid starting structure encountered')
    randomized_smile_orderings = [randomize_smiles(mol) for _ in range(population_size)]
    return randomized_smile_orderings

def randomize_smiles(mol):
    '''Returns a random (dearomatized) SMILES given an rdkit mol object of a molecule.

    Parameters:
    mol (rdkit.Chem.rdchem.Mol) :  RdKit mol object (None if invalid smile string smi)

    Returns:
    mol (rdkit.Chem.rdchem.Mol) : RdKit mol object  (None if invalid smile string smi)
    '''
    if not mol:
        return None

    Chem.Kekulize(mol)

    return rdkit.Chem.MolToSmiles(mol, canonical=False, doRandom=True, isomericSmiles=False, kekuleSmiles=True)