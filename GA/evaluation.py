from GA import utils
import numpy as np
import multiprocessing

from SAS_calculator.sascorer import calculateScore
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity

def get_ECFP4(mol):
    return AllChem.GetMorganFingerprint(mol, 2)

def get_logP(mol):
    '''Calculate logP of a molecule

    Parameters:
    mol (rdkit.Chem.rdchem.Mol) : RdKit mol object, for which logP is to calculates

    Returns:
    float : logP of molecule (mol)
    '''
    return Descriptors.MolLogP(mol)

def get_fp_score(mol, target):
    fp_target = get_ECFP4(target)
    fp_mol = get_ECFP4(mol)
    score = TanimotoSimilarity(fp_mol, fp_target)
    return score

def calc_prop_logP(unseen_smile_ls, property_name, props_collect):
    '''Calculate logP for each molecule in unseen_smile_ls, and record results
       in locked dictionary props_collect
    '''
    for smile in unseen_smile_ls:
        mol, smi_canon, did_convert = utils.sanitize_smiles(smile)
        if did_convert:  # ensure valid smile
            props_collect[property_name][smile] = get_logP(mol)  # Add calculation
        else:
            raise Exception('Invalid smile encountered while atempting to calculate logP')


def calc_prop_SAS(unseen_smile_ls, property_name, props_collect):
    '''Calculate synthetic accesibility score for each molecule in unseen_smile_ls,
       results are recorded in locked dictionary props_collect
    '''
    for smile in unseen_smile_ls:
        mol, smi_canon, did_convert = utils.sanitize_smiles(smile)
        if did_convert:  # ensure valid smile
            props_collect[property_name][smile] = calculateScore(mol)
        else:
            raise Exception('Invalid smile encountered while atempting to calculate SAS')


def calc_prop_RingP(unseen_smile_ls, property_name, props_collect):
    '''Calculate Ring penalty for each molecule in unseen_smile_ls,
       results are recorded in locked dictionary props_collect
    '''
    for smile in unseen_smile_ls:
        mol, smi_canon, did_convert = utils.sanitize_smiles(smile)
        if did_convert:
            cycle_list = mol.GetRingInfo().AtomRings()
            if len(cycle_list) == 0:
                cycle_length = 0
            else:
                cycle_length = max([len(j) for j in cycle_list])
            if cycle_length <= 6:
                cycle_length = 0
            else:
                cycle_length = cycle_length - 6
            props_collect[property_name][smile] = cycle_length
        else:
            raise Exception('Invalid smile encountered while atempting to calculate Ring penalty')

def calc_prop_SIMIL(starting_smile, unseen_smile_ls, property_name, props_collect):
    '''Calculate logP for each molecule in unseen_smile_ls, and record results
       in locked dictionary props_collect
    '''
    target, _, _ = utils.sanitize_smiles(starting_smile)

    for smile in unseen_smile_ls:
        mol, smi_canon, did_convert = utils.sanitize_smiles(smile)
        if did_convert:                                                                # ensure valid smile
            props_collect[property_name][smile] = get_fp_score(mol, target) # Add calculation
        else:
            raise Exception('Invalid smile encountered while atempting to calculate SIMILARITY: ', smile)

# QED 추가하기!!!!!!!!!!!!!!!!!!!!
def create_parr_process(manager, chunks, property_name, target_smi):
    ''' Create parallel processes for calculation of properties
    '''
    # Assign data to each process
    process_collector = []
    collect_dictionaries = []

    for item in chunks:
        props_collect = manager.dict(lock=True)
        smiles_map_ = manager.dict(lock=True)
        props_collect[property_name] = smiles_map_
        collect_dictionaries.append(props_collect)
        if property_name == 'logP':
            process_collector.append(
                multiprocessing.Process(target=calc_prop_logP,
                                        args=(item, property_name, props_collect,)))

        if property_name == 'SAS':
            process_collector.append(
                multiprocessing.Process(target=calc_prop_SAS,
                                        args=(item, property_name, props_collect,)))

        if property_name == 'RingP':
            process_collector.append(
                multiprocessing.Process(target=calc_prop_RingP,
                                        args=(item, property_name, props_collect,)))

        if property_name == 'SIMIL':
            process_collector.append(
                multiprocessing.Process(target=calc_prop_SIMIL,
                                        args=(target_smi, item, property_name, props_collect, )))
    for item in process_collector:
        item.start()

    for item in process_collector:  # wait for all parallel processes to finish
        item.join()

    combined_dict = {}  # collect results from multiple processess
    for i, item in enumerate(collect_dictionaries):
        combined_dict.update(item[property_name])

    return combined_dict

def calc_smilarity(molecules_here, similarity_results):
    Similarity_calculated = []
    for smi in molecules_here:
        Similarity_calculated.append(similarity_results[smi])

    return Similarity_calculated

def calc_standardized_properties(molecules_here, logP_results, SAS_results, ringP_results,
                                     properties_calc_ls):
    ''' Obtain calculated properties of molecules in molecules_here, and standardize
    values base on properties of the Zinc Data set.
    '''
    logP_calculated = []
    SAS_calculated = []
    RingP_calculated = []

    for smi in molecules_here:
        if 'logP' in properties_calc_ls:
            logP_calculated.append(logP_results[smi])
        if 'SAS' in properties_calc_ls:
            SAS_calculated.append(SAS_results[smi])
        if 'RingP' in properties_calc_ls:
            RingP_calculated.append(ringP_results[smi])

    logP_calculated = np.array(logP_calculated)
    SAS_calculated = np.array(SAS_calculated)
    RingP_calculated = np.array(RingP_calculated)

    # Standardize logP based on zinc logP (mean: 2.4729421499641497 & std : 1.4157879815362406)
    logP_norm = (logP_calculated - 2.4729421499641497) / 1.4157879815362406
    logP_norm = logP_norm.reshape((logP_calculated.shape[0], 1))

    # Standardize SAS based on zinc SAS(mean: 3.0470797085649894    & std: 0.830643172314514)
    SAS_norm = (SAS_calculated - 3.0470797085649894) / 0.830643172314514
    SAS_norm = SAS_norm.reshape((SAS_calculated.shape[0], 1))

    # Standardiize RingP based on zinc RingP(mean: 0.038131530820234766 & std: 0.2240274735210179)
    RingP_norm = (RingP_calculated - 0.038131530820234766) / 0.2240274735210179
    RingP_norm = RingP_norm.reshape((RingP_calculated.shape[0], 1))

    return logP_calculated, SAS_calculated, RingP_calculated, logP_norm, SAS_norm, RingP_norm

def evaluation_properties(molcules, target_smi, properties_calc_ls, cpu_count, manager):
    molecules_here_unique = list(set(molcules))
    #molecules_here_unique = list(molcules)
    ratio = len(molecules_here_unique) / cpu_count
    chunks = utils.get_chunks(molecules_here_unique, cpu_count, ratio)
    chunks = [item for item in chunks if len(item) >= 1]

    logP_results, SAS_results, ringP_results, similarity_results = {}, {}, {}, {}

    if 'logP' in properties_calc_ls:
        logP_results = create_parr_process(manager, chunks, 'logP', target_smi)
    if 'SAS' in properties_calc_ls:
        SAS_results = create_parr_process(manager, chunks, 'SAS', target_smi)
    if 'RingP' in properties_calc_ls:
        ringP_results = create_parr_process(manager, chunks, 'RingP', target_smi)
    if 'SIMIL' in properties_calc_ls:
        similarity_results = create_parr_process(manager, chunks, 'SIMIL', target_smi)

    return logP_results, SAS_results, ringP_results, similarity_results, molecules_here_unique

def fitness_plogp(logP_results, SAS_results, ringP_results, properties_calc_ls, all_smi):
    logP_calculated, SAS_calculated, RingP_calculated, \
    logP_norm, SAS_norm, RingP_norm = calc_standardized_properties(all_smi, logP_results, SAS_results, ringP_results, properties_calc_ls)

    fitness = (logP_norm) - (SAS_norm) - (RingP_norm)
    return fitness

def fintenss_cplogp(logP_results, SAS_results, ringP_results, similarity_results, properties_calc_ls, all_smi, delta):
    logP_calculated, SAS_calculated, RingP_calculated, \
    logP_norm, SAS_norm, RingP_norm = calc_standardized_properties(all_smi, logP_results, SAS_results, ringP_results,
                                                                   properties_calc_ls)

    fitness = (logP_norm) - (SAS_norm) - (RingP_norm)

    Similarity_calculated = calc_smilarity(all_smi, similarity_results)
    Similarity_calculated = np.array([0 if x > delta else -10 ** 6 for x in Similarity_calculated])
    Similarity_calculated = Similarity_calculated.reshape((fitness.shape[0], 1))

    fitness = fitness + Similarity_calculated

    return fitness

def fitness_smilarity(similarity_results, all_smi, delta):
    Similarity_calculated = calc_smilarity(all_smi, similarity_results)
    fitness = np.array([0 if x > delta else -10 ** 6 for x in Similarity_calculated])

    return fitness

#def fitness_constrained_satisfaction(molcules):
#    similarity_results =
#    Similarity_calculated = obtain_only_smilarity(molcules, similarity_results)