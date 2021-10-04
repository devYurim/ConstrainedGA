import random
import numpy as np
import GA.crossover as co
import GA.mutation as mu

def apply_generation_cutoff(order, generation_size):
    ''' Return of a list of indices of molecules that are kept (high fitness)
        and a list of indices of molecules that are replaced   (low fitness)

    The cut-off is imposed using a Fermi-Function

    Parameters:
    order (list)          : list of molecule indices arranged in Decreasing order of fitness
    generation_size (int) : number of molecules in a generation

    Returns:
    to_replace (list): indices of molecules that will be replaced by random mutations of
                       molecules in list 'to_keep'
    to_keep    (list): indices of molecules that will be kept for the following generations
    '''
    # Get the probabilities that a molecule with a given fitness will be replaced
    # a fermi function is used to smoothen the transition
    positions = np.array(range(0, len(order))) - 0.2 * float(len(order))
    probabilities = 1.0 / (1.0 + np.exp(-0.02 * generation_size * positions / float(len(order))))

    #    import matplotlib.pyplot as plt
    #    plt.plot(positions, probabilities)
    #    plt.show()

    to_replace = []  # all molecules that are replaced
    to_keep = []  # all molecules that are kept
    for idx in range(0, len(order)):
        tmp = np.random.rand(1)
        if tmp < probabilities[idx]:
            to_replace.append(idx)
        else:
            to_keep.append(idx)

    return to_replace, to_keep

def reproduce_constrined(order, to_replace, to_keep, offspring_size, mutation_rate, constraint_rate, max_molecules_len=81):
    ''' Obtain the next generation of molecules. Bad molecules are replaced by
    mutations of good molecules

    Parameters:
    order (list)            : list of molecule indices arranged in Decreasing order of fitness
    to_replace (list)       : list of indices of molecules to be replaced by random mutations of better molecules
    to_keep (list)          : list of indices of molecules to be kept in following generation
    max_molecules_len (int) : length of largest molecule

    Returns:
    offsprings (list): next generation of mutated molecules as SMILES
    '''

    offspring = [order[idx] for idx in to_keep]

    for iter in range(offspring_size):  # smiles to replace (by better molecules)
        parent_a = random.choice(order)
        parent_b = random.choice(order)

        new_mol = co.crossover(parent_a, parent_b, max_molecules_len)  # 1 offspring generation
        if new_mol == None:
            if len(to_replace) != 0:
                new_mol = order[random.choice(to_replace)]
            else:
                new_mol = random.choice(order)

        new_mol = mu.mutate(new_mol, mutation_rate, constraint_rate, max_molecules_len)
        offspring.append(new_mol)

    return offspring