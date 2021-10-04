import random
import numpy as np
import selfies as sf
from GA import utils


def get_list_from_SELFIES(sel):
    '''
    Split a SELFIES into its symbols
    :param sel(string): SELFIES string
    :return (list) : SELFIES symbol list
    '''
    return list(sf.split_selfies(sel))

def get_symbols():
    '''
    alphabet of SELFIES symbols
    :return: (list) SELFIES symbols
    '''

    alphabet = sf.get_semantic_robust_alphabet()
    alphabet.update([
        '[#Br]', '[#C@@Hexpl]', '[#C@@expl]', '[#C@Hexpl]',
        '[#C@expl]', '[#C]', '[#Cl]', '[#F]', '[#Hexpl]', '[#I]',
        '[#NHexpl]', '[#N]', '[#O]', '[#P]', '[#S]', '[/Br]',
        '[/C@@Hexpl]', '[/C@@expl]', '[/C@Hexpl]', '[/C@expl]', '[/C]',
        '[/Cl]', '[/F]', '[/Hexpl]', '[/I]', '[/NHexpl]', '[/N]',
        '[/O]', '[/P]', '[/S]', '[=Br]', '[=C@@Hexpl]', '[=C@@expl]',
        '[=C@Hexpl]', '[=C@expl]', '[=C]', '[=Cl]', '[=F]', '[=Hexpl]',
        '[=I]', '[=NHexpl]', '[=N]', '[=O]', '[=P]', '[=S]', '[Br]',
        '[Branch1_1]', '[Branch1_2]', '[Branch1_3]', '[Branch2_1]',
        '[Branch2_2]', '[Branch2_3]', '[Branch3_1]', '[Branch3_2]',
        '[Branch3_3]', '[C@@Hexpl]', '[C@@expl]', '[C@Hexpl]',
        '[C@expl]', '[C]', '[Cl]', '[Expl#Ring1]', '[Expl=Ring1]',
        '[F]', '[Hexpl]', '[I]', '[NHexpl]', '[N]', '[O]', '[P]',
        '[Ring1]', '[Ring2]', '[Ring3]', '[S]', '[\\Br]',
        '[\\C@@Hexpl]', '[\\C@@expl]', '[\\C@Hexpl]', '[\\C@expl]',
        '[\\C]', '[\\Cl]', '[\\F]', '[\\Hexpl]', '[\\I]', '[\\NHexpl]',
        '[\\N]', '[\\O]', '[\\P]', '[\\S]', '[epsilon]', '[nop]'
    ])

    return list(alphabet)

def mutate(smi, mutation_rate, constraint_rate, max_molecules_len):
    '''

    :param
    sel:
    mutation_rate:
    :return:
    '''
    symbols = get_symbols()
    '''symbols = ['[Branch1_1]', '[Branch1_2]', '[Branch1_3]', '[epsilon]', '[Ring1]', '[Ring2]', '[Branch2_1]',
               '[Branch2_2]', '[Branch2_3]', '[F]', '[O]', '[=O]', '[N]', '[=N]', '[#N]', '[C]', '[=C]', '[#C]', '[S]',
               '[=S]', '[C][=C][C][=C][C][=C][Ring1][Branch1_1]']'''
    prob = random.random()
    if prob > mutation_rate:
        return smi
    else:
        sel = sf.encoder(smi)
        sel = get_list_from_SELFIES(sel)
        choice_ls = [1, 2, 3]  # 1=Replace; 2=Insert; 3=Delete
        random_choice = np.random.choice(choice_ls, 1)[0]

        if constraint_rate != 1.0:
            mutation_range_front = [i for i in range(0, int(len(sel) * constraint_rate))]
            mutation_range_back = [i for i in range(len(sel), len(sel) - int(len(sel) * constraint_rate), -1)]
            mutation_range = mutation_range_front + mutation_range_back
        else:
            mutation_range = [i for i in range(len(sel))]
        if len(mutation_range) == 0:
            return smi

        if random_choice == 1:  # insert
            mutation_range.append(len(sel) + 1)
            mutated_gene = np.random.choice(mutation_range)
            random_symbol = np.random.choice(symbols, size=1)[0]
            offsp_sel = sel[:mutated_gene] + [random_symbol] + sel[mutated_gene:]

        elif random_choice == 2: # replace
            mutated_gene = np.random.choice(mutation_range)
            random_symbol = np.random.choice(symbols, size=1)[0]
            if mutated_gene == 0:
                offsp_sel = [random_symbol] + sel[mutated_gene + 1:]
            else:
                offsp_sel = sel[:mutated_gene] + [random_symbol] + sel[mutated_gene + 1:]

        elif random_choice == 3: # delete
            mutated_gene = np.random.choice(mutation_range)

            if mutated_gene == 0:
                offsp_sel = sel[mutated_gene+1:]
            else:
                offsp_sel = sel[:mutated_gene] + sel[mutated_gene+1:]
        else:
            raise Exception('Invalid Operation trying to be performed')

        sel = sf.decoder(''.join(sel))
        offsp_smi = sf.decoder(''.join(offsp_sel))
        offsp_mol, offsp_smi, valid = utils.sanitize_smiles(offsp_smi)
        if valid and len(offsp_smi) < max_molecules_len and offsp_smi != "":
            return offsp_smi
        else:
            f = open('selfie_failure_cases.txt','a+')
            f.write('Tried to mutate SELFIE: ' + str(sel)+'To Obtain: '+str(offsp_sel)+'\n')
            f.close()

    return smi
