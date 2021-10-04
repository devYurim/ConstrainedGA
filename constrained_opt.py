import argparse
import numpy as np
import os
import multiprocessing
from multiprocessing import freeze_support
import GA.utils as utils
import GA.evaluation as evaluation
import GA.selection as selection
import GA.reproduce as reproduce
import json
import time
def print_logger(logger, generation, population_scores, gen_time, mol_sec):
    logger.debug(f'{generation} | '
                 f'max: {np.max(population_scores):.3f} | '
                 f'avg: {np.mean(population_scores):.3f} | '
                 f'min: {np.min(population_scores):.3f} | '
                 f'std: {np.std(population_scores):.3f} | '
                 f'sum: {np.sum(population_scores):.3f} | '
                 f'{gen_time:.2f} sec/gen | '
                 f'{mol_sec:.2f} mol/sec')

def constrained_satisfaction(similarity_results, all_smi, delta, population_size):
    fitness = evaluation.fitness_smilarity(similarity_results, all_smi, delta)
    all_smi, scores = selection.penalty_filter(all_smi, fitness, population_size)
    return all_smi, scores

def constrained_optimizer(logP_results, SAS_results, ringP_results, similarity_results, properties_calc_ls ,all_smi, delta, population_size):
    fitness = evaluation.fintenss_cplogp(logP_results, SAS_results, ringP_results, similarity_results, properties_calc_ls, all_smi, delta)
    all_smi, scores = selection.top_k(all_smi, fitness, population_size)

    return all_smi, scores

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--smiles_file', default='data/zinc_dearom.txt')
    parser.add_argument('--target_file', default='data/logp_800.txt')
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--population_size', type=int, default=100)
    parser.add_argument('--offspring_size', type=int, default=1000)
    parser.add_argument('--mutation_rate', type=float, default=0.5)
    parser.add_argument('--constraint_rate', type=float, default=0.1)
    parser.add_argument('--generations', type=int, default=20)
    parser.add_argument('--output_dir', type=str, default='./output/constrained/')
    parser.add_argument('--patience', type=int, default=5)
    parser.add_argument('--delta', type=float, default=0.6) # similarity threshold

    args = parser.parse_args()
    np.random.seed(args.seed)

    if args.output_dir is None:
        args.output_dir = os.path.dirname(os.path.realpath(__file__)) + "./output/constrained/"
    if not os.path.isdir(args.output_dir):
        os.mkdir(args.output_dir)
    with open(os.path.join(args.output_dir, 'constrained_params.json'), 'w') as jf:
        json.dump(vars(args), jf, sort_keys=True, indent=4)

    logger = utils.set_logger(os.path.join(args.output_dir, 'log.json'))
    manager = multiprocessing.Manager()
    lock = multiprocessing.Lock()
    cpu_count = multiprocessing.cpu_count()

    # load dataset (zinc)
    if args.smiles_file != None:
        all_smi = utils.read_dataset(args.smiles_file)
    else:
        all_smi =''
    target_smi = utils.read_dataset(args.target_file)
    simil_calc_ls = ['SIMIL']
    properties_calc_ls = ['logP', 'SAS', 'RingP', 'SIMIL']
    cnt = 1
    for target in target_smi:
        population_scores = []
        population_smi = None
        logger.debug(f'{target}')

        # init population
        if population_smi == None:
            print('selecting initial population...')
            logP_results, SAS_results, ringP_results, similarity_results, all_smi = evaluation.evaluation_properties(all_smi, target, simil_calc_ls, cpu_count, manager)
            population_smi, population_scores = constrained_satisfaction(similarity_results, all_smi, args.delta, args.population_size)

        if len(population_smi) == 0:
            print("Not Found")
            print('random sampling initial population...')
            all_smi = utils.sampling_smi(target, args.population_size)
            logP_results, SAS_results, ringP_results, similarity_results, all_smi = evaluation.evaluation_properties(
                all_smi, target, simil_calc_ls, cpu_count, manager)
            population_smi, population_scores = constrained_satisfaction(similarity_results, all_smi, args.delta,
                                                                         args.population_size)
        patience = 0  # threshold
        t0 = time.time()
        max_molecules = []
        for generation in range(args.generations):
            to_replace, to_keep = reproduce.apply_generation_cutoff(population_smi, len(population_smi))
            offsprings = reproduce.reproduce_constrined(population_smi, to_replace, to_keep, args.offspring_size,
                                                  args.mutation_rate, args.constraint_rate)

            gen_time = time.time() - t0
            mol_sec = args.offspring_size / gen_time
            t0 = time.time()

            #PHASE 1
            print('PHASE 1-------------------------------------')
            logP_results, SAS_results, ringP_results, similarity_results, smiles = evaluation.evaluation_properties(
                offsprings, target, simil_calc_ls, cpu_count, manager)
            smiles, scores = constrained_satisfaction(similarity_results, offsprings, args.delta, args.population_size)
            print("current population size :", len(smiles))

            # PHASE 2
            if len(smiles) != 0:
                print('PHASE 2-------------------------------------')
                logP_results, SAS_results, ringP_results, similarity_results, smiles = evaluation.evaluation_properties(
                    smiles, target, properties_calc_ls, cpu_count, manager)
                print(len(smiles))
                smiles, scores = constrained_optimizer(logP_results, SAS_results, ringP_results, similarity_results,
                                                       properties_calc_ls, smiles, args.delta, args.population_size)
                old_scores = population_scores
                population_smi, population_scores = selection.top_k(smiles, scores, args.population_size)
                print("current population size :", len(population_smi))

                if len(population_smi) != 0:
                    print("max molecules :", population_smi[0], ", score:", population_scores[0])
                    max_molecules.append([population_smi[0], population_scores[0]])  # max save
                else:
                    print("fail")
                    max_molecules.append(['fail', -1000000])

                if population_scores == old_scores:
                    patience += 1
                    print(f'Failed to progress: {patience}')
                    if args.patience == patience:
                        print(f'No more patience, bailing...')
                        print_logger(logger, generation, population_scores, gen_time, mol_sec)
                        break
                else:
                    patience = 0

                print_logger(logger, generation, population_scores, gen_time, mol_sec)

            else:
                print('break')
                print("fail")
                max_molecules.append(['fail', -1000000])
                break

        f = open('{}/{}_similarity{}_molecular(max).txt'.format(args.output_dir, cnt, args.delta), 'w')
        for i, smi_score in enumerate(max_molecules):
            f.write("{},{},{}\n".format(i, smi_score[0], smi_score[1].astype(float)[0]))
        f.close()
        cnt += 1
if __name__ == "__main__":
    freeze_support()
    main()