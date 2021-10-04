
def top_k(smiles, scores, k):
    scored_smiles = list(zip(scores, smiles))
    scored_smiles = sorted(scored_smiles, key=lambda x: x[0], reverse=True)
    #scored_smiles_filtering = []

    #for s in filter(lambda x: x[0] > -100, scored_smiles):
    #    scored_smiles_filtering.append(s)

    return [smile for score, smile in scored_smiles][:k], [score for score, smile in scored_smiles][:k]

def penalty_filter(smiles, scores, population_size):
    scored_smiles = list(zip(scores, smiles))
    scored_smiles = sorted(scored_smiles, key=lambda x: x[0], reverse=True)
    scored_smiles_filtering = []
    for s in filter(lambda x: x[0] > -100, scored_smiles):
        scored_smiles_filtering.append(s)
    return [smile for score, smile in scored_smiles_filtering][:population_size], [score for score, smile in scored_smiles_filtering][:population_size]