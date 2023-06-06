from joblib import Parallel, delayed
import multiprocessing
from itertools import islice, combinations
import os
import math
import numpy as np
from scipy.spatial.distance import pdist
import distance
import random


def chunks(data, num_seqs):
    # Create batches
    it = iter(data)
    for i in range(0, len(data), num_seqs):
        yield {k: data[k] for k in islice(it, num_seqs)}


def chunks_list(data, chunk_size):
    return [data[i:i + chunk_size] for i in range(0, len(data), chunk_size)]

def rnafold(file_name):
    cmd = f'RNAfold -p -d2 --noLP --noPS --noDP < {file_name} > {file_name[:-3]}_fold.out'
    os.system(cmd)


def rnafold_in_parallel(iso2, output_file_stem, ac):
    num_cores = multiprocessing.cpu_count() - 1
    chunk_size = math.ceil(len(iso2.trnas)/num_cores)
    trna_batches = chunks(iso2.trnas, chunk_size)

    file_names = []
    for i, batch in enumerate(trna_batches, 1):
        with open(f'{output_file_stem}_{ac}_{i}.fa', 'w+') as f:
            for name, trna in batch.items():
                f.write(f'>{name}\n{trna.seq[ac]}\n')
        file_names.append(f'{output_file_stem}_{ac}_{i}.fa')

    Parallel(n_jobs=num_cores)(delayed(rnafold)(filename) for filename in file_names)

    outfile_names = [f'{file_name[:-3]}_fold.out' for file_name in file_names]
    with open(f'{output_file_stem}_{ac}_complete_fold.out', 'w') as outfile:
        for f_name in outfile_names:
            with open(f_name) as infile:
                for line in infile:
                    outfile.write(line)


def memo_dist(u, v, memo):
    if frozenset((u[0], v[0])) in memo.keys():
        return memo[frozenset((u[0], v[0]))]
    else:
        dist_ = distance.levenshtein(u[0], v[0])
        memo[frozenset((u[0], v[0]))] = dist_
        return dist_


def max_dist_parallel_memo(part_list, num_seqs, ac, type='part'):
    num_cores = multiprocessing.cpu_count() - 1
    if type == 'part':
        seqs = [[part.seq] for part in part_list]
    elif type == 'tRNA':
        seqs = [[tRNA.seq[ac]] for tRNA in part_list]
    c = [list(x) for x in combinations(seqs, num_seqs)]
    batches = chunks_list(c, math.ceil(len(c) / num_cores))

    def max_dist_memo(batch):
        memo = {}
        distances = []
        for i in batch:
            distances.append(np.min(pdist(i, lambda u, v: memo_dist(u, v, memo))))
        max_dist = max(distances)
        inds = [i for i, x in enumerate(distances) if x == max_dist]
        rows = [batch[ind] for ind in inds]
        return (max_dist, rows)

    processed_list = Parallel(n_jobs=num_cores, backend='loky')(delayed(max_dist_memo)(batch) for batch in batches)
    processed_list = sorted(processed_list, key=lambda tup: tup[0])
    final_trnas = random.choice(processed_list[-1][1])
    dist = processed_list[-1][0]
    final_trnas = [seq for seq_list in final_trnas for seq in seq_list]
    return final_trnas, dist

