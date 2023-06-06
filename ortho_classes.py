import pandas as pd
import numpy as np
import re
from collections import Counter, namedtuple
from sklearn.cluster import AffinityPropagation
import distance
from itertools import product, chain, islice
import time
import ast
from matplotlib import pyplot as plt
import random
import sys
from parallel import max_dist_parallel_memo

# ID elements in the form AlaRS: 2, 3, 5...
id_df = pd.read_excel('identity_elements.xlsx', index_col=0)

# Remove anticodon identity elements - these are changing anyway
id_dict = {index[:3]: {int(val): set() for val in value.strip().split(', ') if int(val) not in {34, 35, 36}}
           for index, value in zip(id_df.index, id_df.iloc[:, 0])}

# df of e. coli tRNAs aligned to 1, 2, 3, 4.... - from Daniele paper
ecoli_df = pd.read_excel('ecoli_tRNAs.xlsx', index_col=0, header=1)
col_dict = {f'Unnamed: {i}': f'VL{j}'
            for i, j in zip(range(51, 68), range(1, 18))}
ecoli_df = ecoli_df.rename(columns={'Unnamed: 1': 'AC', 'Unnamed: 9': '7a', 'Unnamed: 88': '65a'})
ecoli_df = ecoli_df.rename(columns=col_dict)

# Fill id_dict with bases at ID positions
for aa_, ids in id_dict.items():
    for pos_, bases in ids.items():
        bases.update(list(ecoli_df.loc[aa_, pos_]))


def synth_clean(synth_file):
    synth_df = pd.read_excel(synth_file, header=None)
    synth_df = synth_df.iloc[:, :4]
    synth_df.columns = ['synth', 'syn_id', 'trna_id', 'gen_id']
    synth_df.trna_id = synth_df.trna_id.apply(lambda x: x.replace('>', ''))
    return synth_df


def open_dict(dictofdicts):
    new = {}
    for key, val in dictofdicts.items():
        if isinstance(val, dict):
            for k, v in val.items():
                new[k] = v
        else:
            new[key] = val
    return new


# RangeDict data structure can take ranges, or tuples of ranges as keys, and any type as values
# Required to return parts from ID bases
class RangeDict(dict):
    def __getitem__(self, item):
        if not isinstance(item, range):
            for i, key in enumerate(self.keys()):
                if isinstance(key, tuple):
                    if any([item in r for r in key]):
                        return list(self.values())[i]
                if item in key:
                    return self[key]
            raise KeyError(item)
        else:
            return super().__getitem__(item)


base_to_part = RangeDict({range(1, 8): 'tRNA1-7*', range(8, 10): 'tRNA8-9*',
                          range(10, 14): 'tRNA10-13*', range(14, 22): 'tRNA14-21*',
                          range(22, 26): 'tRNA22-25*', range(26, 27): 'tRNA26*',
                          range(27, 32): 'tRNA27-31*', range(32, 39): 'tRNA32-38*',
                          range(39, 44): 'tRNA39-43*', range(44, 49): 'tRNA44-48*',
                          range(49, 54): 'tRNA49-53*', range(54, 61): 'tRNA54-60*',
                          range(61, 66): 'tRNA61-65*', range(66, 73): 'tRNA66-72*',
                          range(73, 77): 'tRNA73-76*'})

part_to_range = {val: key for key, val in base_to_part.items()}

# RangeDict structure to return part name based on single canonical base number
base_to_part_2 = RangeDict({(range(1, 8), range(66, 73)): 'tRNA1-7_66-72*', range(8, 10): 'tRNA8-9*',
                            (range(10, 14), range(22, 26)): 'tRNA10-13_22-25*',
                            (range(14, 22), range(54, 61)): 'tRNA14-21_54-60*',
                            (range(26, 27), range(44, 49)): 'tRNA26_44-48*',
                            (range(27, 32), range(39, 44)): 'tRNA27-31_39-43*',
                            range(32, 39): 'tRNA32-38*', (range(49, 54), range(61, 66)): 'tRNA49-53_61-65*',
                            range(73, 77): 'tRNA73-76*'})

# And reverse to return the ranges from the part name
part_to_range_2 = {val: key for key, val in base_to_part_2.items()}

base_comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
# base_comp = {'A': ['T'], 'T': ['G', 'A'], 'C': ['G'], 'G': ['C', 'T']}

part_order = ['tRNA1-7*', 'tRNA8-9*', 'tRNA10-13*', 'tRNA14-21*', 'tRNA22-25*', 'tRNA26*',
              'tRNA27-31*', 'tRNA32-38*', 'tRNA39-43*', 'tRNA44-48*', 'tRNA49-53*', 'tRNA54-60*',
              'tRNA61-65*', 'tRNA66-72*', 'tRNA73-76*']

# trna_pattern = re.compile(
#     '^\\({5,8}\\.{1,3}\\({4}\\.{5,}\\){3,4}\\.*\\({4,9}\\.{7}\\){4,9}.*\\.\\({5}\\.{2,}\\){5}\\){5,8}\\.{3,}$')
# gly04_pattern = re.compile(
#     '^\\({5,8}\\.{1,3}\\({4}\\.{5,}\\){3,4}\\.*\\({4,9}\\.{3,7}\\){4,9}.*\\.\\({5}\\.{2,}\\){5}\\){5,8}\\.{3,}$')
# ile_pattern = re.compile(
#     '^\\({5,8}\\.{1,3}\\({3,4}\\.{5,}\\){3,4}\\.*\\({4,9}\\.{7}\\){4,9}.*\\.\\({5}\\.{2,}\\){5}\\){5,8}\\.{3,}$')
# trna_pattern_arg08 = re.compile(
#     '^\\({5,8}\\.{1,3}\\({3,4}\\.{2,}\\){3,4}\\.*\\({4,9}\\.{5,7}\\){4,9}.*\\.\\({5}\\.{2,}\\){5}\\){5,8}\\.{3,}$')

allowed = set('A' + 'C' + 'T' + 'G' + '_')


def check(test_str):
    return set(test_str) <= allowed


class Synthetase2(object):
    """Synthetase2 class contains information about the synthetases defined by the user.
    In our case, this is the synthetases we have on hand and want to use in our screens.
    Class contains information about the synthetase (database ID etc.) and importantly which parts
    e.g. tRNA8-9*, contain ID elements as defined by Daniele (which is also a user input above - scalable)
    On initialisation, class instance iterates through these ID parts and stores the sequences in the tRNAs
    corresponding to these parts. e.g. if an ID part is tRNA8-9*, then self.id_seqs will contain the sequence
    found at the corresponding tRNA - "TA\""""

    def __init__(self, name, syn_id, trna_id, gen_id, aa, id_parts, iso, huge_df):

        self.name = name
        self.syn_id = syn_id
        self.trna_id = trna_id
        self.gen_id = gen_id
        self.aa = aa
        self.id_parts = id_parts
        self.iso = iso
        self.huge_df = huge_df
        self.non_id_parts = [part_type for part_type in base_to_part_2.values() if part_type not in self.id_parts]
        self.id_seqs = {}
        self.non_id_seqs = {}

        for part in base_to_part_2.values():
            if part in self.id_parts:
                try:
                    seq = self.huge_df[self.huge_df.seq_id == self.trna_id].iloc[0][part]
                    if part != 'tRNA14-21_54-60*':
                        align = seq
                    else:
                        align = self.huge_df[self.huge_df.seq_id == self.trna_id].iloc[0]['tRNA14-21_54-60* aligned']
                    if isinstance(seq, str):
                        self.id_seqs[part] = Part2(seq, part, self.aa, self.trna_id, align, self.iso)
                    else:
                        raise Exception('ID sequences in databases are not present - check database')

                except IndexError:
                    pass
            else:
                try:
                    seq = self.huge_df[self.huge_df.seq_id == self.trna_id].iloc[0][part]
                    if part != 'tRNA14-21_54-60*':
                        align = seq
                    else:
                        align = self.huge_df[self.huge_df.seq_id == self.trna_id].iloc[0]['tRNA14-21_54-60* aligned']
                    if isinstance(align, str) and check(align):
                        self.non_id_seqs[part] = Part2(seq, part, self.aa, self.trna_id, align, self.iso)

                except IndexError:
                    pass


###################################################################################################################

class Isoacceptor2(object):
    """Main class for pipeline.
    Isoacceptor2 class takes in the synthetase dataframe defined above, the ID dictionary, and the type of
    isoacceptor we want to find tRNAs for e.g. "Arg".
    Class stores all information in custom data structures, in order to get away from slow pandas functions.
    Class contains methods for clustering by part type, and to generate chimeras.
    Will want to have output dataframe that can be analysed by RNAFold. Somehow, will want RNAFold output
    to be handled by this class, as will want to cluster once again to cover largest sequence diversity of
    designed tRNAs."""

    def __init__(self, synth_df, id_dict, aa, huge_df, ac='CTA', comp_arm_strict=True, v_loop_length=9,
                 id_part_change=None, num_iter=1):

        self.num_iter = num_iter

        # id_dict for instance just takes ID elements for defined isoacceptor class
        self.exemplar_parts = None
        self.id_dict = id_dict[aa].keys()
        if id_part_change:
            self.id_part_change = id_part_change
        else:
            self.id_part_change = []

        # id_parts to pass to Synthetase2 instance are parts containing ID bases + acceptor stem
        self.id_parts = {base_to_part_2[base] for base in self.id_dict}
        self.id_parts.update(['tRNA73-76*'])  # self.id_parts.update(['tRNA1-7_66-72*', 'tRNA73-76*'])
        self.id_parts = {part for part in self.id_parts if part not in self.id_part_change}

        # List of positions within id_parts that do not change
        self.non_negotiable_positions = [part_to_range_2[part] for part in self.id_parts]
        self.non_negotiable_positions = [[range_] if isinstance(range_, range) else range_
                                         for range_ in self.non_negotiable_positions]
        for i in range(2):
            self.non_negotiable_positions = list(chain.from_iterable(self.non_negotiable_positions))

        self.non_part_ids = {}

        for num in self.id_dict:
            if num in self.non_negotiable_positions:
                pass
            else:
                if base_to_part_2[num] in self.non_part_ids.keys():
                    self.non_part_ids[base_to_part_2[num]].append(num)
                else:
                    self.non_part_ids[base_to_part_2[num]] = [num]
                self.non_negotiable_positions.append(num)

        # id dict for id positions of other isoacceptors outside of the id positions for this isoacceptor
        # This is to prevent unchangeably high id scores which are constrained by the synthetase chosen
        self.non_id_id_dict = {aa_: {pos: bases for pos, bases in aas_ids.items()
                                     if pos not in self.non_negotiable_positions}
                               for aa_, aas_ids in id_dict.items()}

        self.synth_df = synth_df
        self.aa = aa
        self.ac = ac
        self.huge_df = huge_df
        self.trnas = dict()
        self.used_parts = {}

        # Holds all tRNAs from each iteration
        self.iter_trnas = {}

        self.comp_arm_strict = comp_arm_strict

        # Initialise list of Synthetase2 instances using user input synth_df
        self.synths = [Synthetase2(name, syn_id, trna_id, gen_id, self.aa, self.id_parts, self, self.huge_df)
                       for name, syn_id, trna_id, gen_id in
                       zip(self.synth_df.synth, self.synth_df.syn_id, self.synth_df.trna_id, self.synth_df.gen_id)]

        # Need to retain tRNA_id info - could have list of ids for every unique part sequence
        # Other solution would be to have another tRNA class - could store info here then read back, but may need
        # some work to make this play nice with Part class - looking like will need tRNA class, but ngl this is
        # a lot more work to change around

        # part_tuple is simple namedtuple data structure with a seq and align - only time this is different
        # is for D-loop, which has a sequence and an aligned sequence - just easier to make all parts have both
        part_tuple = namedtuple('Part_tuple', ['seq', 'align'])
        self.all_parts = {}

        # Cut down dataframe to just sequences within isoacceptor class
        amino_df = self.huge_df[self.huge_df['Amino Acid'] == self.aa]

        for part_type in part_to_range_2.keys():
            if part_type != 'tRNA14-21_54-60*':
                self.all_parts[part_type] = [part_tuple(seq, seq) for seq in
                                             amino_df[part_type].unique()
                                             if isinstance(seq, str)]
            else:
                # Might be quicker to do this with a zip command?
                # Although that would require subsetting dataframe for unique sequences first, then iterating,
                # so maybe not - point for pipeline optimisation in future perhaps
                self.all_parts[part_type] = [part_tuple(row['tRNA14-21_54-60*'], row['tRNA14-21_54-60* aligned'])
                                             for index, row in
                                             amino_df.drop_duplicates(subset=['tRNA14-21_54-60*']).iterrows()
                                             if isinstance(row['tRNA14-21_54-60*'], str)]

        # Filters out sequences with other characters e.g. 'R', 'M', 'N'
        allowed = set('A' + 'C' + 'T' + 'G' + '_')

        def check(test_str):
            return set(test_str) <= allowed

        # all_parts contains every part in the isoacceptor class, initialising Part2 instances
        # Will eventually input tRNA_id - or list of tRNA_ids since many contain the same part sequence
        # Actually not very necessary
        self.all_parts = {part_type: [Part2(part_t.seq, part_type, self.aa, 'TODO', part_t.align, self)
                                      for part_t in part_list
                                      if check(part_t.seq)
                                      and isinstance(part_t.align, str)]
                          for part_type, part_list in self.all_parts.items()}

        # Removes ASLs that close up
        self.all_parts['tRNA32-38*'] = [part for part in self.all_parts['tRNA32-38*']
                                        if (part.seq[0] + part.seq[-1] not in ['CG', 'GC'])]

        def comp_checker(part):
            """If comp_arm_strict True, then check that each paired region is complementary"""
            if part.region not in ['tRNA10-13_22-25*', 'tRNA27-31_39-43*', 'tRNA49-53_61-65*']:
                return True
            else:
                part1, part2 = part.seq.split('_')
                # part1_rev = ''.join([base_comp[c] for c in part1])[::-1]
                # Ile does not need to be exact for D-arm
                if self.aa.title() == 'Ile' and part.region == 'tRNA10-13_22-25*':
                    if sum([i in base_comp[j] for i, j in zip(part1, part2[::-1])]) >= 3:
                        return True
                    else:
                        return False
                else:
                    return all([i in base_comp[j] for i, j in zip(part1, part2[::-1])])
                    # if part1_rev == part2:
                    #     return True
                    # else:
                    #     return False

        # Extra step to remove any base paired parts that are not reverse complements (so max 256 seq)
        if self.comp_arm_strict:
            self.all_parts = {part_type: [part for part in part_list if comp_checker(part)]
                              for part_type, part_list in self.all_parts.items()}

        self.all_parts['tRNA26_44-48*'] = [part for part in self.all_parts['tRNA26_44-48*']
                                           if len(part.seq) <= (v_loop_length + 2)]

        self.all_parts['tRNA8-9*'] = [part for part in self.all_parts['tRNA8-9*'] if 'T' in part.seq]

        # Sort list of parts by scores, so clustering can take top scorers
        for part_list in self.all_parts.values():
            part_list.sort()

        self.folded = {}
        self.low_part_types = []
        self.mid_part_types = []

    #######################

    def cluster_parts(self, sample_size, synth_name, clust_id_parts=False, display=False, log_file=None, iteration=1,
                      automatic=False, clust_size_min=None, ham=False):

        """This method uses Levenshtein distance between sequences and affinity propagation to cluster parts.
        Clustering could be optimised. For maximum sample size, parameters have been adjusted to:
        damping=0.9 (parameter for extent to which the current value is maintained relative to incoming values)
        and max_iter=1000 (number of iterations clustering will run for without converging).
        With these parameters, have been able to push algorithm from n=377 to much greater than 2500,
        however for parts it is not clear if I need to push it this far. Default damping is 0.5.
        sample_size argument determines how many of the top parts to take and cluster.
        clust_id_parts takes True or False, and decides whether to also cluster parts containing identity elements
        for the chosen isoacceptor class - since these are kept constant for each synthetase.
        display is for me to show people how cool the clustering output is (will print results for all to see)."""

        now = time.time()
        print('Clustering Parts...')
        clust_dict = self.all_parts
        self.exemplar_parts = {}
        if not clust_size_min:
            clust_size_min = 15

        # Remove ID parts from clust dict
        if not clust_id_parts:
            clust_dict = {part_type: parts for part_type, parts in clust_dict.items()
                          if part_type not in self.id_parts}

        # Make all low numbered part types exemplars
        if iteration == 1:
            for part_type, parts in clust_dict.items():
                if ham:
                    for i, part in enumerate(parts):
                        part.cluster_id = i
                        part.exemplar = 'all'
                    self.low_part_types.append(part_type)
                elif len(parts) < 15:
                    self.low_part_types.append(part_type)
                    for i, part in enumerate(parts):
                        part.cluster_id = i
                        part.exemplar = 'all'
                elif len(parts) < self.num_iter * 3:
                    self.mid_part_types.append(part_type)

        print({part_type: len(part_list) for part_type, part_list in clust_dict.items()})

        # Take top parts
        clust_dict = {part_type: [part for part in parts[:sample_size]]
                      for part_type, parts in clust_dict.items()
                      if len(parts) >= clust_size_min}

        if log_file:
            with open(log_file, 'a') as f:
                f.write('Cluster Dict: \n' +
                        str({part_type: len(part_list) for part_type, part_list in clust_dict.items()}) + '\n')

        for part_type, parts in clust_dict.items():
            # turn sample into numpy array - not used for clustering, but to find parts at end
            sample = np.asarray([part.seq for part in parts])
            # Form square zero matrix
            m = np.zeros((len(parts), len(parts)))
            # Return the indices for the lower-triangle of matrix, with a k=-1 diagonal offset
            tril_idx_rows, tril_idx_cols = np.tril_indices_from(m, k=-1)
            # Fill these indices with the levenshtein distance between the part sequences
            m[(tril_idx_rows, tril_idx_cols)] = [distance.levenshtein(sample[i], sample[j]) for i, j in
                                                 zip(tril_idx_rows, tril_idx_cols)]
            # Form full matrix by adding to transpose
            m = (m + m.T) * -1

            damping = 0.6
            while True:
                try:
                    affprop = AffinityPropagation(affinity="precomputed", damping=damping, max_iter=1000,
                                                  random_state=None)
                    affprop.fit(m)
                    for cluster_id in np.unique(affprop.labels_):
                        # identify exemplar for each cluster (found at centre)
                        parts[affprop.cluster_centers_indices_[cluster_id]].exemplar = iteration
                        for element_list in np.nonzero(affprop.labels_ == cluster_id):
                            for element in element_list:
                                # assign cluster_id to each part in all_parts
                                parts[element].cluster_id = cluster_id
                        if display:
                            exemplar = sample[affprop.cluster_centers_indices_[cluster_id]]
                            cluster = sample[np.nonzero(affprop.labels_ == cluster_id)]
                            cluster_str = ", ".join(cluster)
                            print(f"{part_type}: ID - {cluster_id}: *{exemplar}* {cluster_str}")
                    break
                except IndexError:
                    print(part_type, damping)
                    for part in parts:
                        part.exemplar = False
                    damping += 0.1

            # iterate through clusters

        # self.exemplar_parts = {part_type: [part for part in part_list if part.exemplar in ['all', iteration]]
        #                        for part_type, part_list in self.all_parts.items()
        #                        if part_type not in self.id_parts}

        for part_type, part_list in self.all_parts.items():
            if part_type not in self.id_parts:  # and part_type not in self.exemplar_parts.keys():
                self.exemplar_parts[part_type] = [part for part in part_list if part.exemplar in ['all', iteration]]

        for part_type, part_list in self.exemplar_parts.items():
            if part_type in self.used_parts.keys():
                self.used_parts[part_type] += [part for part in part_list if part.exemplar == iteration]
            else:
                self.used_parts[part_type] = [part for part in part_list if part.exemplar == iteration]
            if part_type in self.mid_part_types:
                self.used_parts[part_type] = []

        est = np.prod([len(part_list) for part_list in self.exemplar_parts.values()])
        print(f'Estimated Chimeras: {est}')

        if not automatic:
            while True:
                size_ans = input('Continue? (y) or (n):\n')
                if size_ans.lower() == 'y':
                    break
                elif size_ans.lower() == 'n':
                    print(f'Consider changing sample size for clustering.\n Current sample size is {sample_size}\n')
                    while True:
                        size_ans_ = input("Enter new sample size or 'q' to quit\n> ")
                        if size_ans_.lower() == 'q':
                            print('Thank you for trying Chi-T!')
                            sys.exit()
                        else:
                            try:
                                new_size = int(size_ans_)
                                self.cluster_parts(new_size, synth_name)
                            except ValueError:
                                print('Inappropriate value!')
                                continue

        print(f'Parts Clustered!...Time elapsed: {time.time() - now}')

    ###########################

    def chimera(self, synth_name, length_filt=79, log_file=None, iteration=1):
        """Generates all possible chimeras from ID parts and exemplar parts"""
        count = 0

        now = time.time()
        print(f'Choosing exemplars...Time Elapsed: {time.time() - now}')
        # print([synth.name for synth in self.synths])
        for part_type, id_seq in [synth for synth in self.synths if synth.name == synth_name][0].id_seqs.items():
            self.exemplar_parts[part_type] = [id_seq]

        for part_type, non_id_seq in [synth for synth in self.synths
                                      if synth.name == synth_name][0].non_id_seqs.items():
            if non_id_seq.seq not in [part.seq for part in self.exemplar_parts[part_type]]:
                self.exemplar_parts[part_type].append(non_id_seq)

        # print(f'ID seqs: {[synth for synth in self.synths if synth.name == synth_name][0].id_seqs.items()}')
        # print(f'non ID seqs: {[synth for synth in self.synths if synth.name == synth_name][0].non_id_seqs.items()}')
        #
        # pprint.pprint(self.exemplar_parts)
        # sys.exit()

        if log_file:
            with open(log_file, 'a') as f:
                f.write('Exemplar Parts: \n')
                f.write(str({part_type: len(part_list) for part_type, part_list in self.exemplar_parts.items()}) + '\n')

        print(f'Mixing parts...Time Elapsed: {time.time() - now}')
        listoflistofparts = [[{part_type: part} for part in part_list] for part_type, part_list in
                             self.exemplar_parts.items()]
        self.trnas = list(product(*listoflistofparts))
        self.trnas = [{part_type: (part.sub_parts if part.sub_parts else part)
                       for part_dict in trna for part_type, part in part_dict.items()
                       }
                      for trna in self.trnas]
        self.trnas = [open_dict(trna) for trna in self.trnas]

        print(f'Joining parts...Time Elapsed: {time.time() - now}')
        self.trnas_ = {}

        for i, trna_ in enumerate(self.trnas):
            self.trnas_.update({f'{synth_name}_iter{iteration}_seq{i}': tRNA(trna_, self.ac)})
            count += 1
            if count % 400000 == 0:
                print(f'{count} chimeras made...Time Elapsed: {time.time() - now}')
        self.trnas = self.trnas_
        if length_filt:
            self.trnas = {name: trna for name, trna in self.trnas.items() if len(trna.seq[self.ac]) < length_filt}
        if log_file:
            with open(log_file, 'a') as f:
                f.write(f'Unfiltered Chimeras: {len(self.trnas_)}\n')
                if length_filt:
                    f.write(f'Length (<{length_filt}) Filtered Chimeras: {len(self.trnas)}\n')

        print(f'{len(self.trnas)} Chimeras Made!...Time Elapsed: {time.time() - now}')

    ##############################

    def cervettini_filter(self, output_dir, synth_name, iteration=1,
                          start_stringency=0.5, min_stringency=0, target=1500000, step_size=0.05, log_file=None):
        """Applies the scoring from Cervettini et al., 2020 then filters.
        Filtering applied iteratively until target sequence number reached.
        start_stringency determines initial filtering threshold for max score across any isoacceptor.
        min_stringency is the lowest threshold that will be applied if target not reached.
        target is the number of trnas below which the filtering must return.
        step_size is the decrease in threshold applied each iteration"""
        now = time.time()
        print(f'Scoring tRNAs...Time Elapsed: {time.time() - now}')
        for trna in self.trnas.values():
            if not trna.cer_score:
                trna.trna_cer_scorer()
        cscores = [max(trna.cer_score.values()) for trna in self.trnas.values()]
        plt.hist(cscores, bins=[i for i in np.arange(-1, 1, 0.05)])
        plt.ylabel('Count')
        plt.xlabel('Cervettini Score')
        plt.savefig(output_dir + f'/plots/cscores_{synth_name}_iter{iteration}.pdf')
        plt.clf()
        print(f'Filtering tRNAs...Time Elapsed: {time.time() - now}')
        stringency = start_stringency
        step_size = round(step_size, 4)
        log_string = 'Cervettini Filtering: \n'
        for i in range(round((start_stringency - min_stringency) / step_size) + 1):
            if len(self.trnas) <= target:
                break

            self.trnas = {name: trna for name, trna in self.trnas.items() if max(trna.cer_score.values()) <= stringency}
            log_string += f'Threshold: {stringency} tRNAs remaining: {len(self.trnas)}\n'
            print(f'Threshold: {stringency} tRNAs remaining: {len(self.trnas)}...Time Elapsed: {time.time() - now}')
            stringency = round(stringency - step_size, 4)

        if log_file:
            with open(log_file, 'a') as f:
                f.write(log_string)
        print(f'tRNAs Filtered!...tRNAs remaining: {len(self.trnas)}...Time elapsed: {time.time() - now}')

    ##############################

    def store_trnas(self, filename):
        """Stores self.trnas to a csv in a form that can retrieve the important information"""

        now = time.time()
        print(f'{len(self.trnas)} being stored at {filename}')
        dict_for_df = {seq_name: [trna.seq, {part_type: (part.seq if isinstance(part, Part2) else part)
                                             for part_type, part in trna.part_dict.items()},
                                  trna.cer_score, {ac: struct for ac, struct in trna.struct.items()},
                                  {ac: div for ac, div in trna.div.items()},
                                  {ac: freq for ac, freq in trna.freq.items()}]
                       for seq_name, trna in self.trnas.items()}
        if dict_for_df:
            trna_df = pd.DataFrame.from_dict(dict_for_df, orient='index')
            trna_df = trna_df.reset_index()
            trna_df.columns = ['name', 'seq', 'part_dict', 'cer_score', 'struct', 'div', 'freq']
            trna_df.to_csv(filename, index=False)
            print(f'tRNAs Stored!...Time elapsed: {time.time() - now}')
        else:
            print('No tRNAs to store!')

    def retrieve_trnas(self, filename):
        """Retrieves said tRNA information from .csv"""

        now = time.time()
        print(f'Retrieving tRNAs fom {filename}')
        trna_df = pd.read_csv(filename)
        trna_df = trna_df.set_index('name')
        trna_dict = trna_df.to_dict('index')
        try:
            self.trnas = {name: tRNA(ast.literal_eval(trna_info['part_dict']), self.ac,
                                     cer_score=ast.literal_eval(trna_info['cer_score']),
                                     struct=ast.literal_eval(trna_info['struct']),
                                     div=ast.literal_eval(trna_info['div']),
                                     freq=ast.literal_eval(trna_info['freq']),
                                     seq=ast.literal_eval(trna_info['seq']))
                          for name, trna_info in trna_dict.items()}
        except:
            self.trnas = {name: tRNA(ast.literal_eval(trna_info['part_dict']), self.ac,
                                     cer_score=None,
                                     struct=ast.literal_eval(trna_info['struct']),
                                     div=ast.literal_eval(trna_info['div']),
                                     freq=ast.literal_eval(trna_info['freq']),
                                     seq=ast.literal_eval(trna_info['seq']))
                          for name, trna_info in trna_dict.items()}
        print(f'{len(self.trnas)} Retrieved!...Time elapsed: {time.time() - now}')

    ###############################

    def designs_2_fa(self, output_file_stem, num_seqs='all', ac='CTA'):

        """Exports current tRNAs to fasta file, which can be input to RNAfold.
        output_file_stem has _anticodon_i.fa added.
        num_seqs denotes how many sequences per fa file - sometimes want it broken up e.g. testing or smaller chunks
        In future, may want to code a way to send outputs through rest of pipeline while RNAfold continues,
        although this would be quite complicated with regards to parallel processing.
        num_seqs = "all" returns one large file
        ac denotes anticodon, important for record keeping"""

        now = time.time()
        print(f'Exporting {len(self.trnas)} tRNAs')

        if num_seqs == 'all':
            num_seqs = len(self.trnas)

        def chunks(data):
            # Creates batches
            it = iter(data)
            for i in range(0, len(data), num_seqs):
                yield {k: data[k] for k in islice(it, num_seqs)}

        trna_batches = chunks(self.trnas)

        for i, batch in enumerate(trna_batches, 1):
            with open(f'{output_file_stem}_{ac}_{i}.fa', 'w+') as f:
                for name, trna in batch.items():
                    f.write(f'>{name}\n{trna.seq[ac]}\n')
        print(f'tRNAs Exported!...Time elapsed: {time.time() - now}')

    ###############################

    def fold_filter(self, ac, fold_file, output_dir, synth_name, pattern, iteration=1, freq_thresh=0.3, div_thresh=10,
                    inplace=True, log_file=None):

        """Filters RNAfold output.
        ac is important to save information to correct dictionary key.
        tRNAs are kept with freq and div above freq_thresh and below div_thresh, respectively.
        inplace=True replaces self.trnas with the filtered trnas, otherwise they are just saved to
        self.well_folded
        Method also removes tRNAs without a T(U) in the first unpaired region, and without an A as the first
        base in the D loop"""

        now = time.time()
        print(f'Filtering folds from {fold_file} with frequency >= {freq_thresh} and diversity <= {div_thresh}')

        with open(fold_file, 'r+') as g:
            i = 0
            for line in g:
                i += 1
                if (i + 5) % 6 == 0:
                    seq_name = line[1:].strip()
                elif (i + 4) % 6 == 0:
                    seq = line.strip()
                elif (i + 3) % 6 == 0:
                    struct = line[:len(seq)]
                elif i % 6 == 0:
                    first, second = line.strip().split(';')
                    freq = re.findall("\d+\.\d+", first)[0]
                    div = re.findall("\d+\.\d+", second)[0]
                    self.trnas[seq_name].struct[ac] = struct
                    self.trnas[seq_name].div[ac] = float(div)
                    self.trnas[seq_name].freq[ac] = float(freq)

        divs = [t.div[ac] for t in self.trnas.values()]
        freqs = [t.freq[ac] for t in self.trnas.values()]
        fig, (ax1, ax2) = plt.subplots(1, 2)
        fig.suptitle(f"Structural output {ac}")
        ax1.hist(divs, bins=np.arange(0, 40, 0.5))
        ax1.set_title('Diversity')
        ax2.hist(freqs, bins=np.arange(0, 1, 0.01))
        ax2.set_title('Frequency')
        fig.savefig(output_dir + f'/structure_{synth_name}_{iteration}_{ac}.pdf')
        plt.clf()

        filt_data = {seq_name: trna for seq_name, trna in self.trnas.items()
                     if pattern.match(trna.struct[ac])
                     and trna.freq[ac] >= freq_thresh
                     and trna.div[ac] <= div_thresh}

        filtered = {}
        for name, trna in filt_data.items():
            first_unpaired = trna.struct[ac].find('.')

            indices = re.compile('\(\\.').finditer(trna.struct[ac])
            next(indices)
            d_loop_start = next(indices).span()[0] + 1

            if 'T' in trna.seq[ac][first_unpaired:first_unpaired + 2] and trna.seq[ac][d_loop_start] == 'A':
                filtered.update({name: trna})

        self.well_folded = {seq_name: trna for seq_name, trna in filtered.items()}
        if log_file:
            with open(log_file, 'a') as f:
                f.write(f'Chimeras: {len(self.well_folded)} Anticodon: {ac} '
                        f'Frequency: >={freq_thresh} Diversity: <={div_thresh}\n'
                        f'Pattern: {pattern}\n')
        print(f'There are {len(self.well_folded)} well-folding tRNA designs!...Time elapsed: {time.time() - now}')

        if inplace:
            self.trnas = self.well_folded

    ##################################

    def change_ac(self, new_anticodons):

        """Method copies the part_dict of each tRNA, alters the anticodon to the new anticodon(s),
        then joins the part_dict to create a new sequence. Method calls designs_2_fa to export new sequences
        to fasta file for RNAfold input."""

        now = time.time()

        for new_ac in new_anticodons:
            print(f'Changing anticodon to {new_ac}')
            for name, trna in self.trnas.items():
                trna.part_dict_ = trna.part_dict.copy()
                if isinstance(trna.part_dict_['tRNA32-38*'], Part2):
                    trna.part_dict_['tRNA32-38*'] = trna.part_dict_['tRNA32-38*'].seq.replace(self.ac, new_ac)
                else:
                    trna.part_dict_['tRNA32-38*'] = trna.part_dict_['tRNA32-38*'].replace(self.ac, new_ac)
                trna.seq[new_ac] = ''.join([(trna.part_dict_[part_type] if isinstance(trna.part_dict_[part_type], str)
                                             else trna.part_dict_[part_type].seq)
                                            for part_type in part_order])

        print(f'Anticodon Changed!...Time elapsed: {time.time() - now}')

    #######################################

    def final_filter(self, freq_thresh=0.3, div_thresh=5, percentile_out=20, log_file=None):

        """Final filtering step. First finds average freq and div across all sequences (with different ac)
        and filters based on the given thresholds.
        Method then calculates levenshtein distance to the E. coli tRNA(s) in the same isoacceptor.
        tRNAs in the top quartile (furthest from the E. coli tRNA(s)) are kept."""

        now = time.time()

        print(f'Final filtering step on {len(self.trnas)} tRNAs')

        for name, trna in self.trnas.items():
            trna.avg_freq = sum([float(freq) for freq in trna.freq.values()]) / len(trna.freq)
            trna.avg_div = sum([float(div) for div in trna.div.values()]) / len(trna.div)
        self.trnas = {name: trna for name, trna in self.trnas.items()
                      if trna.avg_freq > freq_thresh
                      and trna.avg_div < div_thresh}

        if len(self.trnas) > 0:
            ecoli_wt = ecoli_df.loc[[self.aa]]
            ecoli_wt['seq'] = ecoli_wt.iloc[:, 1:].apply(lambda x: ''.join(x), axis=1)
            ecoli_wt['seq'] = ecoli_wt.seq.apply(lambda x: x.replace('-', ''))
            wt_seqs = ecoli_wt.seq
            for name, trna in self.trnas.items():
                trna.ec_lev_dist = sum([distance.levenshtein(trna.seq[self.ac], seq)
                                        for seq in wt_seqs]) / len(wt_seqs)
            upper_quartile = np.percentile([trna.ec_lev_dist for trna in self.trnas.values()], percentile_out)
            self.trnas = {name: trna for name, trna in self.trnas.items() if trna.ec_lev_dist >= upper_quartile}

            self.iter_trnas.update(self.trnas)

        if log_file:
            with open(log_file, 'a') as f:
                f.write(f'Final Filter with freq_thresh: {freq_thresh} '
                        f'div_thresh: {div_thresh} percentile_out: {percentile_out}\nChimeras: {len(self.trnas)}\n')
        print(f'Filtering Complete! {len(self.trnas)} remain...Time elapsed: {time.time() - now}')

    ########################################

    def select(self, synth_name, output_dir, automatic=False, log_file=None):

        """Selection step.
        Method calculates levenshtein distance between the tRNA and the native tRNA for this synthetase.
        Idea is that we want to minimise distance to the native sequence, while maximising distance to coli
        tRNAs. Currently, method outputs plot of distance to native against distance to E. coli.
        Ideas for suggesting tRNAs include taking positive face of graph; taking the tRNA closest to native
        then sampling the sequence space away from this point to maintain diversity.
        The method I think I will use is to bin tRNAs with regard to distance from E. coli, then take
        the tRNA closest to native from each bin i.e. taking the postive face of the graph."""

        self.trnas = self.iter_trnas
        synth_ = [synth for synth in self.synths if synth.name == synth_name][0]
        try:
            native_seq = self.huge_df[self.huge_df.seq_id == synth_.trna_id].loc[:, 'tRNA1-7*':'tRNA73-76*'].apply(
                lambda x: ''.join(x), axis=1)
            native_seq = native_seq.iloc[0]
        except TypeError:
            print('Error in database entry for WT sequence.')
            native_seq = input('Please enter sequence of native sequence: ')
        print(native_seq)
        for name, trna in self.trnas.items():
            trna.nat_lev_dist = distance.levenshtein(trna.seq[self.ac], native_seq)

        dict_for_plot = {name: [trna.ec_lev_dist, trna.nat_lev_dist, trna.seq[self.ac]] for name, trna in
                         self.trnas.items()}
        self.df_for_plot = pd.DataFrame.from_dict(dict_for_plot, orient='index')
        self.df_for_plot.columns = ['to_e_coli', 'to_wt', 'seq']


        plt.plot(self.df_for_plot.to_wt, self.df_for_plot.to_e_coli, 'ro')
        plt.xlabel('Levenshtein Distance to WT')
        plt.ylabel('Levenshtein Distance to E. coli')
        plt.savefig(output_dir + f'/Distance_Plot_{synth_name}.pdf')
        plt.clf()

        if not automatic:
            while True:
                user = input("Cut-off point for distance to native ('x' for exit)\n>  ")
                if user.lower() == 'x':
                    print(f'{len(self.trnas)} tRNAs remaining!')
                    break
                else:
                    try:
                        cutoff = float(user)
                        poss = {name: trna for name, trna in self.trnas.items() if trna.nat_lev_dist <= cutoff}
                        print(f'{len(poss)} tRNAs <= {cutoff}.')
                        user_answer = input("Use this cut-off value? ('y' or 'n')\n> ")
                        if user_answer.lower() == 'y':
                            self.trnas = poss
                            print(f'{len(self.trnas)} tRNAs remaining!')
                            break
                        elif user_answer.lower() == 'n':
                            continue
                        else:
                            print('User must enter y or n!')
                            continue
                    except:
                        print('Inappropriate value!')
                        continue
        else:
            start = 5
            end = 30
            for index, i in enumerate(range(start, end)):
                cutoff = i
                poss = {name: trna for name, trna in self.trnas.items() if trna.nat_lev_dist <= cutoff}
                if len(poss) < 200 or index == 0:
                    n = abs(200 - len(poss))
                else:
                    n_ = abs(200 - len(poss))
                    if n_ <= n:
                        self.trnas = poss
                    else:
                        cutoff = i - 1
                        self.trnas = {name: trna for name, trna in self.trnas.items() if trna.nat_lev_dist <= cutoff}
                    break

        if log_file:
            with open(log_file, 'a') as f:
                f.write(f'Selected {len(self.trnas)} Chimeras with maximum distance to WT {cutoff}\n')
            # self.df_for_plot = self.df_for_plot[self.df_for_plot.to_wt < cutoff]
            # plt.plot(self.df_for_plot.to_wt, self.df_for_plot.to_e_coli, 'ro')
            # plt.pause(0.01)

    #############################################

    def cluster_select(self, num_seqs=4, inplace=False, cluster=True, damp=0.8, convergence_iter=15, preference=-11,
                       log_file=None):

        def cluster_meat(damping):

            sample = np.asarray([trna for trna in self.trnas.values()])
            m = np.zeros((len(sample), len(sample)))
            tril_idx_rows, tril_idx_cols = np.tril_indices_from(m, k=-1)
            m[(tril_idx_rows, tril_idx_cols)] = [distance.levenshtein(sample[i].seq[self.ac], sample[j].seq[self.ac])
                                                 for i, j in zip(tril_idx_rows, tril_idx_cols)]
            m = (m + m.T) * -1

            affprop = AffinityPropagation(affinity="precomputed", damping=damping, max_iter=1000,
                                          convergence_iter=convergence_iter,
                                          preference=preference, random_state=None)
            print(f'Clustering with damping: {damping}')
            affprop.fit(m)

            try:
                for cluster_id in np.unique(affprop.labels_):
                    # identify exemplar for each cluster (found at centre)
                    sample[affprop.cluster_centers_indices_[cluster_id]].exemplar = True
                    for element_list in np.nonzero(affprop.labels_ == cluster_id):
                        for element in element_list:
                            # assign cluster_id to each part in all_parts
                            sample[element].cluster_id = cluster_id
                self.exemplar_trnas = {name: trna for name, trna in self.trnas.items() if trna.exemplar}
                if log_file:
                    with open(log_file, 'a') as f:
                        f.write(f'{len(self.exemplar_trnas)} Cluster Exemplars with Damping: {damping}\n')

            except IndexError:
                damping += 0.1
                if damping >= 1.0:
                    print('Clustering Failed...Taking random sample of 40')
                    self.exemplar_trnas = random.sample(self.trnas.items(), 40)
                    self.exemplar_trnas = {trna_tup[0]: trna_tup[1] for trna_tup in self.exemplar_trnas}
                else:
                    cluster_meat(damping)

        now = time.time()
        print('Clustering tRNAs...')

        for trna in self.trnas.values():
            trna.exemplar = False
            trna.cluster_id = None
        if cluster:
            cluster_meat(damping=damp)
        else:
            self.exemplar_trnas = self.trnas

        read_back = [{name: trna} for name, trna in self.exemplar_trnas.items()]
        print(f'{len(self.exemplar_trnas)} exemplars chosen...Maximising diversity...Time elapsed: {time.time() - now}')

        if num_seqs < len(self.exemplar_trnas):
            diverse_trnas, distances = max_dist_parallel_memo([tRNA for tRNA in self.exemplar_trnas.values()], num_seqs,
                                                              self.ac, type='tRNA')
            self.final_trnas = {trna_name: trna for trna_name, trna in self.exemplar_trnas.items()
                                if trna.seq[self.ac] in diverse_trnas}
        else:
            distances = []
            self.final_trnas = self.exemplar_trnas

        if log_file:
            log_string = ''
            for trna_name, trna in self.final_trnas.items():
                log_string += f">{trna_name}\n{trna.seq[self.ac]}\n{trna.struct[self.ac]}\n" \
                              f"Frequency: {trna.freq[self.ac]} Diversity: {trna.div[self.ac]}\n"
            with open(log_file, 'a') as f:
                f.write(f'{num_seqs} Chimeras Chosen\n')
                if distances:
                    f.write(f'Minimum Distance: {distances}\n')
                f.write(log_string)
        print(f'Designs Finished!...Time Elapsed: {time.time() - now}')

        # if inplace:
        #     self.trnas = {name: trna for name, trna in self.exemplar_trnas.items() if name in chosen_names}


###################################################################################################################


class Part2(object):
    """Part2 class contains information about each part. Initialised with the sequence, part type (region),
    amino acid, tRNA_id (will do in future), and the aligned sequence (mainly for D loop)."""

    def __init__(self, seq, region, aa, trna_id, aligned, iso):
        self.seq = seq

        self.region = region
        self.aa = aa
        self.trna_id = trna_id
        self.aligned = aligned
        self.iso = iso
        self.parent = None
        self.sub_parts = None
        self.cer_score = {}

        # base_range is the base positions of that part e.g. tRNA8-9* part has base_range [8, 9]
        # (or more accurately, range(8, 10))
        try:
            self.base_range = part_to_range_2[self.region]
        except KeyError:
            self.base_range = part_to_range[self.region]

        # merged parts don't have a simple range, but a tuple of ranges in the part_to_range_2 dictionary
        # base_range instead forms a list with e.g. for tRNA1-7_66-72* [1, 2,...,7, '_', 66, 67,..., 72]
        if isinstance(self.base_range, tuple):
            self.base_range_ = self.base_range
            self.base_range = []
            first = True
            for r in self.base_range_:
                if not first:
                    self.base_range.append('_')
                for i in r:
                    self.base_range.append(i)
                first = False

        # Now create seq_dict e.g. for tRNA8-9*, seq_dict = {8: 'T', 9: 'A'}
        if self.region == 'tRNA26_44-48*':
            # Different alignment for variable loop - now I think about this, I could perform this
            # in previous step, and assign the aligned version to the part_tuple.align attribute, simplifying this
            # However this works as is, so don't touch
            # Added 2 to each index compared to old class due to addition of '26base_'
            self.seq_dict = {44: self.seq[2], 48: self.seq[-1], 45: self.seq[3], 26: self.seq[0]}
            if len(self.seq) > 5:
                self.seq_dict[47] = self.seq[-2]
            else:
                self.seq_dict[47] = '-'
            if len(self.seq) > 6:
                self.seq_dict[46] = self.seq[4]
            else:
                self.seq_dict[46] = '-'

        elif self.region == 'tRNA44-48*':
            self.seq_dict = {44: self.seq[0], 48: self.seq[-1], 45: self.seq[1]}
            if len(self.seq) > 3:
                self.seq_dict[47] = self.seq[-2]
            else:
                self.seq_dict[47] = '-'
            if len(self.seq) > 4:
                self.seq_dict[46] = self.seq[2]
            else:
                self.seq_dict[46] = '-'

        else:
            try:
                self.seq_dict = {i: base for i, base in zip(self.base_range, self.aligned)}
            except:
                print(
                    f'seq: {self.seq} align: {self.aligned} region: {self.region} base_range: {self.base_range} id: {self.trna_id}')
                raise Exception

        # id_dict is not associated with the class, but defined outside
        self.id_d = {aa: {pos: base for pos, base in d.items() if pos in self.base_range} for aa, d in id_dict.items()}

        # id_bases is dict of canonical positions within the part with the bases found at them in E. coli tRNAs
        # Feel like I could do this in a one line list comprehension to be more elegant
        self.id_bases = {}
        for pos in self.base_range:
            # Create empty list for each position in base_range within id_bases
            self.id_bases[pos] = []
            # iterate through id_d ('Arg': {2: {'A'}, 10: {'G', 'A'},...},...)
            for aa, d in self.id_d.items():
                if pos in d.keys():
                    # append list of bases found at that position
                    self.id_bases[pos].append(list(d[pos]))
        # For each position in id_bases, use Counter class to create dict-like structure counting the frequency
        # of each base at that position across all E. coli tRNAs
        # This is the way my ID score approach differs from Daniele's
        self.id_bases = {pos: Counter([b for base in base_list for b in base])
                         for pos, base_list in self.id_bases.items()
                         if base_list}

        pos_scores = []
        # Low score more orthogonal
        for pos, counts in self.id_bases.items():
            try:
                # Score for each position is the fraction the base in the part comes up in E. coli tRNAs at that
                # position
                # e.g. our part is a tRNA8-9* region (doesn't matter which isoacceptor)
                # Let's say there are 7 E. coli tRNAs with an ID element at position 8.
                # These bases are ['A', 'A', 'T', 'G', 'C', 'A', 'C']
                # Our part has the sequence 'TA', which will have the seq_dict {8: 'T', 9: 'A'}
                # In this example, no E. coli synthetase has an ID element at position 9
                # Therefore, the id_score for our part is 1/7, since 1 of 7 E. coli tRNAs match the ID element
                # If we could weight the importance of ID elements, this would likely improve our designs, but
                # it is difficult to see how this would be quantified accurately
                pos_score = counts[self.seq_dict[pos]] / sum(counts.values())
            except KeyError:
                # Kept this in here since I was having some bugs, but found the problems
                print(f'{self.region}, {self.seq_dict}, {self.aligned}')
            pos_scores.append(pos_score)
        try:
            self.id_score = sum(pos_scores) / len(pos_scores)
        except ZeroDivisionError:
            #             self.id_score = 'a'
            self.id_score = 0

        if '_' in self.region:
            part_type_1, part_type_2 = self.region.split('_')
            part_type_1 += '*'
            part_type_2 = 'tRNA' + part_type_2

            part_1_seq, part_2_seq = self.seq.split('_')
            if part_type_1 != 'tRNA14-21*':
                part_1_align = part_1_seq
            else:
                part_1_align = self.aligned.split('_')[0]

                # Some manually input parts don't align, so need to do slower dataframe lookup
                if part_1_align is np.nan:
                    part_1_align = self.iso.huge_df[self.iso.huge_df['tRNA14-21*'] == part_1_seq].iloc[0][
                        'tRNA14-21* aligned']

            self.sub_parts = {part_type_1: SubPart(part_1_seq, part_type_1, self.aa,
                                                   self.trna_id, part_1_align, self, self.iso),
                              part_type_2: SubPart(part_2_seq, part_type_2, self.aa,
                                                   self.trna_id, part_2_seq, self, self.iso)}

        self.cluster_id = None
        self.exemplar = False

        for aa_ in id_dict.keys():
            self.cer_score[aa_] = self.cer_scorer(aa_)

    def cer_scorer(self, aa_):
        id_ = 0
        num = 0
        subdict = self.iso.non_id_id_dict[aa_]
        for pos, base in self.seq_dict.items():
            if pos in subdict.keys():
                num += 1
                if base in subdict[pos]:
                    id_ += 1
                else:
                    id_ -= 1
        if num == 0:
            num = 1
        return {'Score': id_, 'Positions': num}

    def __repr__(self):
        """String representation of Part2 object"""
        return f'Seq:{self.seq} ID:{self.trna_id} Score: {self.id_score} Cluster: {self.cluster_id} Exemplar: {self.exemplar}'

    def __lt__(self, other):
        """less than dunder method required for sorting all_parts by score"""
        return self.id_score < other.id_score


###################################################################################################################


class SubPart(Part2):

    def __init__(self, seq, region, aa, trna_id, aligned, parent, iso):
        super().__init__(seq, region, aa, trna_id, aligned, iso)
        self.id_score = None
        self.cluster_id = None
        self.exemplar = None
        self.parent = parent


###################################################################################################################


class tRNA(object):
    """These objects will be tRNA designs.
    Class will have methods to handle and use RNAFold data"""

    def __init__(self, part_dict, ac, cer_score=None, struct=None, div=None, freq=None, seq=None):
        if freq is None:
            freq = {}
        if div is None:
            div = {}
        if struct is None:
            struct = {}
        self.part_dict = part_dict
        self.ac = ac
        self.seq = seq
        self.struct = struct
        self.div = div
        self.freq = freq
        if not seq:
            self.seq = {}
            self.seq[self.ac] = ''.join([self.part_dict[part_type].seq for part_type in part_order])
        else:
            self.seq = seq
        self.cer_score = cer_score

    def trna_cer_scorer(self):
        parent_set = {part.parent if part.parent else part for part in self.part_dict.values()}
        self.cer_score = {aa: sum([part.cer_score[aa]['Score'] for part in parent_set]) / sum(
            [part.cer_score[aa]['Positions'] for part in parent_set])
                          for aa in id_dict.keys()}

    def __lt__(self, other):
        """less than dunder method required for sorting all_parts by score"""
        return self.cer_score < other.cer_score

    def __repr__(self):

        return f'Seq: {self.seq} Cervettini Score: {self.cer_score}'

#############################################################################################################
