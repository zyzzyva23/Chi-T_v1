import argparse
from ortho_classes import Isoacceptor2, id_dict, synth_clean
from parallel import rnafold_in_parallel
import pandas as pd
from itertools import zip_longest
from pathlib import Path
import time
import re
import sys


if __name__ == '__main__':


    print(" ::::::::  :::    ::: :::::::::::               :::::::::::\n"
          ":+:    :+: :+:    :+:     :+:                       :+:\n"
          "+:+        +:+    +:+     +:+                       +:+\n"
          "+#+        +#++:++#++     +#+     +#++:++#++:++     +#+\n"
          "+#+        +#+    +#+     +#+                       +#+\n"
          "#+#    #+# #+#    #+#     #+#                       #+#\n"
          " ########  ###    ### ###########                   ###     ")

    print("""
                                             *  ((((
                                            (((* ** (((.
                                                    .(,
                                                    %%%#
                                            (%%&/%%%
                                            *%#  .(
                                            %%%.%%%#
                                           %%%% %%,
                                            #( *%%(
                                          /%%&#%%%
                                          /%#  .*
                                          %%%.%%%(
                                         %%%#.%%*
       ...  ...                           (* *%%(
    ......  .. ....                (%%%/#%%&#%%%                   ... ....
    ...         .,             (%%%           ,*##/ ##/ ##* ,#*  /.,..    ....
  ....           ((((((((((((,((*              ,%%*.#%#.#&# #&# #%#/         *
   ..            (%* ,&   #* .#(,              ########,#%# #&# #&#         ....
   ...          ,((/*(((/(((/(((,              ###, ,*  *(. *#/ /##        ...
     .... ...,...              ##             ###*                ..../... ...
          ...                 (##(**       *####*
                                 ###%&#######(
                                (##( ,#
                                 ## ####
                               ###%(##(
                               ##   ##
                              (##(%%##
                             #### ##
                         **** ** (##(
                        **,         ****
                        ***          /
                         ***       ,***
                          , ****,***                                            
                             .,                                            """)
    print('Welcome to Chi-T!')
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help='Clean tRNADB-CE file')
    parser.add_argument("synth_file", help='File containing synthetase information')
    parser.add_argument("synth_name", help='Name of synthetase matching to entry in synth_file',
                        nargs='+')
    parser.add_argument("amino_acid", help='Amino Acid specified for tRNA generation')
    parser.add_argument('-o', '--output_directory', help='Directory to store output files', default='')
    parser.add_argument("-cp", "--cluster_parts", type=int,
                        help='Number of parts for each part type to cluster', default=200)
    parser.add_argument("-l", "--length_filt", help='Filter chimeras if 79 nts or longer (or specified value)',
                        nargs='?', const=79, type=int)
    parser.add_argument("-cf", "--cervettini_filt", nargs='+', help='Parameters for Cervettini Filtering in order '
                                                                    'start_stringency, minimum_stringency, target '
                                                                    'number of chimeras, and step size',
                        type=float, default=[0.5, 0.2, 2500000, 0.05])
    parser.add_argument('-a', '--anticodons', nargs='+', help='Anticodons to iterate through')
    parser.add_argument('-F', '--final_frequency', help='Average frequency across anticodons', default=0.35, type=float)
    parser.add_argument('-D', '--final_diversity', help='Average diversity across anticodons', default=9.0, type=float)
    parser.add_argument('-f', '--frequency', help='Frequency threshold for each anticodon', default=0.3, type=float)
    parser.add_argument('-d', '--diversity', help='Diversity threshold for each anticodon', default=10.0, type=float)
    parser.add_argument('-m', '--automatic', help='No user input required', action='store_true')
    parser.add_argument('-i', '--initial', help='If true, save initial chimeras to csv (only if needed for analysis)',
                        action='store_true')
    parser.add_argument('-p', '--pattern', help='Specify a csv file with synth name and regex string column')
    parser.add_argument('-t', '--num_tRNAs', help='Number of designs to output', default=4, type=int)
    parser.add_argument('-ip', '--id_part_change', help='Identity parts that should be chimerified (except ID element)',
                        nargs='+')
    args = parser.parse_args()

    if args.cervettini_filt and len(args.cervettini_filt) > 4:
        raise Exception("Too many filtering parameters!")
    else:
        args.cervettini_filt = [(i if i is not None else j)
                                for i, j in zip_longest(args.cervettini_filt, [0.5, 0, 2500000, 0.05])]
    cf_start, cf_min, cf_targ, cf_ss = args.cervettini_filt

    first_ac = args.anticodons[0]

    Path(f'{args.output_directory}/folding').mkdir(parents=True, exist_ok=True)
    Path(f'{args.output_directory}/plots').mkdir(parents=True, exist_ok=True)

    log_file = f'{args.output_directory}/log_file.txt'
    with open(log_file, 'w') as f:
        f.write('Chi-T\n' +
                str(time.time()) + '\n' +
                str(args) + '\n\n')

    df = pd.read_csv(args.file, usecols=['seq_id', 'Amino Acid', 'tRNA1-7*', 'tRNA8-9*', 'tRNA10-13*', 'tRNA14-21*',
                                         'tRNA22-25*', 'tRNA26*', 'tRNA27-31*', 'tRNA32-38*', 'tRNA39-43*',
                                         'tRNA44-48*', 'tRNA49-53*', 'tRNA54-60*', 'tRNA61-65*', 'tRNA66-72*',
                                         'tRNA73-76*', 'tRNA14-21* aligned', 'tRNA1-7_66-72*', 'tRNA10-13_22-25*',
                                         'tRNA14-21_54-60*', 'tRNA14-21_54-60* aligned', 'tRNA26_44-48*',
                                         'tRNA27-31_39-43*', 'tRNA49-53_61-65*'],
                     dtype={'Amino Acid': 'category', 'tRNA8-9*': 'category', 'tRNA26*': 'category',
                            'tRNA73-76*': 'category'},
                     engine='c')

    trna_pattern = re.compile(
        '^\\({5,8}\\.{1,3}\\({4}\\.{5,}\\){4}\\.*\\({4,9}\\.{7}\\){4,9}.*\\.\\({5}\\.{2,}\\){5}\\){5,8}\\.{3,}$')
    ile_pattern = re.compile(
        '^\\({5,8}\\.{1,3}\\({3}\\.{5,}\\){3}\\.*\\({4,9}\\.{7}\\){4,9}.*\\.\\({5}\\.{2,}\\){5}\\){5,8}\\.{3,}$')
    pat_dict = {synth_name: (trna_pattern if args.amino_acid != 'Ile' else ile_pattern)
                for synth_name in args.synth_name}
    if args.pattern:
        pat_df = pd.read_csv(args.pattern, header=None)
        pat_dict_ = dict(zip(pat_df.iloc[:, 0], pat_df.iloc[:, 1]))
        for k, v in pat_dict_.items():
            pat_dict_[k] = v.replace('\\\\', '\\')
        pat_dict.update(pat_dict_)

    df = df[df['Amino Acid'] == args.amino_acid]
    df['tRNA32-38*'] = df['tRNA32-38*'].apply(lambda x: x[:2] + first_ac + x[-2:] if isinstance(x, str) else x)
    synth_df = synth_clean(args.synth_file)
    synth_df = synth_df[synth_df.synth.isin(args.synth_name)]
    # synth_df_ = synth_df[synth_df.synth in args.synth_name]
    # if not all([seq_id in df.seq_id for seq_id in synth_df_.trna_id]):
    #     raise Exception("One or more tRNA ID in synth file not found in database. Check file.")
    # total_iter = len(args.synth_name)*args.num_iterations

    iso = Isoacceptor2(synth_df, id_dict, args.amino_acid, df, ac=first_ac, id_part_change=args.id_part_change)

    synth_name = args.synth_name[0]
    pattern = re.compile(pat_dict[synth_name])
    iso.iter_trnas = {}
    with open(log_file, 'a') as f:
        f.write(f'##################################################\n\nChi-T Run for {synth_name}\n')

        iso.cluster_parts(args.cluster_parts, synth_name=synth_name, clust_id_parts=False, log_file=log_file,
                          automatic=args.automatic)
        iso.chimera(synth_name, length_filt=args.length_filt, log_file=log_file)
        iso.cervettini_filter(args.output_directory, synth_name, start_stringency=cf_start,
                              min_stringency=cf_min, target=cf_targ, step_size=cf_ss, log_file=log_file)
        if args.initial:
            iso.store_trnas(f'{args.output_directory}/{synth_name}_initial.csv')

        print('Folding...')
        rnafold_in_parallel(iso, f'{args.output_directory}/folding/{synth_name}_para', first_ac)
        iso.fold_filter(first_ac,
                        f'{args.output_directory}/folding/{synth_name}_para_{first_ac}_complete_fold.out',
                        args.output_directory, synth_name, pattern, log_file=log_file, freq_thresh=args.frequency,
                        div_thresh=args.diversity)

        for ac in args.anticodons[1:]:
            if not iso.trnas:
                break
            iso.change_ac([ac])
            rnafold_in_parallel(iso, f'{args.output_directory}/folding/{synth_name}_para', ac)

            iso.fold_filter(ac, f'{args.output_directory}/folding/{synth_name}_para_{ac}_complete_fold.out',
                            args.output_directory, synth_name, pattern, log_file=log_file, freq_thresh=args.frequency,
                            div_thresh=args.diversity)

        if iso.trnas:
            iso.final_filter(freq_thresh=args.final_frequency, div_thresh=args.final_diversity, percentile_out=0,
                             log_file=log_file)
            iso.store_trnas(f'{args.output_directory}/{synth_name}_finalfold.csv', compress=True)

        iso.select(synth_name, args.output_directory, log_file=log_file, automatic=args.automatic)
        iso.store_trnas(f'{args.output_directory}/{synth_name}_selected.csv')

        cluster = len(iso.trnas) > 40
        iso.cluster_select(cluster=cluster, num_seqs=args.num_tRNAs, log_file=log_file)

        if not args.automatic:
            while True:
                accept = input('Use these tRNAs? (y)es or (n)o\n> ')
                if accept.lower() == 'y':
                    iso.trnas = iso.final_trnas
                    iso.store_trnas(f'{args.output_directory}/{synth_name}_final_four.csv')
                    iso.iter_trnas = {}
                    break
                elif accept.lower() == 'n':
                    print("For help, please refer to our manuscript, or the GitHub page https://github.com/zyzzyva23/Chi-T_v1")
                    print("Thank you for trying Chi-T!")
                    sys.exit()
                else:
                    print('Inappropriate value!')
                    continue
        else:
            iso.trnas = iso.final_trnas
            iso.store_trnas(f'{args.output_directory}/{synth_name}_final_designs.csv')