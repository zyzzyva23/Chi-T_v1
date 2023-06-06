import sys

import pandas as pd
from pprint import pprint
import argparse
from cleanup import alignment_dict, d_loop_align, d_loop_extend


def merger(df_):
    # Remove tRNAs whose D-loops don't align
    # May try and align these later
    df = df_
    df.loc[:, 'tRNA1-7_66-72*'] = df.apply(lambda x: x['tRNA1-7*'] + '_' + x['tRNA66-72*'], axis=1)
    df.loc[:, 'tRNA10-13_22-25*'] = df.apply(lambda x: x['tRNA10-13*'] + '_' + x['tRNA22-25*'], axis=1)
    df.loc[:, 'tRNA14-21_54-60*'] = df.apply(lambda x: x['tRNA14-21*'] + '_' + x['tRNA54-60*'], axis=1)
    df.loc[:, 'tRNA14-21_54-60* aligned'] = df.apply(lambda x: x['tRNA14-21* aligned'] + '_' + x['tRNA54-60*'], axis=1)
    df.loc[:, 'tRNA26_44-48*'] = df.apply(lambda x: x['tRNA26*'] + '_' + x['tRNA44-48*'], axis=1)
    df.loc[:, 'tRNA27-31_39-43*'] = df.apply(lambda x: x['tRNA27-31*'] + '_' + x['tRNA39-43*'], axis=1)
    df.loc[:, 'tRNA49-53_61-65*'] = df.apply(lambda x: x['tRNA49-53*'] + '_' + x['tRNA61-65*'], axis=1)
    return df


def aligner(df_, align_dict):
    df_['tRNA14-21* aligned'] = df_['tRNA14-21*'].map(d_al_dict).fillna(df_['tRNA14-21*'])
    df_['tRNA14-21* aligned'] = df_['tRNA14-21* aligned'].apply(d_loop_extend)
    df_['tRNA14-21* aligned'] = df_['tRNA14-21* aligned'].apply(d_loop_align, args=(align_dict,))
    return df_



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("file", help='Clean tRNADB-CE file')
    parser.add_argument('alignment_file', help='File for aligning D-loop')
    parser.add_argument("-n", "--new_name", help='Optional name of new file')
    args = parser.parse_args()

    print(f'Add tRNAs to {args.file}')

    df = pd.read_csv(args.file)
    d_al_dict = alignment_dict(args.alignment_file)

    while True:
        seq_id = input('Name of tRNA: ')
        gen_id = input('Genome ID (optional): ')
        pclass = input('Phylum/Class (optional): ')
        species = input('Species (optional): ')
        aa = input('Amino Acid: ')
        ac = input('Anticodon (optional): ')
        t1_7 = input("Acceptor stem sequence 5' (1-7): ")
        t8_9 = input('Bases 8-9 sequence: ')
        t10_13 = input("D-arm 5' sequence (10-13): ")
        t14_21 = input('D-loop sequence (14-21): ')
        t22_25 = input("D-arm 3' sequence (22-25): ")
        t26 = input('Base 26: ')
        t27_31 = input("Anticodon-stem 5' sequence (27-31): ")
        t32_38 = input('Anticodon-stem loop sequence (32-38): ')
        t39_43 = input("Anticodon-stem 3' sequence (39-43): ")
        t44_48 = input('Variable loop sequence (44-48): ')
        t49_53 = input("T-arm 5' sequence (49-53): ")
        t54_60 = input('T-loop (54-60): ')
        t61_65 = input("T-arm 3' sequence (61-65): ")
        t66_72 = input("Acceptor stem 3' sequence (66-72): ")
        t73_76 = input('Discriminator base + CCA (73-76): ')

        new_line = pd.DataFrame({'seq_id': [seq_id],
                                 'gen_id': [gen_id],
                                 'Phylum/Class': [pclass],
                                 'Species': [species],
                                 'Amino Acid': [aa],
                                 'Anticodon': [ac],
                                 'Upstream seq.': [None],
                                 'tRNA1-7*': [t1_7],
                                 'tRNA8-9*': [t8_9],
                                 'tRNA10-13*': [t10_13],
                                 'tRNA14-21*': [t14_21],
                                 'tRNA22-25*': [t22_25],
                                 'tRNA26*': [t26],
                                 'tRNA27-31*': [t27_31],
                                 'tRNA32-38*': [t32_38],
                                 'tRNA39-43*': [t39_43],
                                 'tRNA44-48*': [t44_48],
                                 'tRNA49-53*': [t49_53],
                                 'tRNA54-60*': [t54_60],
                                 'tRNA61-65*': [t61_65],
                                 'tRNA66-72*': [t66_72],
                                 'tRNA73-76*': [t73_76]})

        new_line_align = aligner(new_line, d_al_dict)
        if pd.isnull(new_line_align['tRNA14-21* aligned'].iloc[0]):
            print('New D-loop did not align!')
            print('Must align manually.')
            print('Use - for empty base.')
            print('Sequence must be 8 nts long.')
            print('GG should be at positions 18 and 19 (5 and 6 in sequence).')
            aligned = input('Aligned D-loop sequence: ')
            new_line_align['tRNA14-21* aligned'] = aligned

        new_line_merged = merger(new_line_align)

        while True:
            pprint(new_line)
            check = input('tRNA reviewed? (y/n): ')
            if check.lower() == 'y':
                df = pd.concat([df, new_line_merged])
                while True:
                    another = input('Add another tRNA? (y/n): ')
                    if another.lower() == 'y':
                        break
                    elif another.lower() == 'n':
                        if args.new_name:
                            df.to_csv(args.new_name, index=False)
                        else:
                            df.to_csv(args.file, index=False)
                        print('Program complete.')
                        sys.exit()
                    else:
                        print('Innapropriate value!')
                        print('Must enter y or n')
                        continue
            elif check.lower() == 'n':
                print('Restarting entry.')
                break
            else:
                print('Innapropriate value!')
                print('Must enter y or n')
                continue
            break




