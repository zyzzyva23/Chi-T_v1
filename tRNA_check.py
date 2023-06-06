import pandas as pd
import argparse
from pprint import pprint


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help='Clean tRNADB-CE file')
    parser.add_argument("tRNA_id", help='ID of tRNA to check')
    args = parser.parse_args()

    df = pd.read_csv(args.file, usecols=['seq_id', 'gen_id', 'Phylum/Class', 'Species', 'Amino Acid',
                                         'Anticodon', 'tRNA1-7*', 'tRNA8-9*', 'tRNA10-13*', 'tRNA14-21*',
                                         'tRNA22-25*', 'tRNA26*', 'tRNA27-31*', 'tRNA32-38*', 'tRNA39-43*',
                                         'tRNA44-48*', 'tRNA49-53*', 'tRNA54-60*', 'tRNA61-65*', 'tRNA66-72*',
                                         'tRNA73-76*'])
    if args.tRNA_id in df.seq_id.values:
        print('True')
        pprint(df[df.seq_id == args.tRNA_id])
    else:
        print('False')