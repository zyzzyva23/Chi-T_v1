import argparse

base_comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def oligo_maker(files, output_name, fwd_5, fwd_3, rev_5, rev_3, ac):

    out = pd.DataFrame()
    for file in files:
        df = pd.read_csv(file)
        df.seq = df.seq.apply(ast.literal_eval)
        df['Sequence 1'] = df.seq.apply(lambda x: fwd_5 + x[ac] + fwd_3)
        df['Sequence 2'] = df.seq.apply(
            lambda x: rev_5 + ''.join([base_comp[base] for base in x[ac][::-1]]) + rev_3)
        df = df.rename(columns={'name': 'Name'})
        df = df.drop(columns=['seq', 'part_dict', 'cer_score', 'struct', 'div', 'freq'])
        out = pd.concat([out, df])

    out.to_csv(output_name + '.csv', index=False)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("design_files", help='File names of final tRNA designs', nargs='+')
    parser.add_argument("output_file", help="Name of output file")
    parser.add_argument("-f5", '--forward_5', help="5' attachment to forward primer", default='GGCCGC')
    parser.add_argument("-f3", '--forward_3', help="3' attachment to forward primer", default='CTGCA')
    parser.add_argument("-r5", '--reverse_5', help="5' attachment to reverse primer", default='GATCTGCAG')
    parser.add_argument("-r3", '--reverse_3', help="3' attachment to reverse primer", default='GC')
    parser.add_argument("-a", '--anticodon', help="Designs to be written to output file (will have this anticodon)",
                        default='CTA')
    args = parser.parse_args()

    oligo_maker(args.design_files, args.output_file, args.forward_5, args.forward_3,
                args.reverse_5, args.reverse_3, args.anticodon)

