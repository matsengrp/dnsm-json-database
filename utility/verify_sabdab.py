import pandas as pd
import os
import sys
from debug_util import *
import sabdab_pipeline
import dnsm_pipeline

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sabdab_true = f"{SCRIPT_DIR}/test_examples/sabdab_summary_2024-01-26_abid_info.tsv"
sabdab_test = f"{SCRIPT_DIR}/_ignore/sabdab_2024-01-26/sabdab_summary_all_2024-01-26_abid_info.human.2000-3000.tsv"


def compare_abid_info_files(true_path, test_path):
    true_df = pd.read_table(true_path)
    test_df = pd.read_table(test_path)
    abids = set(true_df.abid) & set(test_df.abid)

    for abid in abids:
        true_row = true_df[true_df.abid == abid]
        test_row = test_df[test_df.abid == abid]

        assert len(true_row) == 1
        assert len(test_row) == 1

        true_row = true_row.iloc[0]
        test_row = test_row.iloc[0]

        row_cmp = True
        for col in true_df.columns:
            if col == 'pdbfile':
                continue

            true_data = true_row[col]
            test_data = test_row[col]

            col_cmp = (str(true_data) != str(test_data))
            if col_cmp:
                row_cmp = False
                print(
                    f'COL abid: {abid}, col: {col}\n\ttrue: {true_data}, test: {test_data}')

        print(f'ROW abid: {abid}, row_compare: {row_cmp}')


def compare_chainseqs(sabdab_path=sabdab_true):
    df = pd.read_table(sabdab_path)

    for index, row in df.iterrows():
        cprint(f'abid: {row.abid}', color=colors.GREEN)

        dict_data, raw_data, type_counts = dnsm_pipeline.fetch_and_parse_fasta(
            row.pdbid)
        if dict_data == None:
            cprint('error: failed to parse fasta.', color=colors.RED)
            continue

        seq_fa = {}
        for key, data in dict_data.items():
            seq_fa[''.join(data['chainid'])] = data['aa_seq']

        HL = row.abid[4:6]
        seq_sabdab = {
            'a': [HL[1], row.chainseq_a],
            'b': [HL[0], row.chainseq_b]
        }
        for key1, (chainid1, chainseq1) in seq_sabdab.items():
            print(
                f'{key1}_{chainid1}: {len(chainseq1)} {chainseq1[0:5]} {seq_fa.keys()}\n\t', end='')
            for chainid2, chainseq2 in seq_fa.items():
                if chainid2.find(chainid1) != -1:
                    cmp = (chainseq1 == chainseq2)
                    print(
                        f'{chainid2} {cmp} {len(chainseq2)} {chainseq2[0:5]}, ', end='')
            print('\n', end='')
