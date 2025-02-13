import os
import sys
import pandas as pd
import argparse
import pprint as pp
import json
import glob
from sabdab_utility import *


def list_to_string(my_list):
    my_str = ','.join(my_list)
    return my_str


def dict_to_string(my_dict):
    my_str = []
    for key, val in my_dict.items():
        my_str.append(f'{key}={val}')
    return list_to_string(my_str)


def string_to_list(my_str):
    my_list = [x.split('=') for x in my_str.split(',')]
    return my_list


def string_to_dict(my_str):
    my_list = string_to_list(my_str)
    my_dict = {}
    for key, val in my_list:
        my_dict[key] = val
    return my_dict


# default args
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_PATH = f'{SCRIPT_DIR}/../_output/sabdab_summary_for_dnsm.tsv'
OUTPUT_PATH = f'{SCRIPT_DIR}/../_output/sabdab_summary_for_dnsm.json'
OG_SABDAB_1_PATH = f'{SCRIPT_DIR}/../_output/sabdab_temp/sabdab/sabdab_summary.1.tsv'
OG_SABDAB_2_PATH = f'{SCRIPT_DIR}/../_output/sabdab_temp/sabdab/sabdab_summary.2.tsv'
OG_SABDAB_PATH = OG_SABDAB_2_PATH
LOG_PATH = f'{SCRIPT_DIR}/../_output/dnsm_temp/log.dnsm_pipeline.txt'
JSON_DIR = f'{SCRIPT_DIR}/../_output/dnsm_output'

COLUMNS = ['pdbid', 'abid', 'organism', 'ja', 'jb', 'va', 'vb']
RENAME = {'ja': 'jl', 'jb': 'jh', 'va': 'vl', 'vb': 'vh'}
# orientation: (split | records | index | columns | values)
ORIENT = 'records'

COLUMNS = list_to_string(COLUMNS)
RENAME = dict_to_string(RENAME)

# file data
LOG_HEADER = ['id', 'jobid', 'pdbid', 'abid', 'status', 'match', 'miss']
# ACCEPTABLE_STATUSES = ['success']
ACCEPTABLE_STATUSES = ['success', 'error:seq_structure_mismatch']


def parse_args():
    parser = argparse.ArgumentParser(
        description='Converts sabdab summary tsv file to a json file.')
    parser.add_argument('--input-path', '-i',
                        default=INPUT_PATH, help='input sabdab summary file')
    parser.add_argument('--output-path', '-o',
                        default=OUTPUT_PATH, help='output json file')
    parser.add_argument('--log-path', default=LOG_PATH,
                        help='path to dnsm logfile')
    parser.add_argument('--og-sabdab-path', default=OG_SABDAB_PATH,
                        help='path to original sabdab summary file')
    parser.add_argument('--json-dir', default=JSON_DIR,
                        help='path to dnsm output directory')
    parser.add_argument(
        '--column', '-c', default=COLUMNS, help='list of fields to output. default: ALL')
    parser.add_argument(
        '--filter', '-f', help='list of filter queries. e.g, organism=human_ig')
    parser.add_argument(
        '--rename', '-r', default=RENAME, help='list of columns to rename. e.g, va=vh')
    parser.add_argument(
        '--orient', default=ORIENT, help='orientation of json output. [split, records, index, columns, values]')
    args = parser.parse_args()

    if args.column:
        args.columns = args.column.split(',')
    if args.filter:
        args.filters = string_to_list(args.filters)
    if args.rename:
        args.rename = string_to_dict(args.rename)
    return args


def main(args):
    gene_dict = {
        'va': [],
        'vb': [],
        'ja': [],
        'jb': []
    }

    df = pd.read_csv(args.input_path, sep='\t')
    df_init = df.copy()
    print(f'FINAL:: len: {len(df)}')
    print(f'FINAL:: pdbids: {len(set(df.pdbid))}, abids: {len(set(df.abid))}')

    # get number of different values
    organisms = set(df.organism)
    pdbids = set(df.pdbid)
    abids = set(df.abid)
    print(f'organisms: {len(organisms)} {organisms}')

    # get all types of genes
    for gene in gene_dict:
        for gene_list in df[gene]:
            gene_dict[gene] += gene_list.split(',')
        gene_dict[gene] = set(gene_dict[gene])
        print(f'{gene}_genes: {len(gene_dict[gene])} {gene_dict[gene]}')

    # get all types of errors
    if args.log_path:
        log_df = pd.read_csv(args.log_path, names=LOG_HEADER, sep=' ')
        statuses = set(log_df.status)
        print(f'log_statuses: {statuses}')

    # filter by queries
    if args.filter:
        for filter in args.filters:
            df = df[df[filter[0]] == filter[1]]
    print(f'filter_by_query: {len(df)}')

    # filter by json that dont exist in output
    json_paths = glob.glob(f'{args.json_dir}/????-combined.ALL.json')
    pdbids_json = [x.replace(
        f'{args.json_dir}/', '').replace(f'-combined.ALL.json', '') for x in json_paths]
    df = df[df.pdbid.isin(pdbids_json)]
    print(f'filter_by_jsons: {len(df)}')

    # filter by log status
    if args.log_path:
        log_df = log_df[log_df.status.isin(ACCEPTABLE_STATUSES)]
        pbdids_status = set(log_df.pdbid)
        df = df[df.pdbid.isin(pbdids_status)]
        print(f'filter_by_status: {len(df)}')

    # assert that all pdbids have same genes
    pdbids = set(df.pdbid)
    for pdbid in pdbids:
        filter_df = df[df.pdbid == pdbid]
        for col in ['organism', 'ja', 'jb', 'va', 'vb']:
            if col == 'abid':
                continue
            col_values = set([','.join(sorted(x.split(',')))
                              for x in filter_df[col]])
            if len(col_values) != 1:
                cprint(
                    f'col_values do not match for pdbid: {pdbid} {col} {len(col_values)} {col_values}', color=colors.RED)
                # exit()

    # filter down columns
    df = df[args.columns]
    # rename columns
    if args.rename:
        df = df.rename(columns=args.rename)

    # drop duplicate pdbids
    # drop_df = df.drop_duplicates(subset='pdbid', keep='first')
    # assert len(drop_df) == len(df)

    # short example
    json_data = json.loads(df.to_json(orient=args.orient))
    print(f'json_data: {len(json_data)}')
    json_data = json.loads(df[0:5].to_json(orient=args.orient))
    print(json.dumps(json_data, indent=2))

    # compare to original sabdab file
    og_df1 = pd.read_table(args.og_sabdab_path)
    print(f'SABDAB_1:: pdbids: {len(set(og_df1.pdbid))}')

    print(f'FINAL:: pdbids: {len(set(df.pdbid))}, abids: {len(set(df.abid))}')
    human_df = df[df.organism == 'human_ig']
    print(
        f'FINAL_HUMAN:: pdbids: {len(set(human_df.pdbid))}, abids: {len(set(human_df.abid))}')
    mouse_df = df[df.organism == 'mouse_ig']
    print(
        f'FINAL_MOUSE:: pdbids: {len(set(mouse_df.pdbid))}, abids: {len(set(mouse_df.abid))}')

    # write final output
    if args.output_path:
        df.to_json(args.output_path, orient=args.orient)


if __name__ == "__main__":
    args = parse_args()

    main(args)
