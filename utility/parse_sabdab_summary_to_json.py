import os
import sys
import pandas as pd
import argparse
import pprint as pp
import json
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
COLUMNS = ['pdbid', 'abid', 'organism', 'ja', 'jb', 'va', 'vb']
RENAME = {'ja': 'jl', 'jb': 'jh', 'va': 'vl', 'vb': 'vh'}
# orientation: (split | records | index | columns | values)
ORIENT = 'records'

COLUMNS = list_to_string(COLUMNS)
RENAME = dict_to_string(RENAME)

# file data
LOGFILE_HEADER = ['']
ACCEPTABLE_STATUS = ['success', 'error:']


def parse_args():
    parser = argparse.ArgumentParser(
        description='Converts sabdab summary tsv file to a json file.')
    parser.add_argument('--input-path', '-i', help='input sabdab summary file')
    parser.add_argument('--output-path', '-o', help='output json file')
    parser.add_argument('--log-path', '-l', help='logfile ')
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

    # filter by queries
    if args.filter:
        for filter in args.filters:
            df = df[df[filter[0]] == filter[1]]

    # filter down columns
    df = df[args.columns]
    print(f'sabdf: {len(df)}')

    # get all types of genes
    for gene in gene_dict:
        for gene_list in df[gene]:
            gene_dict[gene] += gene_list.split(',')
        gene_dict[gene] = set(gene_dict[gene])
        print(f'gene: {gene} {len(gene_dict[gene])} {gene_dict[gene]}')

    # assert that all pdbids have same genes
    pdbids = set(df.pdbid)
    print(f'pdbids: {len(pdbids)}')
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

    # rename columns
    if args.rename:
        print('rename:', args.rename)
        df = df.rename(columns=args.rename)

    # clean up pdbids fail cases from log
    # if args.log_path:
    log_df = pd.read_csv()
    statuses = set(log_df.status)
    print('log_statuses: {statuses}')

    # drop duplicate pdbids
    df = df.drop_duplicates(subset='pdbid', keep='first')

    # examples
    json_data = json.loads(df.to_json(orient=args.orient))
    print(f'json_data: {len(json_data)}')
    json_data = json.loads(df[0:5].to_json(orient=args.orient))
    print(json.dumps(json_data, indent=2))

    if args.output_file:
        df.to_json(args.output_path, orient=args.orient)


if __name__ == "__main__":
    args = parse_args()

    main(args)
