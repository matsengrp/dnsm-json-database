import os
import sys
import pandas as pd
import argparse
import pprint as pp
import json

parser = argparse.ArgumentParser(
    description='Converts a tsv file to a json file.')
parser.add_argument('--input_file', '-i')
parser.add_argument('--output_file', '-o')
parser.add_argument(
    '--column', '-c', help='list of fields to output. default: ALL')
parser.add_argument(
    '--filter', '-f', help='list of filter queries. e.g, organism=human_ig')
parser.add_argument(
    '--rename', '-r', help='list of columns to rename. e.g, va=vh')
parser.add_argument('--orient', help='orientation of json output.')
args = parser.parse_args()

# defaults
input = 'dnsm-json-database/sabdab_summary_2024-01-26_abid_info.tsv'
columns = ['pdbid', 'abid', 'organism', 'ja', 'jb', 'va', 'vb']
orient = 'records'  # split, records, index, columns, values

df = pd.read_csv(args.input_file, sep='\t')
if args.filter:
    filters = [x.split('=') for x in args.filter.split(',')]
    for filter in filters:
        df = df[df[filter[0]] == filter[1]]
if args.column:
    columns = args.column.split(',')
df = df[columns]
if args.rename:
    renames = [x.split('=') for x in args.rename.split(',')]
    rename_dict = {}
    for rename in renames:
        rename_dict[rename[0]] = rename[1]
    print('rename:', rename_dict)
    df = df.rename(columns=rename_dict)

json_data = json.loads(df[1:5].to_json(orient=orient))
print(json.dumps(json_data, indent=2))

if args.output_file:
    df.to_json(args.output_file, orient=orient)
