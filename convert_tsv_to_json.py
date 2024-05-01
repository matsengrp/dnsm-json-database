import os
import sys
import pandas as pd
import argparse
import pprint as pp
import json

parser = argparse.ArgumentParser(
    description="Converts a tsv file to a json file.")
parser.add_argument('--input_file', '-i')
parser.add_argument('--output_file', '-o')
parser.add_argument(
    '--column', '-c', help='list of fields to output. default: ALL')
parser.add_argument(
    '--filter', '-f', help='list of filter terms. e.g, organism=human')
args = parser.parse_args()

input = 'dnsm-json-database/sabdab_summary_2024-01-26_abid_info.tsv'
columns = ['pdbid', 'abid', 'ja', 'jb', 'va', 'vb']

df = pd.read_csv(args.input_file, sep='\t')
if args.filter:
    filters = [x.split("=") for x in args.filter.split(",")]
    for filter in filters:
        df = df[df[filter[0]] == filter[1]]
if args.column:
    columns = args.column.split(",")
df = df[columns]

json_data = json.loads(df[1:5].to_json())
print(json.dumps(json_data, indent=2))

if args.output_file:
    df.to_json(args.output_file)
