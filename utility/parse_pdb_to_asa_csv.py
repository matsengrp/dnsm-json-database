import sys
import os
import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
max_asa_path = f'{SCRIPT_DIR}/max_asa.csv'
max_asa_df = pd.read_csv(max_asa_path)

ASA_METHODS = ["Tien et al. 2013 (theor.)", "Tien et al. 2013 (emp.)",
               "Miller et al. 1987", "Rose et al. 1985"]
ASA_METHOD = ASA_METHODS[1]


def pdb_to_dssp_csv(pdb_id, input_pdb_path, output_dssp_csv_path, chain_id=None):
    parser = PDBParser()
    structure = parser.get_structure(pdb_id, input_pdb_path)
    model = structure[0]
    rel_dssp = DSSP(model, input_pdb_path, dssp='mkdssp', acc_array='Wilke')
    abs_dssp = DSSP(model, input_pdb_path, dssp='mkdssp', acc_array='Sander')
    rel_dssp_keys = list(rel_dssp.keys())

    filter_keys = rel_dssp_keys
    if chain_id:
        filter_keys = list(filter(lambda x: x[0] == chain_id, rel_dssp_keys))

    full_range = list(range(1, len(filter_keys) + 1))
    protein_sites = [f'{x[1][0]}{x[1][1]}{x[1][2]}'.replace(
        ' ', '') for x in filter_keys]
    sequential_sites = range(1, len(protein_sites)+1)
    site_aas = [rel_dssp[key][1] for key in filter_keys]
    rel_asas = [rel_dssp[key][3] for key in filter_keys]
    chain_ids = [key[0] for key in filter_keys]
    abs_asas = [abs_dssp[key][3] for key in filter_keys]
    max_asas = []
    asas = []

    for i in range(len(filter_keys)):
        site_aa = site_aas[i]
        rel_asa = rel_asas[i]
        max_asa = max_asa_df[max_asa_df['Abbr']
                             == site_aa][ASA_METHOD].values[0]
        max_asas.append(max_asa)
        asa = rel_asa * max_asa
        asas.append(asa)

    df = pd.DataFrame({
        'wildtype': site_aas,
        'site': sequential_sites,
        'mutant': '-',
        'rel_asa_authH': rel_asas,
        'asa_authH': asas,
        'abs_asa_authH': abs_asas,
        'max_asa_authH': max_asas,
        'protein_site': protein_sites,
        'chain_ids': chain_ids
    })
    df.to_csv(output_dssp_csv_path, index=False)

    return df


if __name__ == "__main__":
    print("[BEGIN] pdb_to_dssp_csv")

    if len(sys.argv) not in [4, 5]:
        print(
            f'Error: incorrect number of args (expected [4,5], got {len(sys.argv)})')
        print('usage: <i:pdb_path> <o:asa_csv_path> <pdb_id> <optional:chain_id>')
        exit()

    input_pdb_path = sys.argv[1]
    output_dssp_csv_path = sys.argv[2]
    pdb_id = sys.argv[3]
    chain_id = sys.argv[4]
    if len(sys.argv) >= 5:
        chain_id = sys.argv[4]

    df = pdb_to_dssp_csv(pdb_id, input_pdb_path,
                         output_dssp_csv_path, chain_id)

    print(df)

    print("[END] pdb_to_dssp_csv")
