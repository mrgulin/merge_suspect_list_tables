import os.path
import numpy as np
import cirpy  # https://cirpy.readthedocs.io/_/downloads/en/latest/pdf/
import molmass
import pandas as pd
from cirpy import Molecule
from multiprocessing import Pool
import time
from concurrent.futures import ThreadPoolExecutor

FOLDER = "combined"

RESULT_LIST = []

class StatusOptions:
    SUCCESS = "success"
    FAILED = "fail"
    NO_RESULT = 'no result'


def obtain_info_about_molecule(value: str, resolvers=('name_by_opsin', 'name_by_cir')):
    output_dictionary = {'name_in_list': value}
    try:
        mol_obj = Molecule(value, resolvers)

        request_dictionary = {
            "cas_number": "cas",
            "inchikey": 'stdinchikey',
            'inchi': "stdinchi",
            "iupac": 'iupac_name',
            "smiles": "smiles",
            "formula": "formula",
            "approximate_mass": 'mw'
        }

        for table_name, request_name in request_dictionary.items():
            response = mol_obj.__getattribute__(request_name)
            if response is not None:
                output_dictionary[table_name] = response
            else:
                output_dictionary[table_name] = ''

        if output_dictionary["formula"] != "":
            output_dictionary["monoisotopic_mass"] = molmass.Formula(output_dictionary['formula']).isotope.mass
            # Isotope composed of most abundant elemental isotopes
            output_dictionary['status'] = StatusOptions.SUCCESS
        else:
            output_dictionary = {'name_in_list': value, 'status': StatusOptions.NO_RESULT}
    except Exception as e:
        print(e)
        output_dictionary = {'name_in_list': value, 'status': StatusOptions.FAILED}
    return output_dictionary


def resolve_list_of_compounds(input_list: list[str], timeout_sec: int) -> pd.DataFrame:
    # with Pool(20) as p:
    #     result = p.map(obtain_info_about_molecule, input_list)
    # return pd.DataFrame(result, index=input_list)

    final_list = []

    with ThreadPoolExecutor(20) as executor:
        threads = [executor.submit(obtain_info_about_molecule, compound_name) for compound_name in input_list]

        while any([thread.running() for thread in threads]):
            time.sleep(1)
        exceptions = [thread.exception() for thread in threads if thread.exception()]
        print(exceptions)

        final_list = [i.result() for i in threads]
    final_list = [i for i in final_list if i is not None]
    df = pd.DataFrame(final_list)
    df = df.set_index(df['name_in_list'])
    return df


def get_new_excel_file(path: str, row: pd.Series) -> pd.DataFrame:
    rename_dict = {row['internal_index']: 'internal_index',
                   row['cas_number']: 'cas_number_old', row['name_in_list']: 'name_in_list',
                   row['inchikey']: 'inchikey_old', row['inchi']: 'inchi_old', row['iupac']: 'iupac_old',
                   row['smiles']: 'smiles_old', row['formula']: 'formula_old',
                   row['monoisotopic_mass']: 'monoisotopic_mass_old', row['comment']: 'comment_old',
                   row['approximate_mass']: 'approximate_mass_old'}
    rename_dict = {key: value for key, value in rename_dict.items() if not pd.isna(key)}

    df = pd.read_excel(path, sheet_name=0)
    df.rename(columns=rename_dict, inplace=True)
    df['status'] = np.nan
    return df


def handle_inperfect_excel_files(recalculate_on=(StatusOptions.FAILED,), timeout_seconds=600):
    config_df = pd.read_csv("config_fill_files.csv")

    for row_index, row in config_df.iterrows():
        if row["recalculate"] is False:
            continue
        original_file_path = f"{FOLDER}/{row['source']}"
        filled_file_path = f"{FOLDER}/{row['source'][:-5]}_filled.xlsx"
        if os.path.isfile(filled_file_path):
            df = pd.read_excel(filled_file_path)
        else:
            df = get_new_excel_file(original_file_path, row)

        df = df.set_index(df['name_in_list'])
        df = df[~df.index.duplicated(keep='first')]
        unresolved_names = df.loc[df['status'].isin(recalculate_on) | df['status'].isna(),
                                  'name_in_list'].tolist()
        if len(unresolved_names) == 0:
            continue
        resolved_structures_df = resolve_list_of_compounds(unresolved_names, timeout_sec=timeout_seconds)
        resolved_structures_df['cas_number'] = resolved_structures_df['cas_number'].apply(
            lambda x: "\n".join(x) if isinstance(x, list) else x)
        for column in resolved_structures_df.columns:
            if column not in df.columns:
                df[column] = np.nan
        df.update(resolved_structures_df)
        df.to_excel(f"{FOLDER}/{row['source'][:-5]}_filled.xlsx", index=False)


if __name__ == "__main__":
    handle_inperfect_excel_files()
    # obtain_info_about_molecule("Octocrylene")
