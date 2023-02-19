import os

import numpy as np
import pandas as pd

import fill_missing_fields

# os.listdir('combined')

FOLDER = "combined"


def merge_compounds():
    config_df = pd.read_csv("config_file.csv")
    print(config_df)

    df_tuple = []
    for row_index, row in config_df.iterrows():
        print(row['source'])
        partial_df = pd.read_excel(f"{FOLDER}/{row['source']}", sheet_name=0)

        rename_dict = {  # row['source']: 'source',
            row['internal_index']: 'internal_index',
            row['cas_number']: 'cas_number', row['name_in_list']: 'name_in_list',
            row['inchikey']: 'inchikey', row['inchi']: 'inchi', row['iupac']: 'iupac',
            row['smiles']: 'smiles', row['formula']: 'formula',
            row['monoisotopic_mass']: 'monoisotopic_mass', row['comment']: 'comment',
            row['approximate_mass']: 'approximate_mass'}
        partial_df['source'] = row['source']
        important_columns = list(rename_dict.keys())

        important_columns_filtered = [col_name for col_name in important_columns if not pd.isna(col_name)]

        partial_df = partial_df[important_columns_filtered]

        partial_df.rename(columns=rename_dict, inplace=True)

        df_tuple.append(partial_df)

    merged_df = pd.concat(df_tuple)
    return merged_df


def deduplicate_df(df: pd.DataFrame, unique_columns=('inchikey', 'cas_number', 'name_in_list')):
    print(f"When merging all list we obtain {len(df)} rows.")

    df_filtered = df
    for column in unique_columns:
        df_filtered = df_filtered[(~df_filtered[column].duplicated()) | (df_filtered[column].isnull())]
        print(f"After dropping {column} we are left with {len(df_filtered)} rows.")

    return df_filtered


def sort_on_mass(df):
    def to_float(value):
        try:
            ret = float(value)
        except ValueError:
            ret = 0
        return ret
    df = df.copy()
    df.fillna("", inplace=True)
    df.loc[df["approximate_mass"] == '', "approximate_mass"] = df.loc[df["approximate_mass"] == '', "monoisotopic_mass"]
    df["approximate_mass"] = df["approximate_mass"].apply(to_float)

    df.sort_values(by=['approximate_mass'], inplace=True)
    return df


def main():
    # fill_missing_fields.handle_inperfect_excel_files()

    df_all = merge_compounds()
    df_filtered = deduplicate_df(df_all)
    df_sorted = sort_on_mass(df_filtered)
    df_sorted.to_csv('result1.csv')


if __name__ == "__main__":
    main()
