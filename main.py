import os

import numpy as np
import pandas as pd
# os.listdir('combined')

FOLDER = "combined"

config_df = pd.read_csv("config_file.csv")
print(config_df)

df_tuple = []
for row_index, row in config_df.iterrows():
    print(row['source'])
    partial_df = pd.read_excel(f"{FOLDER}/{row['source']}")

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
print(f"When merging all list we obtain {len(merged_df)} rows.")

merged_df_filtered = merged_df.drop_duplicates(subset=['inchikey'])
print(f"After dropping InChikey we are left with {len(merged_df_filtered)} rows.")

merged_df_filtered = merged_df_filtered.drop_duplicates(subset=['cas_number'])
print(f"Next, after dropping CAS number we are left with {len(merged_df_filtered)} rows.")

merged_df_filtered = merged_df_filtered.drop_duplicates(subset=['name_in_list'])
print(f"Next, after dropping 'preferred name' we are left with {len(merged_df_filtered)} rows.")

def to_float(value):
    try:
        ret = float(value)
    except:
        ret = 0
    return ret

merged_df_filtered.fillna("", inplace=True)
merged_df_filtered["approximate_mass"] = merged_df_filtered['monoisotopic_mass'].astype(str) +\
                                         merged_df_filtered['approximate_mass'].astype(str)
merged_df_filtered["approximate_mass"] = merged_df_filtered["approximate_mass"].apply(to_float)


merged_df_filtered.sort_values(by=['approximate_mass'], inplace=True)

merged_df_filtered.to_csv('result1.csv')
