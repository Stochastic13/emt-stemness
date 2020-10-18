import pandas as pd
import os
helpers = __import__('helpers')

# Normalize the collated_files according to the oe0 version in all cases except uncoupled


paths = ['GRHL2/full/oe0', 'GRHL2/full/oe10', 'GRHL2/full/de10',
         'GRHL2/partial/oe0', 'GRHL2/partial/oe10', 'GRHL2/partial/de10',
         'NRF2/oe0', 'NRF2/oe10', 'NRF2/de10',
         'OVOL/oe0', 'OVOL/oe10', 'OVOL/de10',
         'base/oe0']

raw_data_paths = ['collated_data/' + x for x in paths]  # where the csv files are stored
save_paths = ['collated_data/' + x for x in paths]  # where the new normalized DataFrames are stored

reference_df = []  # to store the DataFrames corresponding to oe0 runs of the currently considered circuit
for path_to_folder, save_path in zip(raw_data_paths, save_paths):
    print('Starting ' + path_to_folder)
    files = os.listdir(path_to_folder)
    circuit_names = [f.split('_')[0] for f in files if 'collated.csv' in f]  # names of the replicates in the folder
    circuit_names = sorted(circuit_names)  # arrange in alphabetical order
    oe_level = path_to_folder.split('/')[-1]  # the over/down-expression level (from among oe0, oe10, de10)
    reference = False
    if oe_level == 'oe0':  # set the new reference DataFrame
        reference_df = []  # reset the reference_df list
        reference = True
    # assuming no name has '_' in it
    for circuit_name in circuit_names:
        df = pd.read_csv(path_to_folder + '/' + circuit_name + '_collated.csv')
        if reference:
            reference_df.append(df)
        run_index = int(circuit_name[-1])  # last letter in the circuit name, the index of the run/replicate
        df2 = helpers.column_normalize_foo(df, reference_df[run_index - 1])
        df2.to_csv(save_path + '/' + circuit_name + '_oe0Normalized.csv', index=False)
