import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import iqr
import os

# Extra analysis for multi-stability and stemness

# Boilerplate code taken from tables_and_numbers_2.py
paths = ['GRHL2/full/oe0', 'GRHL2/full/oe10', 'GRHL2/full/de10',
         'GRHL2/partial/oe0', 'GRHL2/partial/oe10', 'GRHL2/partial/de10',
         'NRF2/oe0', 'NRF2/oe10', 'NRF2/de10',
         'OVOL/oe0', 'OVOL/oe10', 'OVOL/de10',
         'base/oe0']

raw_data_paths = ['collated_data/' + x for x in paths]
save_paths = ['figures/tables_numbers/additional/'] * 13

df_tables_numbers_means = pd.DataFrame()
df_tables_numbers_sds = pd.DataFrame()
parameters = ['p1_e', 'p1_he', 'p1_hm', 'p1_m',
              'p2_e', 'p2_he', 'p2_hm', 'p2_m']

df_tables_numbers_means['parameter_name'] = parameters
df_tables_numbers_sds['parameter_name'] = parameters

all_runs = [pd.read_csv('collated_data/base/oe0/base' + str(x) + '_collated.csv')['LIN28'] for x in range(1, 6)]
all_runs = [x.to_numpy() for x in all_runs]
avg_median = np.mean([np.median(x) for x in all_runs])
avg_iqr = np.mean([iqr(x) for x in all_runs])
unnormalized_lin28_range = (round(avg_median - avg_iqr, 2), round(avg_median + avg_iqr, 2))
rng = unnormalized_lin28_range[1] - unnormalized_lin28_range[0]
mn = np.mean(unnormalized_lin28_range)
unnormalized_stemness_window = (round(mn - 15 / 100 * rng, 2), round(mn + 15 / 100 * rng, 2))

stem_window_znorm = []
for path_to_folder, save_path in zip(raw_data_paths, save_paths):
    print('\nStarting ' + path_to_folder)
    files = os.listdir(path_to_folder)
    circuit_names = [f.split('_')[0] for f in files if 'oe0Normalized.csv' in f]
    oe_level = path_to_folder.split('/')[-1]
    df_list = []
    for circuit_name in circuit_names:
        df_list.append(pd.read_csv(path_to_folder + '/' + circuit_name + '_clusData.csv'))
        df_list[-1] = df_list[-1].loc[df_list[-1]['num_states'] > 1, :]
        # this is the extra line of code to filter out the monostable solutions
    psfname = circuit_names[0][:-1]

    append_vector_means = []
    append_vector_sds = []

    if oe_level == 'oe0':
        df_list_2 = [pd.read_csv(path_to_folder + '/' + psfname + str(i) + '_collated.csv') for i in range(1, 6)]
        references = [(np.mean(df_list_2[i]['LIN28']), np.std(df_list_2[i]['LIN28'], ddof=1)) for i in range(5)]
        stem_window_znorm = []
        for i in range(5):
            lim1 = (unnormalized_stemness_window[0] - references[i][0]) / references[i][1]
            lim2 = (unnormalized_stemness_window[1] - references[i][0]) / references[i][1]
            stem_window_znorm.append((lim1, lim2))

    # p1 (probability of high-stemness given a phenotype)
    for p in ['e', 'he', 'hm', 'm']:
        prob_vec = []
        for i in range(5):
            df = df_list[i]
            bool_mask = (df['phenotype'].to_numpy() == p)
            lin28_values = df.loc[bool_mask, 'LIN28'].to_numpy()
            hs = np.logical_and(lin28_values >= stem_window_znorm[i][0], lin28_values <= stem_window_znorm[i][1])
            prob_vec.append(np.sum(hs) / np.sum(bool_mask))
        append_vector_means.append(np.mean(prob_vec))
        append_vector_sds.append(np.std(prob_vec, ddof=1))

    # p2 (probability of a phenotype given a high-stemness solution)
    hs_phenotypes = []
    for i in range(5):
        df = df_list[i]
        hs_boolmask = np.logical_and(df['LIN28'] >= stem_window_znorm[i][0], df['LIN28'] <= stem_window_znorm[i][1])
        hs_phenotypes.append(df.loc[hs_boolmask, 'phenotype'].to_list())
    for p in ['e', 'he', 'hm', 'm']:
        prob_vec = [hs_phenotypes[i].count(p) / len(hs_phenotypes[i]) for i in range(5)]
        append_vector_means.append(np.mean(prob_vec))
        append_vector_sds.append(np.std(prob_vec, ddof=1))

    df_tables_numbers_means[psfname + oe_level] = append_vector_means
    df_tables_numbers_sds[psfname + oe_level] = append_vector_sds

df_tables_numbers_means.to_csv(save_paths[0] + '/stemness_means_only_multistable.csv', index=False)
df_tables_numbers_sds.to_csv(save_paths[0] + '/stemness_sds_only_multistable.csv', index=False)


def quick(a):  # to handle rounding below
    if isinstance(a, str):
        return a
    else:
        return str(round(a, 2))

# remove the knock-down circuit (GRHL2-KD)
df_tables_numbers_means = df_tables_numbers_means.drop(
    columns=[x for x in df_tables_numbers_means.columns if 'grhl2partial' in x])
# printing directly in latex-usable format of a table
print(' & '.join(df_tables_numbers_means.columns) + '\\\\')
for row in df_tables_numbers_means.itertuples(index=False):
    print(' & '.join([quick(x) for x in row]) + '\\\\')
