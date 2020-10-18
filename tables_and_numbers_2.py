import numpy as np
import pandas as pd
import os
from scipy.stats import iqr

#  To compute proportions (and the respective standard deviations) for multiple stemness probabilities for the different
#  phenotypes

paths = ['GRHL2/full/oe0', 'GRHL2/full/oe10', 'GRHL2/full/de10',
         'GRHL2/partial/oe0', 'GRHL2/partial/oe10', 'GRHL2/partial/de10',
         'NRF2/oe0', 'NRF2/oe10', 'NRF2/de10',
         'OVOL/oe0', 'OVOL/oe10', 'OVOL/de10',
         'base/oe0']

raw_data_paths = ['collated_data/' + x for x in paths]  # where the csv files are stored
save_paths = ['figures/tables_numbers/'] * 13  # where the plots will be stored

df_tables_numbers_means = pd.DataFrame()  # to store all the calculated proportions (means) to later use for plotting
df_tables_numbers_sds = pd.DataFrame()  # to store all the calculated proportions (sds) to later use for plotting
# the columns are the different circuits and oe-levels
# the rows are different parameters
parameters = ['p1_e', 'p1_he', 'p1_hm', 'p1_m',
              'p2_e', 'p2_he', 'p2_hm', 'p2_m',
              'p1_hyb', 'p1_non_hyb',
              'p2_hyb', 'p2_non_hyb']

df_tables_numbers_means['parameter_name'] = parameters
df_tables_numbers_sds['parameter_name'] = parameters

# creating a stemness window
all_runs = [pd.read_csv('collated_data/base/oe0/base' + str(x) + '_collated.csv')['LIN28'] for x in range(1, 6)]
all_runs = [x.to_numpy() for x in all_runs]  # only for uniformity and clarity
# load all the un-normalized base runs (only LIN28 values)
avg_median = np.mean([np.median(x) for x in all_runs])  # median averaged across runs
avg_iqr = np.mean([iqr(x) for x in all_runs])  # IQR averaged across runs
unnormalized_lin28_range = (round(avg_median - avg_iqr, 2), round(avg_median + avg_iqr, 2))  # "biological LIN28 range"
# selecting middle 30% of the range as the stemness window
rng = unnormalized_lin28_range[1] - unnormalized_lin28_range[0]
mn = np.mean(unnormalized_lin28_range)  # Same as the avg_median above, with difference only due to rounding
unnormalized_stemness_window = (round(mn - 15 / 100 * rng, 2), round(mn + 15 / 100 * rng, 2))
print('Unnormalized LIN28 range: ', unnormalized_lin28_range)
print('Unnormalized stemness window: ', unnormalized_stemness_window)

stem_window_znorm = []  # ignore the declaration (to avoid IDE suggestion flags *sigh*)
for path_to_folder, save_path in zip(raw_data_paths, save_paths):
    print('\nStarting ' + path_to_folder)
    files = os.listdir(path_to_folder)
    circuit_names = [f.split('_')[0] for f in files if 'oe0Normalized.csv' in f]  # names of the runs in the folder
    oe_level = path_to_folder.split('/')[-1]  # the over/down-expression level (from among oe0, oe10, de10)
    # assuming no name has '_' in it
    df_list = []  # stores the DataFrames for all the runs/replicates
    df_param_list = []
    for circuit_name in circuit_names:
        df_list.append(pd.read_csv(path_to_folder + '/' + circuit_name + '_clusData.csv'))
        df_param_list.append(pd.read_csv(path_to_folder + '/' + circuit_name + '_parameters_formatted.csv'))
    psfname = circuit_names[0][:-1]  # from among ovol, grhl2full, grhl2partial, base, nrf2

    append_vector_means = []  # to store all parameters and then finally add it to the DataFrame
    append_vector_sds = []

    # calculating reference mean/sd for Z transformation to get the circuit specific normalized stemness window
    if oe_level == 'oe0':
        df_list_2 = [pd.read_csv(path_to_folder + '/' + psfname + str(i) + '_collated.csv') for i in range(1, 6)]
        # load the un-normalized data frames to get the transformation values
        references = [(np.mean(df_list_2[i]['LIN28']), np.std(df_list_2[i]['LIN28'], ddof=1)) for i in range(5)]
        stem_window_znorm = []  # Z-normalized stemness window for each of the runs
        for i in range(5):
            lim1 = (unnormalized_stemness_window[0] - references[i][0]) / references[i][1]
            lim2 = (unnormalized_stemness_window[1] - references[i][0]) / references[i][1]
            stem_window_znorm.append((lim1, lim2))
        # each of the replicates of the reference (oe0) is considered separately. For all other expression levels,
        # replicate number x is compared with the stemness window normalized according to replicate number x of the
        # reference set
        print(psfname + ' normalized stemness windows:')
        print(stem_window_znorm)  # to get an approximate idea of the normalized stemness window

    # p1 (probability of high-stemness given a phenotype)
    for p in ['e', 'he', 'hm', 'm']:
        prob_vec = []
        for i in range(5):
            df = df_list[i]
            bool_mask = (df['phenotype'].to_numpy() == p)  # boolean mask. True for all p-phenotype solutions
            lin28_values = df.loc[bool_mask, 'LIN28'].to_numpy()
            hs = np.logical_and(lin28_values >= stem_window_znorm[i][0], lin28_values <= stem_window_znorm[i][1])
            # hs is True if the solution lies int he stemness window
            prob_vec.append(np.sum(hs) / np.sum(bool_mask))
        append_vector_means.append(np.mean(prob_vec))
        append_vector_sds.append(np.std(prob_vec, ddof=1))

    # p2 (probability of a phenotype given a high-stemness solution)
    hs_phenotypes = []  # hs = high_stemness
    for i in range(5):
        df = df_list[i]
        hs_boolmask = np.logical_and(df['LIN28'] >= stem_window_znorm[i][0], df['LIN28'] <= stem_window_znorm[i][1])
        hs_phenotypes.append(df.loc[hs_boolmask, 'phenotype'].to_list())
    for p in ['e', 'he', 'hm', 'm']:
        prob_vec = [hs_phenotypes[i].count(p) / len(hs_phenotypes[i]) for i in range(5)]
        append_vector_means.append(np.mean(prob_vec))
        append_vector_sds.append(np.std(prob_vec, ddof=1))

    # p1 (hybrid vs non hybrid)
    for p in ['he/hm', 'e/m']:
        prob_vec = []
        for i in range(5):
            df = df_list[i]
            bool_mask = (df['phenotype'].isin(p.split('/')))
            lin28_values = df.loc[bool_mask, 'LIN28'].to_numpy()
            hs = np.logical_and(lin28_values >= stem_window_znorm[i][0], lin28_values <= stem_window_znorm[i][1])
            prob_vec.append(np.sum(hs) / np.sum(bool_mask))
        append_vector_means.append(np.mean(prob_vec))
        append_vector_sds.append(np.std(prob_vec, ddof=1))

    # p2 (hybrid vs non-hybrod)
    hs_phenotypes = []  # hs = high_stemness
    for i in range(5):
        df = df_list[i]
        hs_boolmask = np.logical_and(df['LIN28'] >= stem_window_znorm[i][0], df['LIN28'] <= stem_window_znorm[i][1])
        hs_phenotypes.append(df.loc[hs_boolmask, 'phenotype'].to_list())
    for p in ['he/hm', 'e/m']:
        val = [hs_phenotypes[i].count(p.split('/')[0]) + hs_phenotypes[i].count(p.split('/')[1]) for i in range(5)]
        prob_vec = [val[i] / len(hs_phenotypes[i]) for i in range(5)]
        append_vector_means.append(np.mean(prob_vec))
        append_vector_sds.append(np.std(prob_vec, ddof=1))

    df_tables_numbers_means[psfname + oe_level] = append_vector_means
    df_tables_numbers_sds[psfname + oe_level] = append_vector_sds

df_tables_numbers_means.to_csv(save_paths[0] + '/stemness_means.csv', index=False)
df_tables_numbers_sds.to_csv(save_paths[0] + '/stemness_sds.csv', index=False)
