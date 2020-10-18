import numpy as np
import pandas as pd
from scipy.stats import iqr
import os
from math import log10, floor


def round_sig(x, sig=3):  # helper function to round a number to "sig" significant digits
    x = float(x)
    return round(x, sig - int(floor(log10(abs(x)))) - 1)


#  To compute proportions (and the respective standard deviations) for number of states, and different
#  multi-stable distributions of the phases

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
parameters = ['num_states_1', 'num_states_2', 'num_states_3', 'num_states_4', 'num_states_5', 'num_states_morethan_5',
              'overall_e', 'overall_he', 'overall_hm', 'overall_m',
              'monostable_e', 'monostable_he', 'monostable_hm', 'monostable_m',
              'bistable_e/he', 'bistable_e/hm', 'bistable_e/m', 'bistable_he/m', 'bistable_hm/m', 'bistable_he/hm',
              'tristable_e/he/m', 'tristable_e/hm/m', 'tristable_e/he/hm', 'tristable_he/hm/m',
              'bistable_same_cluster', 'tristable_same_cluster_2', 'tristable_same_cluster_3']

df_tables_numbers_means['parameter_name'] = parameters
df_tables_numbers_sds['parameter_name'] = parameters

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

    # num_states (upto 6, and the remaining ones clubbed into a single bin)
    num_states = [[df_param_list[i]['num_states'].to_list().count(j) / 10000 for j in range(1, 6)] for i in range(5)]
    num_states_mean = [np.mean([x[i] for x in num_states]) for i in range(5)]
    num_states_mean += [np.mean([1 - np.sum(x) for x in num_states])]  # add the >5 value
    num_states_sd = [np.std([x[i] for x in num_states], ddof=1) for i in range(5)]
    num_states_sd += [np.std([1 - np.sum(x) for x in num_states], ddof=1)]  # add the >5 value

    append_vector_means += num_states_mean
    append_vector_sds += num_states_sd

    # overall proportions/probabilities
    for p in ['e', 'he', 'hm', 'm']:
        vals = [df_list[i]['phenotype'].to_list().count(p) / df_list[i].shape[0] for i in range(5)]
        append_vector_means.append(np.mean(vals))
        append_vector_sds.append(np.std(vals, ddof=1))

    # monostable proportions/probabilities
    df_special_list = [x.loc[x['num_states'] == 1, :] for x in df_list]
    for p in ['e', 'he', 'hm', 'm']:
        vals = [df_special_list[i]['phenotype'].to_list().count(p) / df_special_list[i].shape[0] for i in range(5)]
        append_vector_means.append(np.mean(vals))
        append_vector_sds.append(np.std(vals, ddof=1))

    extras_mean = []  # models with more than one solution in the same cluster (removed from calculations otherwise)
    extras_sd = []

    # bistable probabilities
    phens = [x.split('/') for x in ['e/he', 'e/hm', 'e/m', 'he/m', 'hm/m', 'he/hm']]  # all bistable phases
    df_special_list = [x.loc[x['num_states'] == 2, :] for x in df_list]
    adjusted_prob_vec = []  # adjustment after removing the ">1 solutions in the same cluster" systems
    common_cluster_vec = []  # the removed systems
    for df in df_special_list:
        common_cluster = 0
        subns = [0, 0, 0, 0, 0, 0]  # ordered in correspondence with the "phens" list
        totalmodels = len(np.unique(df['model_id'].to_numpy()))
        for m_id in np.unique(df['model_id'].to_numpy()):  # for all model_ids, check the phenotypes of the system
            subdf = df.loc[df['model_id'] == m_id, 'phenotype'].to_numpy()
            if len(np.unique(subdf)) == 1:  # only one unique phenotype among the two states
                common_cluster += 1
                continue
            membership = [(subdf[0] in x, subdf[1] in x) for x in phens]  # check where each phenotype belongs
            assert membership.count((True, True)) == 1  # should match with only one phase
            subns[membership.index((True, True))] += 1
        adjustedsubprob = [subns[i] / np.sum(subns) for i in range(6)]  # the adjusted probabilities for one run
        adjusted_prob_vec.append(adjustedsubprob)
        common_cluster_vec.append(common_cluster / totalmodels)
    adjusted_prob_means = [np.mean([adjusted_prob_vec[j][i] for j in range(5)]) for i in range(6)]
    adjusted_prob_sds = [np.std([adjusted_prob_vec[j][i] for j in range(5)], ddof=1) for i in range(6)]
    extras_mean.append(np.mean(common_cluster_vec))
    extras_sd.append(np.std(common_cluster_vec, ddof=1))

    append_vector_means += adjusted_prob_means
    append_vector_sds += adjusted_prob_sds

    # tristable probabilities
    phens = [x.split('/') for x in ['e/he/m', 'e/hm/m', 'e/he/hm', 'he/hm/m']]  # all tristable phases
    df_special_list = [x.loc[x['num_states'] == 3, :] for x in df_list]
    adjusted_prob_vec = []
    common_cluster_vec_2 = []  # systems with two solutions in one of the phenotypes
    common_cluster_vec_3 = []  # systems with all solutions in one of the phenotypes
    for df in df_special_list:
        common_cluster_2 = 0
        common_cluster_3 = 0
        subns = [0, 0, 0, 0]
        totalmodels = len(np.unique(df['model_id'].to_numpy()))
        for m_id in np.unique(df['model_id'].to_numpy()):
            subdf = df.loc[df['model_id'] == m_id, 'phenotype'].to_numpy()
            if len(np.unique(subdf)) == 1:
                common_cluster_3 += 1
                continue
            elif len(np.unique(subdf)) == 2:
                common_cluster_2 += 1
                continue
            membership = [(subdf[0] in x, subdf[1] in x, subdf[2] in x) for x in phens]
            assert membership.count((True, True, True)) == 1
            subns[membership.index((True, True, True))] += 1
        adjustedsubprob = [subns[i] / np.sum(subns) for i in range(4)]
        adjusted_prob_vec.append(adjustedsubprob)
        common_cluster_vec_2.append(common_cluster_2 / totalmodels)
        common_cluster_vec_3.append(common_cluster_3 / totalmodels)
    adjusted_prob_means = [np.mean([adjusted_prob_vec[j][i] for j in range(5)]) for i in range(4)]
    adjusted_prob_sds = [np.std([adjusted_prob_vec[j][i] for j in range(5)], ddof=1) for i in range(4)]
    extras_mean.append(np.mean(common_cluster_vec_2))
    extras_sd.append(np.std(common_cluster_vec_2, ddof=1))
    extras_mean.append(np.mean(common_cluster_vec_3))
    extras_sd.append(np.std(common_cluster_vec_3, ddof=1))

    append_vector_means += adjusted_prob_means
    append_vector_sds += adjusted_prob_sds

    append_vector_means += extras_mean
    append_vector_sds += extras_sd

    df_tables_numbers_means[psfname + oe_level] = append_vector_means
    df_tables_numbers_sds[psfname + oe_level] = append_vector_sds

    if psfname == 'base':  # get replicate-wise ZEB/LIN28 means to illustrate replicate stability
        # the output is directly in form of a latex-table format making it easier to compile into a table
        count = 0
        for df in df_list:
            zeb_e = np.median(df.loc[df['phenotype'] == 'e', 'ZEB']), iqr(df.loc[df['phenotype'] == 'e', 'ZEB'])
            zeb_he = np.median(df.loc[df['phenotype'] == 'he', 'ZEB']), iqr(df.loc[df['phenotype'] == 'he', 'ZEB'])
            zeb_hm = np.median(df.loc[df['phenotype'] == 'hm', 'ZEB']), iqr(df.loc[df['phenotype'] == 'hm', 'ZEB'])
            zeb_m = np.median(df.loc[df['phenotype'] == 'm', 'ZEB']), iqr(df.loc[df['phenotype'] == 'm', 'ZEB'])
            zeb_list = [zeb_e, zeb_he, zeb_hm, zeb_m]
            zeb_list = [str(round_sig(x[0])) + ' (' + str(round_sig(x[1])) + ')' for x in zeb_list]
            print(str(count + 1) + ' (ZEB) & ' + ' & '.join(zeb_list) + '\\\\')
            count += 1
            print('\\hline')
        print('\\hline\\hline')
        count = 0
        for df in df_list:
            lin28_e = np.median(df.loc[df['phenotype'] == 'e', 'LIN28']), iqr(df.loc[df['phenotype'] == 'e', 'LIN28'])
            lin28_he = np.median(df.loc[df['phenotype'] == 'he', 'LIN28']), iqr(df.loc[df['phenotype'] == 'he', 'LIN28'])
            lin28_hm = np.median(df.loc[df['phenotype'] == 'hm', 'LIN28']), iqr(df.loc[df['phenotype'] == 'hm', 'LIN28'])
            lin28_m = np.median(df.loc[df['phenotype'] == 'm', 'LIN28']), iqr(df.loc[df['phenotype'] == 'm', 'LIN28'])
            lin28_list = [lin28_e, lin28_he, lin28_hm, lin28_m]
            lin28_list = [str(round_sig(x[0])) + ' (' + str(round_sig(x[1])) + ')' for x in lin28_list]
            print(str(count + 1) + ' (LIN28) & ' + ' & '.join(lin28_list) + '\\\\')
            count += 1
            print('\\hline')

df_tables_numbers_means.to_csv(save_paths[0] + '/base_dists_means.csv', index=False)
df_tables_numbers_sds.to_csv(save_paths[0] + '/base_dists_sds.csv', index=False)
