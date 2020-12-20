import numpy as np
import pandas as pd
from scipy.stats import iqr
import os

# Extra analysis for multi-stability and stemness

# Boilerplate code taken from tables_and_numbers_2.py
paths = ['base/oe0', 'GRHL2/full/oe0', 'NRF2/oe0', 'OVOL/oe0']

raw_data_paths = ['collated_data/' + x for x in paths]
save_paths = ['figures/tables_numbers/additional/'] * 13

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
    df_param_list = []
    for circuit_name in circuit_names:
        df_list.append(pd.read_csv(path_to_folder + '/' + circuit_name + '_clusData.csv'))
        df_param_list.append(pd.read_csv(path_to_folder + '/' + circuit_name + '_parameters_formatted.csv'))
    psfname = circuit_names[0][:-1]

    if oe_level == 'oe0':
        df_list_2 = [pd.read_csv(path_to_folder + '/' + psfname + str(i) + '_collated.csv') for i in range(1, 6)]
        references = [(np.mean(df_list_2[i]['LIN28']), np.std(df_list_2[i]['LIN28'], ddof=1)) for i in range(5)]
        stem_window_znorm = []
        for i in range(5):
            lim1 = (unnormalized_stemness_window[0] - references[i][0]) / references[i][1]
            lim2 = (unnormalized_stemness_window[1] - references[i][0]) / references[i][1]
            stem_window_znorm.append((lim1, lim2))

    # bistable_phases
    bistable_phases = [('e', 'he'), ('e', 'hm'), ('e', 'm'), ('he', 'hm'), ('he', 'm'), ('hm', 'm')]
    stemness_poss = [(True, True), (True, False), (False, True), (False, False)]
    main_map = np.zeros((6, 4))
    for i in range(5):
        sub_map = np.zeros((6, 4))
        df = df_list[i]
        hs_boolmask = np.logical_and(df['LIN28'] >= stem_window_znorm[i][0], df['LIN28'] <= stem_window_znorm[i][1])
        df['stemness_bool'] = hs_boolmask
        for m_id in df_param_list[i]['model_id']:
            subdf = df.loc[df['model_id'] == m_id, :]
            if subdf.shape[0] == 2:
                phens = tuple(sorted(subdf['phenotype'].tolist()))
                if phens in bistable_phases:
                    val1 = subdf.loc[subdf['phenotype'] == phens[0], 'stemness_bool'].to_numpy()[0]
                    val2 = subdf.loc[subdf['phenotype'] == phens[1], 'stemness_bool'].to_numpy()[0]
                    sub_map[bistable_phases.index(phens), stemness_poss.index((val1, val2))] += 1
        main_map += sub_map / np.sum(sub_map, axis=1)[:, None] * 100  # append percentages
    main_map /= 5  # average across the replicates
    # outputting in latex-usable format
    fstring = '{:^10s}&{:^10s}&{:^10s}&{:^10s}&{:^10s}\\\\'
    print(fstring.format(*(['phase'] + [str(x).replace('alse', '').replace('rue', '') for x in stemness_poss])))
    print('\\hline\\hline')
    for i in range(6):
        print(fstring.format(*(['/'.join(bistable_phases[i])] + [str(round(x, 2)) for x in main_map[i, :].tolist()])))
        print('\\hline')

    tristable_phases = [('e', 'he', 'm'), ('e', 'hm', 'm'), ('he', 'hm', 'm'), ('e', 'he', 'hm')]
    stemness_poss = [(True, True, True), (True, True, False), (True, False, True), (False, True, True),
                     (True, False, False), (False, True, False), (False, False, True), (False, False, False)]
    main_map = np.zeros((4, 8))
    for i in range(5):
        sub_map = np.zeros((4, 8))
        df = df_list[i]
        hs_boolmask = np.logical_and(df['LIN28'] >= stem_window_znorm[i][0], df['LIN28'] <= stem_window_znorm[i][1])
        df['stemness_bool'] = hs_boolmask
        for m_id in df_param_list[i]['model_id']:
            subdf = df.loc[df['model_id'] == m_id, :]
            if subdf.shape[0] == 3:
                phens = tuple(sorted(subdf['phenotype'].tolist()))
                if phens in tristable_phases:
                    val1 = subdf.loc[subdf['phenotype'] == phens[0], 'stemness_bool'].to_numpy()[0]
                    val2 = subdf.loc[subdf['phenotype'] == phens[1], 'stemness_bool'].to_numpy()[0]
                    val3 = subdf.loc[subdf['phenotype'] == phens[2], 'stemness_bool'].to_numpy()[0]
                    sub_map[tristable_phases.index(phens), stemness_poss.index((val1, val2, val3))] += 1
        main_map += sub_map / np.sum(sub_map, axis=1)[:, None] * 100
    main_map /= 5
    fstring = '{:^10s}&{:^13s}&{:^13s}&{:^13s}&{:^13s}&{:^13s}&{:^13s}&{:^13s}&{:^13s}\\\\'
    print(fstring.format(*(['phase'] + [str(x).replace('alse', '').replace('rue', '') for x in stemness_poss])))
    print('\\hline\\hline')
    for i in range(4):
        print(fstring.format(*(['/'.join(tristable_phases[i])] + [str(round(x, 2)) for x in main_map[i, :].tolist()])))
        print('\\hline')
