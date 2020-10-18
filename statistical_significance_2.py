import pandas as pd
import numpy as np
from scipy.stats import ttest_ind_from_stats
import statsmodels.stats.multitest as mt

helpers = __import__('helpers')

# Includes the significance testing done in the PSF circuits

df_dists_means = pd.read_csv('figures/tables_numbers/base_dists_means.csv')
df_dists_sds = pd.read_csv('figures/tables_numbers/base_dists_sds.csv')
df_stemness_means = pd.read_csv('figures/tables_numbers/stemness_means.csv')
df_stemness_sds = pd.read_csv('figures/tables_numbers/stemness_sds.csv')

parameters = ['monostable_e', 'monostable_he', 'monostable_hm', 'monostable_m',
              'bistable_e/he', 'bistable_e/hm', 'bistable_e/m', 'bistable_he/m', 'bistable_hm/m', 'bistable_he/hm',
              'tristable_e/he/m', 'tristable_e/hm/m', 'tristable_e/he/hm', 'tristable_he/hm/m',
              'p1_e', 'p1_he', 'p1_hm', 'p1_m',
              'p2_e', 'p2_he', 'p2_hm', 'p2_m']

df_significance = pd.DataFrame({'parameters': parameters})

for circ_name in ['grhl2full', 'grhl2partial', 'ovol', 'nrf2']:
    pvalues = []
    effect_sizes = []
    for i in parameters[:4]:
        m1 = df_dists_means.loc[df_dists_means['parameter_name'] == i, circ_name + 'de10'].to_numpy(dtype=float)[0]
        sd1 = df_dists_sds.loc[df_dists_sds['parameter_name'] == i, circ_name + 'de10'].to_numpy(dtype=float)[0]
        m2 = df_dists_means.loc[df_dists_means['parameter_name'] == i, circ_name + 'oe10'].to_numpy(dtype=float)[0]
        sd2 = df_dists_sds.loc[df_dists_sds['parameter_name'] == i, circ_name + 'oe10'].to_numpy(dtype=float)[0]
        t, p = ttest_ind_from_stats(m1, sd1, 5, m2, sd2, 5, equal_var=False)
        pvalues.append(p)
        effect_sizes.append(m2 / m1)
    for i in parameters[4:10]:
        m1 = df_dists_means.loc[df_dists_means['parameter_name'] == i, circ_name + 'de10'].to_numpy(dtype=float)[0]
        sd1 = df_dists_sds.loc[df_dists_sds['parameter_name'] == i, circ_name + 'de10'].to_numpy(dtype=float)[0]
        m2 = df_dists_means.loc[df_dists_means['parameter_name'] == i, circ_name + 'oe10'].to_numpy(dtype=float)[0]
        sd2 = df_dists_sds.loc[df_dists_sds['parameter_name'] == i, circ_name + 'oe10'].to_numpy(dtype=float)[0]
        t, p = ttest_ind_from_stats(m1, sd1, 5, m2, sd2, 5, equal_var=False)
        pvalues.append(p)
        effect_sizes.append(m2 / m1)
    for i in parameters[10:14]:
        m1 = df_dists_means.loc[df_dists_means['parameter_name'] == i, circ_name + 'de10'].to_numpy(dtype=float)[0]
        sd1 = df_dists_sds.loc[df_dists_sds['parameter_name'] == i, circ_name + 'de10'].to_numpy(dtype=float)[0]
        m2 = df_dists_means.loc[df_dists_means['parameter_name'] == i, circ_name + 'oe10'].to_numpy(dtype=float)[0]
        sd2 = df_dists_sds.loc[df_dists_sds['parameter_name'] == i, circ_name + 'oe10'].to_numpy(dtype=float)[0]
        t, p = ttest_ind_from_stats(m1, sd1, 5, m2, sd2, 5, equal_var=False)
        pvalues.append(p)
        effect_sizes.append(m2 / m1)
    for i in parameters[14:18]:
        m1 = df_stemness_means.loc[df_stemness_means['parameter_name'] == i, circ_name + 'de10'].to_numpy(dtype=float)[
            0]
        sd1 = df_stemness_sds.loc[df_stemness_sds['parameter_name'] == i, circ_name + 'de10'].to_numpy(dtype=float)[0]
        m2 = df_stemness_means.loc[df_stemness_means['parameter_name'] == i, circ_name + 'oe10'].to_numpy(dtype=float)[
            0]
        sd2 = df_stemness_sds.loc[df_stemness_sds['parameter_name'] == i, circ_name + 'oe10'].to_numpy(dtype=float)[0]
        t, p = ttest_ind_from_stats(m1, sd1, 5, m2, sd2, 5, equal_var=False)
        pvalues.append(p)
        effect_sizes.append(m2 / m1)
    for i in parameters[18:]:
        m1 = df_stemness_means.loc[df_stemness_means['parameter_name'] == i, circ_name + 'de10'].to_numpy(dtype=float)[
            0]
        sd1 = df_stemness_sds.loc[df_stemness_sds['parameter_name'] == i, circ_name + 'de10'].to_numpy(dtype=float)[0]
        m2 = df_stemness_means.loc[df_stemness_means['parameter_name'] == i, circ_name + 'oe10'].to_numpy(dtype=float)[
            0]
        sd2 = df_stemness_sds.loc[df_stemness_sds['parameter_name'] == i, circ_name + 'oe10'].to_numpy(dtype=float)[0]
        t, p = ttest_ind_from_stats(m1, sd1, 5, m2, sd2, 5, equal_var=False)
        pvalues.append(p)
        effect_sizes.append(m2 / m1)
    pvalues = np.array(pvalues)
    pvalues[np.array([0, 14, 18])] = mt.multipletests(pvalues[np.array([0, 14, 18])], method='h')[1].tolist()
    pvalues[np.array([1, 15, 19])] = mt.multipletests(pvalues[np.array([1, 15, 19])], method='h')[1].tolist()
    pvalues[np.array([2, 16, 20])] = mt.multipletests(pvalues[np.array([2, 16, 20])], method='h')[1].tolist()
    pvalues[np.array([3, 17, 21])] = mt.multipletests(pvalues[np.array([3, 17, 21])], method='h')[1].tolist()
    df_significance[circ_name + '_p'] = pvalues
    df_significance[circ_name + '_es'] = effect_sizes

df_significance.to_csv('figures/tables_numbers/stat_significance2.csv', index=False)
