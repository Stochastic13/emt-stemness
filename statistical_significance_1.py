import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, ttest_ind, ttest_ind_from_stats
import statsmodels.stats.multitest as mt

helpers = __import__('helpers')
# Includes the significance testing done for the plots in the base circuit

df_dists_means = pd.read_csv('figures/tables_numbers/base_dists_means.csv')
df_dists_sds = pd.read_csv('figures/tables_numbers/base_dists_sds.csv')
df_stemness_means = pd.read_csv('figures/tables_numbers/stemness_means.csv')
df_stemness_sds = pd.read_csv('figures/tables_numbers/stemness_sds.csv')
df_links_means = pd.read_csv('figures/tables_numbers/pca_link_strength_means.csv')
df_links_sds = pd.read_csv('figures/tables_numbers/pca_link_strength_sds.csv')

param_names = []
effect_sizes = []
effect_sizes2 = []
pvalues = []

dfbase = pd.read_csv('collated_data/base/oe0/base1_clusData.csv')

# significant differences of gene levels in clusters (base)
phen = ['e', 'he', 'hm', 'm']
for grp1 in range(4):
    for grp2 in range(grp1 + 1, 4):
        subpvalues = []  # pvalues of each group to be separately corrected for multiple testing
        for g in ['ZEB', 'LIN28', 'u200', 'let7', 'SNAIL', 'NFkB']:
            e = dfbase.loc[dfbase['phenotype'] == 'e', g].to_numpy()
            he = dfbase.loc[dfbase['phenotype'] == 'he', g].to_numpy()
            hm = dfbase.loc[dfbase['phenotype'] == 'hm', g].to_numpy()
            m = dfbase.loc[dfbase['phenotype'] == 'm', g].to_numpy()
            groups = {'e': e, 'he': he, 'hm': hm, 'm': m}
            u, p = mannwhitneyu(groups[phen[grp1]], groups[phen[grp2]], use_continuity=True, alternative='two-sided')
            param_names.append('base_' + g + '_' + phen[grp1] + '_' + phen[grp2])
            effect_sizes.append(np.median(groups[phen[grp1]]) / np.median(groups[phen[grp2]]))
            effect_sizes2.append(np.abs(np.median(groups[phen[grp1]]) - np.median(groups[phen[grp2]])))
            subpvalues.append(p)
        pvalues += mt.multipletests(subpvalues, method='h')[1].tolist()  # Holm aka Holm-Bonferroni correction

# significance cross gene
subpvalues = []
for phen in ['e', 'he', 'hm', 'm']:
    u200 = dfbase.loc[dfbase['phenotype'] == phen, 'u200'].to_numpy()
    zeb = dfbase.loc[dfbase['phenotype'] == phen, 'ZEB'].to_numpy()
    lin28 = dfbase.loc[dfbase['phenotype'] == phen, 'LIN28'].to_numpy()
    let7 = dfbase.loc[dfbase['phenotype'] == phen, 'let7'].to_numpy()
    u, p = mannwhitneyu(u200, zeb, use_continuity=True, alternative='two-sided')
    param_names.append('base_' + phen + '_u200_zeb')
    effect_sizes.append(np.median(u200) / np.median(zeb))
    effect_sizes2.append(np.abs(np.median(u200) - np.median(zeb)))
    subpvalues.append(p)
    u, p = mannwhitneyu(lin28, let7, use_continuity=True, alternative='two-sided')
    param_names.append('base_' + phen + '_lin28_let7')
    effect_sizes.append(np.median(lin28) / np.median(let7))
    effect_sizes2.append(np.abs(np.median(lin28) - np.median(let7)))
    subpvalues.append(p)
pvalues += subpvalues
subpvalues = []

param_names += ['base_monostable_e_m', 'base_p1_e_m', 'base_p2_e_m', 'base_zeb/u200_asym_e_m',
                'base_lin28/let7_asym_e_m', 'base_monostable_he_hm', 'base_p1_he_hm', 'base_p2_he_hm',
                'base_coupling_he_hm', 'base_p1_e_he', 'base_p2_e_he', 'base_zeb/u200_asym_e_he', 'base_coupling_e_he',
                'base_p1_hm_m', 'base_p2_hm_m', 'base_zeb/u200_asym_hm_m', 'base_coupling_hm_m']
# comparisons between e and m
m1 = df_dists_means.loc[df_dists_means['parameter_name'] == 'monostable_e', 'baseoe0'].to_numpy(dtype=float)[0]
sd1 = df_dists_sds.loc[df_dists_sds['parameter_name'] == 'monostable_e', 'baseoe0'].to_numpy(dtype=float)[0]
m2 = df_dists_means.loc[df_dists_means['parameter_name'] == 'monostable_m', 'baseoe0'].to_numpy(dtype=float)[0]
sd2 = df_dists_sds.loc[df_dists_sds['parameter_name'] == 'monostable_m', 'baseoe0'].to_numpy(dtype=float)[0]
t, p = ttest_ind_from_stats(m1, sd1, 5, m2, sd2, 5, equal_var=False)
subpvalues.append(p)
effect_sizes.append(m1 / m2)
effect_sizes2.append(np.abs(m1 - m2))
m1 = df_stemness_means.loc[df_stemness_means['parameter_name'] == 'p1_e', 'baseoe0'].to_numpy(dtype=float)[0]
sd1 = df_stemness_sds.loc[df_stemness_sds['parameter_name'] == 'p1_e', 'baseoe0'].to_numpy(dtype=float)[0]
m2 = df_stemness_means.loc[df_stemness_means['parameter_name'] == 'p1_m', 'baseoe0'].to_numpy(dtype=float)[0]
sd2 = df_stemness_sds.loc[df_stemness_sds['parameter_name'] == 'p1_m', 'baseoe0'].to_numpy(dtype=float)[0]
t, p = ttest_ind_from_stats(m1, sd1, 5, m2, sd2, 5, equal_var=False)
subpvalues.append(p)
effect_sizes.append(m1 / m2)
effect_sizes2.append(np.abs(m1 - m2))
m1 = df_stemness_means.loc[df_stemness_means['parameter_name'] == 'p2_e', 'baseoe0'].to_numpy(dtype=float)[0]
sd1 = df_stemness_sds.loc[df_stemness_sds['parameter_name'] == 'p2_e', 'baseoe0'].to_numpy(dtype=float)[0]
m2 = df_stemness_means.loc[df_stemness_means['parameter_name'] == 'p2_m', 'baseoe0'].to_numpy(dtype=float)[0]
sd2 = df_stemness_sds.loc[df_stemness_sds['parameter_name'] == 'p2_m', 'baseoe0'].to_numpy(dtype=float)[0]
t, p = ttest_ind_from_stats(m1, sd1, 5, m2, sd2, 5, equal_var=False)
subpvalues.append(p)
effect_sizes.append(m1 / m2)
effect_sizes2.append(np.abs(m1 - m2))
m1 = df_links_means.loc[df_links_means['parameter_name'] == 'zeb_u200_asym_e', 'baseoe0'].to_numpy(dtype=float)[0]
sd1 = df_links_sds.loc[df_links_sds['parameter_name'] == 'zeb_u200_asym_e', 'baseoe0'].to_numpy(dtype=float)[0]
m2 = df_links_means.loc[df_links_means['parameter_name'] == 'zeb_u200_asym_m', 'baseoe0'].to_numpy(dtype=float)[0]
sd2 = df_links_sds.loc[df_links_sds['parameter_name'] == 'zeb_u200_asym_m', 'baseoe0'].to_numpy(dtype=float)[0]
t, p = ttest_ind_from_stats(m1, sd1, 5, m2, sd2, 5, equal_var=False)
subpvalues.append(p)
effect_sizes.append(m1 / m2)
effect_sizes2.append(np.abs(m1 - m2))
m1 = df_links_means.loc[df_links_means['parameter_name'] == 'lin28_let7_asym_e', 'baseoe0'].to_numpy(dtype=float)[0]
sd1 = df_links_sds.loc[df_links_sds['parameter_name'] == 'lin28_let7_asym_e', 'baseoe0'].to_numpy(dtype=float)[0]
m2 = df_links_means.loc[df_links_means['parameter_name'] == 'lin28_let7_asym_m', 'baseoe0'].to_numpy(dtype=float)[0]
sd2 = df_links_sds.loc[df_links_sds['parameter_name'] == 'lin28_let7_asym_m', 'baseoe0'].to_numpy(dtype=float)[0]
t, p = ttest_ind_from_stats(m1, sd1, 5, m2, sd2, 5, equal_var=False)
subpvalues.append(p)
effect_sizes.append(m1 / m2)
effect_sizes2.append(np.abs(m1 - m2))
pvalues += mt.multipletests(subpvalues, method='h')[1].tolist()
subpvalues = []

# comparisons between he and hm
m1 = df_dists_means.loc[df_dists_means['parameter_name'] == 'monostable_he', 'baseoe0'].to_numpy(dtype=float)[0]
sd1 = df_dists_sds.loc[df_dists_sds['parameter_name'] == 'monostable_he', 'baseoe0'].to_numpy(dtype=float)[0]
m2 = df_dists_means.loc[df_dists_means['parameter_name'] == 'monostable_hm', 'baseoe0'].to_numpy(dtype=float)[0]
sd2 = df_dists_sds.loc[df_dists_sds['parameter_name'] == 'monostable_hm', 'baseoe0'].to_numpy(dtype=float)[0]
t, p = ttest_ind_from_stats(m1, sd1, 5, m2, sd2, 5, equal_var=False)
subpvalues.append(p)
effect_sizes.append(m1 / m2)
effect_sizes2.append(np.abs(m1 - m2))
m1 = df_stemness_means.loc[df_stemness_means['parameter_name'] == 'p1_he', 'baseoe0'].to_numpy(dtype=float)[0]
sd1 = df_stemness_sds.loc[df_stemness_sds['parameter_name'] == 'p1_he', 'baseoe0'].to_numpy(dtype=float)[0]
m2 = df_stemness_means.loc[df_stemness_means['parameter_name'] == 'p1_hm', 'baseoe0'].to_numpy(dtype=float)[0]
sd2 = df_stemness_sds.loc[df_stemness_sds['parameter_name'] == 'p1_hm', 'baseoe0'].to_numpy(dtype=float)[0]
t, p = ttest_ind_from_stats(m1, sd1, 5, m2, sd2, 5, equal_var=False)
subpvalues.append(p)
effect_sizes.append(m1 / m2)
effect_sizes2.append(np.abs(m1 - m2))
m1 = df_stemness_means.loc[df_stemness_means['parameter_name'] == 'p2_he', 'baseoe0'].to_numpy(dtype=float)[0]
sd1 = df_stemness_sds.loc[df_stemness_sds['parameter_name'] == 'p2_he', 'baseoe0'].to_numpy(dtype=float)[0]
m2 = df_stemness_means.loc[df_stemness_means['parameter_name'] == 'p2_hm', 'baseoe0'].to_numpy(dtype=float)[0]
sd2 = df_stemness_sds.loc[df_stemness_sds['parameter_name'] == 'p2_hm', 'baseoe0'].to_numpy(dtype=float)[0]
t, p = ttest_ind_from_stats(m1, sd1, 5, m2, sd2, 5, equal_var=False)
subpvalues.append(p)
effect_sizes.append(m1 / m2)
effect_sizes2.append(np.abs(m1 - m2))
m1 = df_links_means.loc[df_links_means['parameter_name'] == 'total_coupling_he', 'baseoe0'].to_numpy(dtype=float)[0]
sd1 = df_links_sds.loc[df_links_sds['parameter_name'] == 'total_coupling_he', 'baseoe0'].to_numpy(dtype=float)[0]
m2 = df_links_means.loc[df_links_means['parameter_name'] == 'total_coupling_hm', 'baseoe0'].to_numpy(dtype=float)[0]
sd2 = df_links_sds.loc[df_links_sds['parameter_name'] == 'total_coupling_hm', 'baseoe0'].to_numpy(dtype=float)[0]
t, p = ttest_ind_from_stats(m1, sd1, 5, m2, sd2, 5, equal_var=False)
subpvalues.append(p)
effect_sizes.append(m1 / m2)
effect_sizes2.append(np.abs(m1 - m2))
pvalues += mt.multipletests(subpvalues, method='h')[1].tolist()
subpvalues = []

# comparisons between he and e
m1 = df_stemness_means.loc[df_stemness_means['parameter_name'] == 'p1_e', 'baseoe0'].to_numpy(dtype=float)[0]
sd1 = df_stemness_sds.loc[df_stemness_sds['parameter_name'] == 'p1_e', 'baseoe0'].to_numpy(dtype=float)[0]
m2 = df_stemness_means.loc[df_stemness_means['parameter_name'] == 'p1_he', 'baseoe0'].to_numpy(dtype=float)[0]
sd2 = df_stemness_sds.loc[df_stemness_sds['parameter_name'] == 'p1_he', 'baseoe0'].to_numpy(dtype=float)[0]
t, p = ttest_ind_from_stats(m1, sd1, 5, m2, sd2, 5, equal_var=False)
subpvalues.append(p)
effect_sizes.append(m1 / m2)
effect_sizes2.append(np.abs(m1 - m2))
m1 = df_stemness_means.loc[df_stemness_means['parameter_name'] == 'p2_e', 'baseoe0'].to_numpy(dtype=float)[0]
sd1 = df_stemness_sds.loc[df_stemness_sds['parameter_name'] == 'p2_e', 'baseoe0'].to_numpy(dtype=float)[0]
m2 = df_stemness_means.loc[df_stemness_means['parameter_name'] == 'p2_he', 'baseoe0'].to_numpy(dtype=float)[0]
sd2 = df_stemness_sds.loc[df_stemness_sds['parameter_name'] == 'p2_he', 'baseoe0'].to_numpy(dtype=float)[0]
t, p = ttest_ind_from_stats(m1, sd1, 5, m2, sd2, 5, equal_var=False)
subpvalues.append(p)
effect_sizes.append(m1 / m2)
effect_sizes2.append(np.abs(m1 - m2))
m1 = df_links_means.loc[df_links_means['parameter_name'] == 'zeb_u200_asym_e', 'baseoe0'].to_numpy(dtype=float)[0]
sd1 = df_links_sds.loc[df_links_sds['parameter_name'] == 'zeb_u200_asym_e', 'baseoe0'].to_numpy(dtype=float)[0]
m2 = df_links_means.loc[df_links_means['parameter_name'] == 'zeb_u200_asym_he', 'baseoe0'].to_numpy(dtype=float)[0]
sd2 = df_links_sds.loc[df_links_sds['parameter_name'] == 'zeb_u200_asym_he', 'baseoe0'].to_numpy(dtype=float)[0]
t, p = ttest_ind_from_stats(m1, sd1, 5, m2, sd2, 5, equal_var=False)
subpvalues.append(p)
effect_sizes.append(m1 / m2)
effect_sizes2.append(np.abs(m1 - m2))
m1 = df_links_means.loc[df_links_means['parameter_name'] == 'total_coupling_e', 'baseoe0'].to_numpy(dtype=float)[0]
sd1 = df_links_sds.loc[df_links_sds['parameter_name'] == 'total_coupling_e', 'baseoe0'].to_numpy(dtype=float)[0]
m2 = df_links_means.loc[df_links_means['parameter_name'] == 'total_coupling_he', 'baseoe0'].to_numpy(dtype=float)[0]
sd2 = df_links_sds.loc[df_links_sds['parameter_name'] == 'total_coupling_he', 'baseoe0'].to_numpy(dtype=float)[0]
t, p = ttest_ind_from_stats(m1, sd1, 5, m2, sd2, 5, equal_var=False)
subpvalues.append(p)
effect_sizes.append(m1 / m2)
effect_sizes2.append(np.abs(m1 - m2))
pvalues += mt.multipletests(subpvalues, method='h')[1].tolist()
subpvalues = []

# comparisons between hm and m
m1 = df_stemness_means.loc[df_stemness_means['parameter_name'] == 'p1_hm', 'baseoe0'].to_numpy(dtype=float)[0]
sd1 = df_stemness_sds.loc[df_stemness_sds['parameter_name'] == 'p1_hm', 'baseoe0'].to_numpy(dtype=float)[0]
m2 = df_stemness_means.loc[df_stemness_means['parameter_name'] == 'p1_m', 'baseoe0'].to_numpy(dtype=float)[0]
sd2 = df_stemness_sds.loc[df_stemness_sds['parameter_name'] == 'p1_m', 'baseoe0'].to_numpy(dtype=float)[0]
t, p = ttest_ind_from_stats(m1, sd1, 5, m2, sd2, 5, equal_var=False)
subpvalues.append(p)
effect_sizes.append(m1 / m2)
effect_sizes2.append(np.abs(m1 - m2))
m1 = df_stemness_means.loc[df_stemness_means['parameter_name'] == 'p2_hm', 'baseoe0'].to_numpy(dtype=float)[0]
sd1 = df_stemness_sds.loc[df_stemness_sds['parameter_name'] == 'p2_hm', 'baseoe0'].to_numpy(dtype=float)[0]
m2 = df_stemness_means.loc[df_stemness_means['parameter_name'] == 'p2_m', 'baseoe0'].to_numpy(dtype=float)[0]
sd2 = df_stemness_sds.loc[df_stemness_sds['parameter_name'] == 'p2_m', 'baseoe0'].to_numpy(dtype=float)[0]
t, p = ttest_ind_from_stats(m1, sd1, 5, m2, sd2, 5, equal_var=False)
subpvalues.append(p)
effect_sizes.append(m1 / m2)
effect_sizes2.append(np.abs(m1 - m2))
m1 = df_links_means.loc[df_links_means['parameter_name'] == 'zeb_u200_asym_hm', 'baseoe0'].to_numpy(dtype=float)[0]
sd1 = df_links_sds.loc[df_links_sds['parameter_name'] == 'zeb_u200_asym_hm', 'baseoe0'].to_numpy(dtype=float)[0]
m2 = df_links_means.loc[df_links_means['parameter_name'] == 'zeb_u200_asym_m', 'baseoe0'].to_numpy(dtype=float)[0]
sd2 = df_links_sds.loc[df_links_sds['parameter_name'] == 'zeb_u200_asym_m', 'baseoe0'].to_numpy(dtype=float)[0]
t, p = ttest_ind_from_stats(m1, sd1, 5, m2, sd2, 5, equal_var=False)
subpvalues.append(p)
effect_sizes.append(m1 / m2)
effect_sizes2.append(np.abs(m1 - m2))
m1 = df_links_means.loc[df_links_means['parameter_name'] == 'total_coupling_hm', 'baseoe0'].to_numpy(dtype=float)[0]
sd1 = df_links_sds.loc[df_links_sds['parameter_name'] == 'total_coupling_hm', 'baseoe0'].to_numpy(dtype=float)[0]
m2 = df_links_means.loc[df_links_means['parameter_name'] == 'total_coupling_m', 'baseoe0'].to_numpy(dtype=float)[0]
sd2 = df_links_sds.loc[df_links_sds['parameter_name'] == 'total_coupling_m', 'baseoe0'].to_numpy(dtype=float)[0]
t, p = ttest_ind_from_stats(m1, sd1, 5, m2, sd2, 5, equal_var=False)
subpvalues.append(p)
effect_sizes.append(m1 / m2)
effect_sizes2.append(np.abs(m1 - m2))
pvalues += mt.multipletests(subpvalues, method='h')[1].tolist()
subpvalues = []

param_names += ['base_lin28/let7_asym_e_hm', 'base_lin28/let7_asym_he_m']
# comparisons between e/hm and he/m (not needing multiple correction)
m1 = df_links_means.loc[df_links_means['parameter_name'] == 'lin28_let7_asym_e', 'baseoe0'].to_numpy(dtype=float)[0]
sd1 = df_links_sds.loc[df_links_sds['parameter_name'] == 'lin28_let7_asym_e', 'baseoe0'].to_numpy(dtype=float)[0]
m2 = df_links_means.loc[df_links_means['parameter_name'] == 'lin28_let7_asym_hm', 'baseoe0'].to_numpy(dtype=float)[0]
sd2 = df_links_sds.loc[df_links_sds['parameter_name'] == 'lin28_let7_asym_hm', 'baseoe0'].to_numpy(dtype=float)[0]
t, p = ttest_ind_from_stats(m1, sd1, 5, m2, sd2, 5, equal_var=False)
subpvalues.append(p)
effect_sizes.append(m1 / m2)
effect_sizes2.append(np.abs(m1 - m2))
m1 = df_links_means.loc[df_links_means['parameter_name'] == 'lin28_let7_asym_he', 'baseoe0'].to_numpy(dtype=float)[0]
sd1 = df_links_sds.loc[df_links_sds['parameter_name'] == 'lin28_let7_asym_he', 'baseoe0'].to_numpy(dtype=float)[0]
m2 = df_links_means.loc[df_links_means['parameter_name'] == 'lin28_let7_asym_m', 'baseoe0'].to_numpy(dtype=float)[0]
sd2 = df_links_sds.loc[df_links_sds['parameter_name'] == 'lin28_let7_asym_m', 'baseoe0'].to_numpy(dtype=float)[0]
t, p = ttest_ind_from_stats(m1, sd1, 5, m2, sd2, 5, equal_var=False)
subpvalues.append(p)
effect_sizes.append(m1 / m2)
effect_sizes2.append(np.abs(m1 - m2))
pvalues += subpvalues

outputdf = pd.DataFrame(
    data={'params': param_names, 'effect_sizes': effect_sizes, 'p': pvalues, 'effect_sizes2': effect_sizes2})
outputdf.to_csv('figures/tables_numbers/stat_significance.csv', index=False)
