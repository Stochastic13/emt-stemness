import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

helper = __import__('helpers')

df_tables_means = pd.read_csv('figures/tables_numbers/pca_link_strength_means.csv')
df_tables_sds = pd.read_csv('figures/tables_numbers/pca_link_strength_sds.csv')
df_significance = pd.read_csv('figures/tables_numbers/stat_significance.csv')

psfnames = ['base', 'grhl2full', 'grhl2partial', 'ovol', 'nrf2']
oe = ['de10', 'oe0', 'oe10']
paths = ['collated_data/base/oe0/base1_clusData.csv', 'collated_data/GRHL2/full/oe0/grhl2full1_clusData.csv',
         'collated_data/GRHL2/partial/oe0/grhl2partial1_clusData.csv', 'collated_data/OVOL/oe0/ovol1_clusData.csv',
         'collated_data/NRF2/oe0/nrf21_clusData.csv']
dfs = [pd.read_csv(i) for i in paths]

sc1 = ['#7a5195', '#003f5c', '#ef5675', '#ffa600']
sc10 = ['#011627', '#ff9f1c', '#e71d36', '#2ec486']
gr1 = ['#1e3d58', '#057dcd', '#43b0f1']

plt.rc('axes', labelsize=15)
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
plt.rc('legend', fontsize=18)

# link strength
fig, ax = plt.subplots(figsize=(7, 7))
ax.bar([0, 1, 2, 3], np.array(df_tables_means['baseoe0'].to_list()[4:8], dtype=float),
       yerr=df_tables_sds['baseoe0'].to_list()[4:8], capsize=5, ecolor='black', color=gr1[1],
       error_kw={'elinewidth': 2, 'capthick': 2})
helper.significance_annotate((0, 1.5), (1, 1.5), df_significance['p'].to_numpy()[df_significance['params'] == 'base_zeb/u200_asym_e_he'], ax, h=0.2, h2=0.1, h3=0.05)
helper.significance_annotate((2, 5.2), (3, 5.2), df_significance['p'].to_numpy()[df_significance['params'] == 'base_zeb/u200_asym_hm_m'], ax, h=0.2, h2=0.1, h3=0.05)
helper.significance_annotate((0, 5.6), (3, 5.6), df_significance['p'].to_numpy()[df_significance['params'] == 'base_zeb/u200_asym_e_m'], ax, h=0.2, h2=0.1, h3=0.05)
ax.set_xticks([0, 1, 2, 3])
ax.set_xticklabels(['e', 'he', 'hm', 'm'])
ax.set_xlabel('ZEB/u200 Asymmetry')
fig.savefig('figures/plots/zeb_u200_asym.png', dpi=500, bbox_inches='tight')
plt.close()

fig, ax = plt.subplots(figsize=(7, 7))
ax.bar([0, 1, 2, 3], np.array(df_tables_means['baseoe0'].to_list()[8:12], dtype=float),
       yerr=df_tables_sds['baseoe0'].to_list()[8:12], capsize=5, ecolor='black', color=gr1[1],
       error_kw={'elinewidth': 2, 'capthick': 2})
helper.significance_annotate_inverted((0, -1.2), (2, -1.2), df_significance['p'].to_numpy()[df_significance['params'] == 'base_lin28/let7_asym_e_hm'], ax, h=0.2, h2=0.1, h3=0.05)
helper.significance_annotate((1, 4.2), (3, 4.1), df_significance['p'].to_numpy()[df_significance['params'] == 'base_lin28/let7_asym_he_m'], ax, h=0.2, h2=0.1, h3=0.05)
helper.significance_annotate((0, 4.6), (3, 4.5), df_significance['p'].to_numpy()[df_significance['params'] == 'base_lin28/let7_asym_e_m'], ax, h=0.2, h2=0.1, h3=0.05)
ax.set_xticks([0, 1, 2, 3])
ax.set_xticklabels(['e', 'he', 'hm', 'm'])
ax.set_xlabel('LIN28/let7 Asymmetry')
fig.savefig('figures/plots/lin28_let7_asym.png', dpi=500, bbox_inches='tight')
plt.close()

fig, ax = plt.subplots(figsize=(7, 7))
ax.bar([0, 1, 2, 3], np.array(df_tables_means['baseoe0'].to_list()[12:16], dtype=float),
       yerr=df_tables_sds['baseoe0'].to_list()[12:16], capsize=5, ecolor='black', color=gr1[1],
       error_kw={'elinewidth': 2, 'capthick': 2})
helper.significance_annotate((0, 20.5), (1, 20.2), df_significance['p'].to_numpy()[df_significance['params'] == 'base_coupling_e_he'], ax, h=0.1, h2=0.05, h3=0.01)
helper.significance_annotate((2, 19.7), (3, 18.8), df_significance['p'].to_numpy()[df_significance['params'] == 'base_coupling_hm_m'], ax, h=0.1, h2=0.05, h3=0.01)
helper.significance_annotate((1, 19.2), (2, 17.5), df_significance['p'].to_numpy()[df_significance['params'] == 'base_coupling_he_hm'], ax, h=0.1, h2=0.05, h3=0.01)
ax.set_xticks([0, 1, 2, 3])
ax.set_xticklabels(['e', 'he', 'hm', 'm'])
ax.set_xlabel('Total Coupling Strength')
ax.set_ylim(18, 21)
fig.savefig('figures/plots/total_coupling_strength.png', dpi=500, bbox_inches='tight')
plt.close()

for psf, df in zip(psfnames, dfs):
    if not (psf == 'base'):
        continue
    pca_coeff_1 = [float(x) for x in df_tables_means[psf + 'oe0'].to_numpy()[0].split(' ')]
    pca_coeff_2 = [float(x) for x in df_tables_means[psf + 'oe0'].to_numpy()[1].split(' ')]
    pca1 = np.zeros(df.shape[0])
    count = 0
    for c in pca_coeff_1:
        pca1 += c * df[df.columns[3 + count]].to_numpy()
        count += 1
    pca2 = np.zeros(df.shape[0])
    count = 0
    for c in pca_coeff_2:
        pca2 += c * df[df.columns[3 + count]].to_numpy()
        count += 1

    xlims = np.percentile(pca1, [0.1, 99.9])  # filter the outliers for better visualization
    ylims = np.percentile(pca2, [0.1, 99.9])

    fig, ax = plt.subplots(figsize=(8, 7))
    count = 0
    for p in ['e', 'he', 'hm', 'm']:
        ax.scatter(pca1[df['phenotype'] == p], pca2[df['phenotype'] == p], color=sc1[count], s=5, alpha=0.5)
        count += 1
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    ax.set_xticks([round(x, 2) for x in np.linspace(xlims[0], xlims[1], 5)])
    ax.set_yticks([round(x, 2) for x in np.linspace(ylims[0], ylims[1], 5)])
    fig.savefig('figures/plots/' + psf + '_run1_oe0_pca_scatter.png', dpi=500, bbox_inches='tight')
    plt.close()

    # transform dataset to make visualization better (capping all values > +2.5 and < -2.5)
    genes = ['u200', 'ZEB', 'LIN28', 'let7', 'SNAIL', 'NFkB', 'GRHL2', 'OVOL', 'KEAP1', 'ECad', 'NRF2']
    for i in [x for x in df.columns[3:] if x in genes]:
        df.loc[df[i] <= -2.5, i] = -2.5
        df.loc[df[i] >= 2.5, i] = 2.5
    col_vec = np.zeros(df.shape[0], dtype=np.int8)  # to color the rows according to the assigned clusters (KMeans)
    col_vec[df['phenotype'] == 'he'] = 1
    col_vec[df['phenotype'] == 'hm'] = 2
    col_vec[df['phenotype'] == 'm'] = 3
    col_vec = np.array(sc1)[col_vec]

    genes = ['u200', 'ZEB', 'LIN28', 'let7']
    fig = sns.clustermap(df.loc[:, genes], method='ward', col_cluster=False, cmap='mako', figsize=(8, 7),
                         dendrogram_ratio=0.1, tree_kws={'linewidths': 2}, cbar_pos=None, row_colors=col_vec)
    fig.savefig('figures/plots/' + psf + '_run1_heatmap_coreGenes.png', dpi=500, bbox_inches='tight')
    plt.close()
