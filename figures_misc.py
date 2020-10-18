import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans, AgglomerativeClustering
from scipy.cluster import hierarchy

helpers = __import__('helpers')

psfnames = ['base', 'grhl2full', 'grhl2partial', 'ovol', 'nrf2']
oe = ['de10', 'oe0', 'oe10']
paths = ['collated_data/base/oe0/base1_clusData.csv', 'collated_data/GRHL2/full/oe0/grhl2full1_clusData.csv',
         'collated_data/GRHL2/partial/oe0/grhl2partial1_clusData.csv', 'collated_data/OVOL/oe0/ovol1_clusData.csv',
         'collated_data/NRF2/oe0/nrf21_clusData.csv']
dfs = [pd.read_csv(i) for i in paths]
df_significance = pd.read_csv('figures/tables_numbers/stat_significance.csv')

sc1 = ['#7a5195', '#003f5c', '#ef5675', '#ffa600']
sc10 = ['#011627', '#ff9f1c', '#e71d36', '#2ec486']
gr1 = ['#1e3d58', '#057dcd', '#43b0f1']

plt.rc('axes', labelsize=15)
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
plt.rc('legend', fontsize=18)

# shifted hill function plots with varying parameters
fig, ax = plt.subplots(2, 3, figsize=(16, 12))
b0 = 20
b = np.linspace(0, 70, 100)
n = 4
ax[0, 0].plot(b, [helpers.shifted_hill_foo(2, b0, x, n) for x in b], color=gr1[0], linewidth=2, label='$\lambda = 2$')
ax[0, 0].plot(b, [helpers.shifted_hill_foo(4, b0, x, n) for x in b], color=gr1[1], linewidth=2, label='$\lambda = 4$')
ax[0, 0].plot(b, [helpers.shifted_hill_foo(6, b0, x, n) for x in b], color=gr1[2], linewidth=2, label='$\lambda = 6$')
ax[0, 0].legend()

b0 = 20
b = np.linspace(0, 70, 100)
y = 4
ax[0, 1].plot(b, [helpers.shifted_hill_foo(y, b0, x, 2) for x in b], color=gr1[0], linewidth=2, label='$n = 2$')
ax[0, 1].plot(b, [helpers.shifted_hill_foo(y, b0, x, 4) for x in b], color=gr1[1], linewidth=2, label='$n = 4$')
ax[0, 1].plot(b, [helpers.shifted_hill_foo(y, b0, x, 6) for x in b], color=gr1[2], linewidth=2, label='$n = 6$')
ax[0, 1].legend()

n = 4
b = np.linspace(0, 70, 100)
y = 4
ax[0, 2].plot(b, [helpers.shifted_hill_foo(y, 10, x, n) for x in b], color=gr1[0], linewidth=2, label='$\mu = 10$')
ax[0, 2].plot(b, [helpers.shifted_hill_foo(y, 20, x, n) for x in b], color=gr1[1], linewidth=2, label='$\mu = 20$')
ax[0, 2].plot(b, [helpers.shifted_hill_foo(y, 30, x, n) for x in b], color=gr1[2], linewidth=2, label='$\mu = 30$')
ax[0, 2].legend()

b0 = 20
b = np.linspace(0, 70, 100)
n = 4
ax[1, 0].plot(b, [helpers.shifted_hill_foo(0.75, b0, x, n) for x in b], color=gr1[0], linewidth=2,
              label='$\lambda = 0.75$')
ax[1, 0].plot(b, [helpers.shifted_hill_foo(0.5, b0, x, n) for x in b], color=gr1[1], linewidth=2,
              label='$\lambda = 0.5$')
ax[1, 0].plot(b, [helpers.shifted_hill_foo(0.25, b0, x, n) for x in b], color=gr1[2], linewidth=2,
              label='$\lambda = 0.25$')
ax[1, 0].legend()

b0 = 20
b = np.linspace(0, 70, 100)
y = 0.5
ax[1, 1].plot(b, [helpers.shifted_hill_foo(y, b0, x, 2) for x in b], color=gr1[0], linewidth=2, label='$n = 2$')
ax[1, 1].plot(b, [helpers.shifted_hill_foo(y, b0, x, 4) for x in b], color=gr1[1], linewidth=2, label='$n = 4$')
ax[1, 1].plot(b, [helpers.shifted_hill_foo(y, b0, x, 6) for x in b], color=gr1[2], linewidth=2, label='$n = 6$')
ax[1, 1].legend()

n = 4
b = np.linspace(0, 70, 100)
y = 0.5
ax[1, 2].plot(b, [helpers.shifted_hill_foo(y, 10, x, n) for x in b], color=gr1[0], linewidth=2, label='$\mu = 10$')
ax[1, 2].plot(b, [helpers.shifted_hill_foo(y, 20, x, n) for x in b], color=gr1[1], linewidth=2, label='$\mu = 20$')
ax[1, 2].plot(b, [helpers.shifted_hill_foo(y, 30, x, n) for x in b], color=gr1[2], linewidth=2, label='$\mu = 30$')
ax[1, 2].legend()

for i in ax[0, :].ravel():
    i.set_xlabel('node expression levels')
    i.set_ylabel('$H^{S+}$')
for i in ax[1, :].ravel():
    i.set_xlabel('node expression levels')
    i.set_ylabel('$H^{S-}$')
plt.subplots_adjust(wspace=0.3)
fig.savefig('figures/plots/shifted_hill_function.png', dpi=500, bbox_inches='tight')
plt.close()

# base histograms
# base histograms for ZEB, LIN28, u200, ZEB, SNAIL, NFkB
for gene in ['u200', 'ZEB', 'LIN28', 'let7', 'SNAIL', 'NFkB']:
    fig, ax = plt.subplots(figsize=(7, 7))  # stacked plots not used for the report
    ax.hist([dfs[0].loc[dfs[0]['phenotype'] == x, gene].to_numpy() for x in ['e', 'he', 'hm', 'm']], bins=100,
            stacked=True, color=sc1)
    ax.set_xlabel(gene)
    ax.set_ylabel('frequency')
    if gene == 'LIN28':
        # for these numbers, run tables_and_numbers_3.py
        ax.axvline(-0.543, linestyle='--', color='black', linewidth=2)
        ax.axvline(0.525, linestyle='--', color='black', linewidth=2)
    fig.savefig('figures/plots/base_run1_' + gene + '_stacked.png', dpi=500, bbox_inches='tight')
    plt.close()

    fig, ax = plt.subplots(figsize=(7, 7))
    for x in range(4):
        ax.hist(dfs[0].loc[dfs[0]['phenotype'] == ['e', 'he', 'hm', 'm'][x], gene].to_numpy(), bins=100,
                stacked=False, histtype='step', linewidth=3, color=sc1[x])
        ax.hist(dfs[0].loc[dfs[0]['phenotype'] == ['e', 'he', 'hm', 'm'][x], gene].to_numpy(), bins=100, alpha=0.4,
                color=sc1[x])
    ax.set_xlabel(gene)
    ax.set_ylabel('frequency')
    if gene == 'ZEB':  # annotation for significant differences
        helpers.significance_annotate((np.median(dfs[0].loc[dfs[0]['phenotype'] == 'e', gene].to_numpy()), 210),
                                      (np.median(dfs[0].loc[dfs[0]['phenotype'] == 'he', gene].to_numpy()), 200),
                                      df_significance['p'].to_numpy()[df_significance['params'] == 'base_ZEB_e_he'],
                                      ax, h=10, h2=5, h3=2)
        helpers.significance_annotate((np.median(dfs[0].loc[dfs[0]['phenotype'] == 'hm', gene].to_numpy()), 310),
                                      (np.median(dfs[0].loc[dfs[0]['phenotype'] == 'm', gene].to_numpy()), 310),
                                      df_significance['p'].to_numpy()[df_significance['params'] == 'base_ZEB_hm_m'],
                                      ax, h=10, h2=5, h3=2)
    elif gene == 'u200':
        helpers.significance_annotate((np.median(dfs[0].loc[dfs[0]['phenotype'] == 'e', gene].to_numpy()), 200),
                                      (np.median(dfs[0].loc[dfs[0]['phenotype'] == 'he', gene].to_numpy()), 200),
                                      df_significance['p'].to_numpy()[
                                          df_significance['params'] == 'base_u200_e_he'],
                                      ax, h=10, h2=5, h3=2)
        helpers.significance_annotate((np.median(dfs[0].loc[dfs[0]['phenotype'] == 'hm', gene].to_numpy()), 240),
                                      (np.median(dfs[0].loc[dfs[0]['phenotype'] == 'm', gene].to_numpy()), 240),
                                      df_significance['p'].to_numpy()[
                                          df_significance['params'] == 'base_u200_hm_m'],
                                      ax, h=10, h2=5, h3=2)
    elif gene == 'let7':
        helpers.significance_annotate((np.median(dfs[0].loc[dfs[0]['phenotype'] == 'e', gene].to_numpy()), 280),
                                      (np.median(dfs[0].loc[dfs[0]['phenotype'] == 'hm', gene].to_numpy()), 280),
                                      df_significance['p'].to_numpy()[
                                          df_significance['params'] == 'base_let7_e_hm'],
                                      ax, h=10, h2=5, h3=2)
        helpers.significance_annotate((np.median(dfs[0].loc[dfs[0]['phenotype'] == 'he', gene].to_numpy()), 280),
                                      (np.median(dfs[0].loc[dfs[0]['phenotype'] == 'm', gene].to_numpy()), 280),
                                      df_significance['p'].to_numpy()[
                                          df_significance['params'] == 'base_let7_he_m'],
                                      ax, h=10, h2=5, h3=2)
    elif gene == 'LIN28':
        # for these numbers, run tables_and_numbers_3.py
        ax.axvline(-0.543, linestyle='--', color='black', linewidth=2)
        ax.axvline(0.525, linestyle='--', color='black', linewidth=2)
        helpers.significance_annotate((np.median(dfs[0].loc[dfs[0]['phenotype'] == 'e', gene].to_numpy()), 210),
                                      (np.median(dfs[0].loc[dfs[0]['phenotype'] == 'hm', gene].to_numpy()), 220),
                                      df_significance['p'].to_numpy()[
                                          df_significance['params'] == 'base_LIN28_e_hm'],
                                      ax, h=10, h2=5, h3=2)
        helpers.significance_annotate((np.median(dfs[0].loc[dfs[0]['phenotype'] == 'he', gene].to_numpy()), 340),
                                      (np.median(dfs[0].loc[dfs[0]['phenotype'] == 'm', gene].to_numpy()), 340),
                                      df_significance['p'].to_numpy()[
                                          df_significance['params'] == 'base_LIN28_he_m'],
                                      ax, h=10, h2=5, h3=2)

    fig.savefig('figures/plots/base_run1_' + gene + '_overlapping.png', dpi=500, bbox_inches='tight')
    plt.close()

# PSF level histograms
for psf, path in zip(psfnames, paths):
    old_psf = psf
    if psf == 'base':
        continue
    df_oe0 = pd.read_csv(path)
    df_de10 = pd.read_csv(path.replace('oe0', 'de10'))
    df_oe10 = pd.read_csv(path.replace('oe0', 'oe10'))
    if 'grhl2' in psf:  # for both the partial and full circuit
        psf = 'GRHL2'
    else:
        psf = psf.upper()
    fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    for x in range(4):  # overlapping version not used for the report
        ax[0].hist(df_de10.loc[df_de10['phenotype'] == ['e', 'he', 'hm', 'm'][x], psf].to_numpy(), bins=100,
                   stacked=False, histtype='step', linewidth=3, color=sc1[x])
        ax[0].hist(df_de10.loc[df_de10['phenotype'] == ['e', 'he', 'hm', 'm'][x], psf].to_numpy(), bins=100, alpha=0.4,
                   color=sc1[x])
        ax[1].hist(df_oe0.loc[df_oe0['phenotype'] == ['e', 'he', 'hm', 'm'][x], psf].to_numpy(), bins=100,
                   stacked=False, histtype='step', linewidth=3, color=sc1[x])
        ax[1].hist(df_oe0.loc[df_oe0['phenotype'] == ['e', 'he', 'hm', 'm'][x], psf].to_numpy(), bins=100, alpha=0.4,
                   color=sc1[x])
        ax[2].hist(df_oe10.loc[df_oe10['phenotype'] == ['e', 'he', 'hm', 'm'][x], psf].to_numpy(), bins=100,
                   stacked=False, histtype='step', linewidth=3, color=sc1[x])
        ax[2].hist(df_oe10.loc[df_oe10['phenotype'] == ['e', 'he', 'hm', 'm'][x], psf].to_numpy(), bins=100, alpha=0.4,
                   color=sc1[x])
    fig.savefig('figures/plots/' + old_psf + '_psflevel_overlapping.png', dpi=500, bbox_inches='tight')

    fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    ax[0].hist([df_de10.loc[df_de10['phenotype'] == x, psf].to_numpy() for x in ['e', 'he', 'hm', 'm']], bins=100,
               stacked=True, color=sc1)
    ax[1].hist([df_oe0.loc[df_oe0['phenotype'] == x, psf].to_numpy() for x in ['e', 'he', 'hm', 'm']], bins=100,
               stacked=True, color=sc1)
    ax[2].hist([df_oe10.loc[df_oe10['phenotype'] == x, psf].to_numpy() for x in ['e', 'he', 'hm', 'm']], bins=100,
               stacked=True, color=sc1)
    fig.savefig('figures/plots/' + old_psf + '_psflevel_stacked.png', dpi=500, bbox_inches='tight')

# base custom boxplot
fig, ax = plt.subplots(1, 4, figsize=(16, 7))
ax = ax.ravel()
phen = ['e', 'he', 'hm', 'm']
for p in phen:
    count = 0
    for g in ['u200', 'ZEB', 'LIN28', 'let7']:
        percs = np.percentile(dfs[0].loc[dfs[0]['phenotype'] == p, g].to_numpy(), [5, 15, 35, 50, 65, 85, 95])
        med = percs[3]
        percs = percs[:3].tolist() + percs[4:].tolist()  # remove the median from the list
        helpers.custom_boxplot(ax[phen.index(p)], percs, count, sc10[count], med)
        count += 1
    helpers.significance_annotate((0, 1.5), (1, 1.5), df_significance['p'].to_numpy()[
        df_significance['params'] == 'base_' + p + '_u200_zeb'], ax[phen.index(p)], h=0.1, h2=0.1, h3=0.02)
    helpers.significance_annotate((2, 1.5), (3, 1.5), df_significance['p'].to_numpy()[
        df_significance['params'] == 'base_' + p + '_lin28_let7'], ax[phen.index(p)], h=0.1, h2=0.1, h3=0.02)
for i in ax:
    i.set_ylim(-2, 2)
    i.set_xticks([])
for i in ax[1:]:
    i.set_yticks([])
ax[0].set_yticks(np.linspace(-2, 2, 5))
plt.subplots_adjust(wspace=0.1)
fig.savefig('figures/plots/base_customBox.png', dpi=500, bbox_inches='tight')
plt.close()

# combined cluster center plot
fig, ax = plt.subplots(figsize=(7, 7))
count = 0
for psf, df in zip(psfnames, dfs):
    if psf == 'grhl2partial':
        continue
    for ps in ['e', 'he', 'hm', 'm']:
        bool_mask = (df['phenotype'].to_numpy() == ps)
        mx = np.median(df['ZEB'].to_numpy()[bool_mask])
        my = np.median(df['LIN28'].to_numpy()[bool_mask])
        a2 = ax.scatter(mx, my, marker='*', c=sc1[count], s=70)
        c = plt.Circle((mx, my), np.sum(bool_mask) * 1.3 / len(df), fill=False, ls='-', color=sc1[count], linewidth=3)
        ax.add_artist(c)
    count += 1
ax.set_xticks(np.linspace(-2, 2, 5))
ax.set_yticks(np.linspace(-2, 2, 5))
ax.set_ylabel('LIN28')
ax.set_xlabel('ZEB')
fig.savefig('figures/plots/combined_cluster_centers.png', dpi=500, bbox_inches='tight')
plt.close()

# comparison with uncoupled circuit
# normalization using the base1 run data (to be able to compare it with base): usually not recommended
dfbase = pd.read_csv('collated_data/base/oe0/base1_collated.csv')
dfbase2 = pd.read_csv('collated_data/base/oe0/base1_clusData.csv')
dfuncoupled_emt = pd.read_csv('collated_data/uncoupled/emt_collated.csv')
dfuncoupled_emt = helpers.column_normalize_foo(dfuncoupled_emt, dfbase)
dfuncoupled_stem = pd.read_csv('collated_data/uncoupled/stemness_collated.csv')
dfuncoupled_stem = helpers.column_normalize_foo(dfuncoupled_stem, dfbase)
df_list = {'emt': dfuncoupled_emt, 'stem': dfuncoupled_stem}
genes = {'emt': ['u200', 'ZEB'], 'stem': ['LIN28', 'let7']}

# clustering at multiple levels
final_labels = []
for k in ['emt', 'stem']:  # only the plot for emt is used for the report
    main_data = df_list[k].loc[:, genes[k]].to_numpy()
    major_labs_base = [np.array([]), np.array([])]  # to allow indexing directly by c (number of clusters)
    inertia_vec = []  # store KMeans inertia (one of the cluster metrics)
    for c in range(2, 9):
        kmeans_obj = KMeans(c, init='k-means++', n_init=50, random_state=1, algorithm='full', max_iter=500,
                            tol=1e-5)
        kmeans_obj.fit(main_data)
        major_labs_base.append(np.array(kmeans_obj.labels_, copy=True))
        inertia_vec.append(kmeans_obj.inertia_)
        if c == 2:  # save the 2 cluster solution
            final_labels.append(np.array(kmeans_obj.labels_, copy=True))
    fig, ax = plt.subplots(2, 2, figsize=(12, 12))
    ssvec, chvec, dbvec = [], [], []  # silhouette score, calinski-harabasz criteria, davies-bouldin score
    for c in range(2, 9):
        ss, ch, db = helpers.overall_scores_foo(main_data, major_labs_base[c])
        ssvec.append(ss)
        chvec.append(ch)
        dbvec.append(db)
    print(k, ssvec)
    ax[0, 0].plot([x for x in range(2, 9)], ssvec)
    ax[1, 0].plot([x for x in range(2, 9)], chvec)
    ax[0, 1].plot([x for x in range(2, 9)], dbvec)
    ax[1, 1].plot([x for x in range(2, 9)], inertia_vec)
    ax[0, 0].set_xlabel('Cluster Numbers')
    ax[1, 0].set_xlabel('Cluster Numbers')
    ax[0, 1].set_xlabel('Cluster Numbers')
    ax[1, 1].set_xlabel('Cluster Numbers')
    fig.savefig('figures/plots/uncoupled_' + k + '_clusScores.png', dpi=500, bbox_inches='tight')
    plt.close()


fig, ax = plt.subplots(figsize=(7, 7))
for ps in [(0, 0), (0, 1), (1, 0), (1, 1)]:  # all possible pairings of the clusters in the two circuits
    # 14237 entries in the smaller one
    bool_mask = np.logical_and(final_labels[0][:14237] == ps[0], final_labels[1][:14237] == ps[1])
    mx = np.median(dfuncoupled_emt['ZEB'].to_numpy()[:14237][bool_mask])
    my = np.median(dfuncoupled_stem['LIN28'].to_numpy()[:14237][bool_mask])
    a1 = ax.scatter(mx, my, marker='P', c=sc1[-2], s=70)
    c = plt.Circle((mx, my), np.sum(bool_mask) * 1.3 / 14237, fill=False, ls='-', color=sc1[-2], linewidth=3)
    ax.add_artist(c)
for ps in ['e', 'he', 'hm', 'm']:
    bool_mask = (dfbase2['phenotype'].to_numpy() == ps)
    mx = np.median(dfbase2['ZEB'].to_numpy()[bool_mask])
    my = np.median(dfbase2['LIN28'].to_numpy()[bool_mask])
    a2 = ax.scatter(mx, my, marker='*', c=sc1[1], s=70)
    c = plt.Circle((mx, my), np.sum(bool_mask) * 1.3 / len(dfbase), fill=False, ls='-', color=sc1[1], linewidth=3)
    ax.add_artist(c)
ax.set_xlabel('ZEB')
ax.set_ylabel('LIN28')
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)
ax.set_xticks(np.linspace(-2, 2, 5))
ax.set_yticks(np.linspace(-2, 2, 5))
fig.savefig('figures/plots/uncoupled_clusterCenters.png', dpi=500, bbox_inches='tight')
plt.close()

# dendrogram
main_data = dfs[0].loc[:, ['u200', 'ZEB', 'LIN28', 'let7']].to_numpy()
hc = AgglomerativeClustering(affinity='euclidean', linkage='ward', distance_threshold=0, n_clusters=None)
hc = hc.fit(main_data)  # include data for all merges
fig, ax = plt.subplots(figsize=(7, 7))
hierarchy.dendrogram(helpers.custom_dendrogram_foo(hc, helpers.count_foo(hc), 50), truncate_mode='lastp', p=13, ax=ax)
plt.ylabel('ward-distance')
plt.xlabel('observations')
plt.savefig('figures/plots/base_run1_dendrogram.png', dpi=500, bbox_inches='tight')
plt.close()
