import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

helper = __import__('helpers')

psfnames = ['base', 'grhl2full', 'grhl2partial', 'ovol', 'nrf2']  # grhl2partial: GRHL2-KD
oe = ['de10', 'oe0', 'oe10']

df_cluster = pd.read_csv('figures/clustering/cluster_indices.csv')
df_tables_means = pd.read_csv('figures/tables_numbers/base_dists_means.csv')
df_tables_sds = pd.read_csv('figures/tables_numbers/base_dists_sds.csv')
df_significance = pd.read_csv('figures/tables_numbers/stat_significance.csv')

# color palettes
sc1 = ['#7a5195', '#003f5c', '#ef5675', '#ffa600']  # not using all of these palettes. Kept here for the kind reader :D
sc2 = ['#8a27a3', '#27a38a', '#a38a27']
sc3 = ['#626520', '#3d91e0', '#952b60', '#f1b620']
sc4 = ['#a92420', '#d1913e', '#0e3160']
sc5 = ['#afd275', '#7e685a', '#e7717d']
sc6 = ['#123c69', '#ac3b61', '#edc7b7', '#bab2b5']
sc8 = ['#8222dd', '#ff0064', '#f38600', '#7ddd22']
sc9 = ['#0f1f38', '#005f84', '#00a383', '#57a60d']
sc10 = ['#011627', '#ff9f1c', '#e71d36', '#2ec486']
gr1 = ['#1e3d58', '#057dcd', '#43b0f1']
gr2 = ['#b1d8b7', '#76b947', '#2f5233']
gr3 = ['#d3b1c2', '#c197d2', '#613659']

# change sizes. The axes labels are edited after creating the graph, hence the size kept low here
plt.rc('axes', labelsize=15)
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
plt.rc('legend', fontsize=18)


def staggered_barplot(bars, ebs, ax, labels, ticklabels, colors):  # for custom barplot with multiple bars per location
    n = len(bars)
    assert n == len(ebs)
    w = 0.8 / n  # almost arbitrary value of 0.8
    start = [i - n * w / 2 + w / 2 for i in range(1, len(bars[0]) + 1)]  # where the first kind of barplot starts
    # this starting points are progressively shifted as we plot the remaining kinds of barplots. Here, n=3 almost always
    for i in range(n):
        ax.bar(start, bars[i], color=colors[i], width=w, yerr=ebs[i], ecolor='black', label=labels[i], capsize=6 - n,
               error_kw={'elinewidth': 2, 'capthick': 2})
        start = [i + w for i in start]
    ax.set_yticks([round(x, 2) for x in np.linspace(0, np.max(np.array(bars).flatten()), 5)])
    ax.set_xticks([i for i in range(1, len(bars[0]) + 1)])
    ax.set_xticklabels(ticklabels)


# individual circuit plot for all the indices combined (the ones for the PSF circuits not shown in the paper)
for psf in psfnames:
    fig, ax = plt.subplots(2, 2, figsize=(16, 16))
    ax = ax.ravel()
    ax[0].plot([i for i in range(2, 9)], df_cluster[psf + 'oe0sils_means'], color=gr1[1], label='oe0', linewidth=2)
    ax[0].errorbar([i for i in range(2, 9)], df_cluster[psf + 'oe0sils_means'], yerr=df_cluster[psf + 'oe0sils_sd'],
                   color=gr1[1], fmt='.', ecolor=gr1[1], capsize=4, capthick=2, elinewidth=2)
    ax[0].set_xticks([i for i in range(2, 9)])
    ax[0].set_xlabel('clusters')
    ax[1].plot([i for i in range(2, 9)], df_cluster[psf + 'oe0ch_means'], color=gr1[1], label='oe0', linewidth=2)
    ax[1].errorbar([i for i in range(2, 9)], df_cluster[psf + 'oe0ch_means'], yerr=df_cluster[psf + 'oe0ch_sd'],
                   color=gr1[1], fmt='.', ecolor=gr1[1], capsize=4, capthick=2, elinewidth=2)
    ax[1].set_xticks([i for i in range(2, 9)])
    ax[1].set_xlabel('clusters')
    ax[2].plot([i for i in range(2, 9)], df_cluster[psf + 'oe0db_means'], color=gr1[1], label='oe0', linewidth=2)
    ax[2].errorbar([i for i in range(2, 9)], df_cluster[psf + 'oe0db_means'], yerr=df_cluster[psf + 'oe0db_sd'],
                   color=gr1[1], fmt='.', ecolor=gr1[1], capsize=4, capthick=2, elinewidth=2)
    ax[2].set_xticks([i for i in range(2, 9)])
    ax[2].set_xlabel('clusters')
    ax[3].plot([i for i in range(2, 9)], df_cluster[psf + 'oe0inert_means'], color=gr1[1], label='oe0', linewidth=2)
    ax[3].errorbar([i for i in range(2, 9)], df_cluster[psf + 'oe0inert_means'], yerr=df_cluster[psf + 'oe0inert_sd'],
                   color=gr1[1], fmt='.', ecolor=gr1[1], capsize=4, capthick=2, elinewidth=2)
    ax[3].set_xticks([i for i in range(2, 9)])
    ax[3].set_xlabel('clusters')

    if not (psf == 'base'):  # then plot the graphs for de10 and oe10 (not relevant for base)
        ax[0].plot([i for i in range(2, 9)], df_cluster[psf + 'de10sils_means'], color=gr1[2], label='de10',
                   linewidth=2)
        ax[0].errorbar([i for i in range(2, 9)], df_cluster[psf + 'de10sils_means'],
                       yerr=df_cluster[psf + 'de10sils_sd'],
                       color=gr1[2], fmt='.', ecolor=gr1[2], capthick=2, elinewidth=2, capsize=4)
        ax[0].plot([i for i in range(2, 9)], df_cluster[psf + 'oe10sils_means'], color=gr1[0], label='oe10',
                   linewidth=2)
        ax[0].errorbar([i for i in range(2, 9)], df_cluster[psf + 'oe10sils_means'],
                       yerr=df_cluster[psf + 'oe10sils_sd'],
                       color=gr1[0], fmt='.', ecolor=gr1[0], capsize=4, capthick=2, elinewidth=2)
        ax[1].plot([i for i in range(2, 9)], df_cluster[psf + 'de10ch_means'], color=gr1[2], label='de10',
                   linewidth=2)
        ax[1].errorbar([i for i in range(2, 9)], df_cluster[psf + 'de10ch_means'],
                       yerr=df_cluster[psf + 'de10ch_sd'],
                       color=gr1[2], fmt='.', ecolor=gr1[2], capthick=2, elinewidth=2, capsize=4)
        ax[1].plot([i for i in range(2, 9)], df_cluster[psf + 'oe10ch_means'], color=gr1[0], label='oe10',
                   linewidth=2)
        ax[1].errorbar([i for i in range(2, 9)], df_cluster[psf + 'oe10ch_means'],
                       yerr=df_cluster[psf + 'oe10ch_sd'],
                       color=gr1[0], fmt='.', ecolor=gr1[0], capsize=4, capthick=2, elinewidth=2)
        ax[2].plot([i for i in range(2, 9)], df_cluster[psf + 'de10db_means'], color=gr1[2], label='de10',
                   linewidth=2)
        ax[2].errorbar([i for i in range(2, 9)], df_cluster[psf + 'de10db_means'],
                       yerr=df_cluster[psf + 'de10db_sd'],
                       color=gr1[2], fmt='.', ecolor=gr1[2], capthick=2, elinewidth=2, capsize=4)
        ax[2].plot([i for i in range(2, 9)], df_cluster[psf + 'oe10db_means'], color=gr1[0], label='oe10',
                   linewidth=2)
        ax[2].errorbar([i for i in range(2, 9)], df_cluster[psf + 'oe10db_means'],
                       yerr=df_cluster[psf + 'oe10db_sd'],
                       color=gr1[0], fmt='.', ecolor=gr1[0], capsize=4, capthick=2, elinewidth=2)
        ax[3].plot([i for i in range(2, 9)], df_cluster[psf + 'de10inert_means'], color=gr1[2], label='de10',
                   linewidth=2)
        ax[3].errorbar([i for i in range(2, 9)], df_cluster[psf + 'de10inert_means'],
                       yerr=df_cluster[psf + 'de10inert_sd'],
                       color=gr1[2], fmt='.', ecolor=gr1[2], capthick=2, elinewidth=2, capsize=4)
        ax[3].plot([i for i in range(2, 9)], df_cluster[psf + 'oe10inert_means'], color=gr1[0], label='oe10',
                   linewidth=2)
        ax[3].errorbar([i for i in range(2, 9)], df_cluster[psf + 'oe10inert_means'],
                       yerr=df_cluster[psf + 'oe10inert_sd'],
                       color=gr1[0], fmt='.', ecolor=gr1[0], capsize=4, capthick=2, elinewidth=2)
    plt.subplots_adjust()
    fig.savefig('figures/plots/' + psf + '_individual_all_indices.png', dpi=500, bbox_inches='tight')
    plt.close()

# combined circuit plots for silhouettes (only reference circuits)
fig, ax = plt.subplots(figsize=(7, 7))
ax.plot([i for i in range(2, 9)], df_cluster['baseoe0sils_means'], color=sc1[0], label='base', linewidth=2)
ax.errorbar([i for i in range(2, 9)], df_cluster['baseoe0sils_means'], yerr=df_cluster['baseoe0sils_sd'], color=sc1[0],
            fmt='.', ecolor=sc1[0], capsize=4, capthick=2, elinewidth=2)
ax.plot([i for i in range(2, 9)], df_cluster['grhl2fulloe0sils_means'], color=sc1[1], label='GRHL2', linewidth=2)
ax.errorbar([i for i in range(2, 9)], df_cluster['grhl2fulloe0sils_means'], yerr=df_cluster['grhl2fulloe0sils_sd'],
            color=sc1[1], fmt='.', ecolor=sc1[1], capsize=4, capthick=2, elinewidth=2)
ax.plot([i for i in range(2, 9)], df_cluster['ovoloe0sils_means'], color=sc1[2], label='OVOL', linewidth=2)
ax.errorbar([i for i in range(2, 9)], df_cluster['ovoloe0sils_means'], yerr=df_cluster['ovoloe0sils_sd'], color=sc1[2],
            fmt='.', ecolor=sc1[2], capsize=4, capthick=2, elinewidth=2)
ax.plot([i for i in range(2, 9)], df_cluster['nrf2oe0sils_means'], color=sc1[3], label='NRF2', linewidth=2)
ax.errorbar([i for i in range(2, 9)], df_cluster['nrf2oe0sils_means'], yerr=df_cluster['nrf2oe0sils_sd'], color=sc1[3],
            fmt='.', ecolor=sc1[3], capsize=4, capthick=2, elinewidth=2)
ax.set_ylabel('avg silhouettes')
ax.set_xlabel('clusters')
ax.set_xticks([i for i in range(2, 9)])
ax.set_yticks([round(x, 2) for x in np.linspace(0.25, 0.47, 5)])
fig.savefig('figures/plots/combined_silhouettes_ebars.png', dpi=500, bbox_inches='tight')
plt.close()

# number of states across all circuits
fig, ax = plt.subplots(figsize=(7, 7))
ax.bar([i for i in range(6)], df_tables_means['baseoe0'].to_numpy()[:6], color=sc1[0], label='base',
       yerr=df_tables_sds['baseoe0'].to_numpy()[:6], ecolor=sc1[0], capsize=4, alpha=0.7,
       error_kw={'elinewidth': 2, 'capthick': 2})
ax.plot([i for i in range(6)], df_tables_means['grhl2fulloe0'].to_numpy()[:6], color=sc1[1], label='GRHL2', linewidth=2)
ax.errorbar([i for i in range(6)], df_tables_means['grhl2fulloe0'].to_numpy()[:6],
            yerr=df_tables_sds['grhl2fulloe0'].to_numpy()[:6], color=sc1[1], fmt='.', ecolor=sc1[1], capsize=4,
            capthick=2, elinewidth=2)
ax.plot([i for i in range(6)], df_tables_means['ovoloe0'].to_numpy()[:6], color=sc1[2], label='OVOL', linewidth=2)
ax.errorbar([i for i in range(6)], df_tables_means['ovoloe0'].to_numpy()[:6],
            yerr=df_tables_sds['ovoloe0'].to_numpy()[:6], color=sc1[2], fmt='.', ecolor=sc1[2], capsize=4, capthick=2,
            elinewidth=2)
ax.plot([i for i in range(6)], df_tables_means['nrf2oe0'].to_numpy()[:6], color=sc1[3], label='NRF2', linewidth=2)
ax.errorbar([i for i in range(6)], df_tables_means['nrf2oe0'].to_numpy()[:6],
            yerr=df_tables_sds['nrf2oe0'].to_numpy()[:6], color=sc1[3], fmt='.', ecolor=sc1[3], capsize=4, capthick=2,
            elinewidth=2)
ax.set_xlabel('number of states')
ax.set_ylabel('proportion')
ax.set_xticks([i for i in range(6)])
ax.set_xticklabels(['1', '2', '3', '4', '5', '>5'])
ax.set_yticks([round(x, 2) for x in np.linspace(0, 0.4, 5)])
fig.savefig('figures/plots/num_states_all.png', dpi=500, bbox_inches='tight')
plt.close()

# monostable states separately for each circuit
for psf in psfnames:
    fig, ax = plt.subplots(figsize=(7, 7))
    if psf == 'base':  # staggered barplots not needed
        ax.bar([i for i in range(4)], df_tables_means[psf + 'oe0'].to_numpy()[10:14],
               yerr=df_tables_sds[psf + 'oe0'].to_numpy()[10:14], color=gr1[1], ecolor='black', capsize=4,
               error_kw={'elinewidth': 2, 'capthick': 2})
        ax.set_xticks([i for i in range(4)])
        ax.set_xticklabels(['e', 'he', 'hm', 'm'])
        ax.set_yticks([round(x, 2) for x in np.linspace(0, np.max(df_tables_means[psf + 'oe0'].to_numpy()[10:14]), 4)])
        helper.significance_annotate((0, 0.43), (3, 0.43), df_significance['p'].to_numpy()[
            df_significance['params'] == 'base_monostable_e_m'], ax, h=0.015, h2=0.05, h3=0.01)
        helper.significance_annotate((1, 0.11), (2, 0.11), df_significance['p'].to_numpy()[
            df_significance['params'] == 'base_monostable_he_hm'], ax, h=0.015, h2=0.05, h3=0.01)
    else:
        staggered_barplot([df_tables_means[psf + x].to_list()[10:14] for x in oe],
                          [df_tables_sds[psf + x].to_list()[10:14] for x in oe], ax, oe, ['e', 'he', 'hm', 'm'], gr1)
    ax.set_xlabel('monostable states')
    ax.set_ylabel('proportion')
    fig.savefig('figures/plots/' + psf + '_monostable.png', dpi=500, bbox_inches='tight')
    plt.close()

# bistable
order_states = ['e/he', 'e/hm', 'e/m', 'he/m', 'hm/m', 'he/hm']
ylabels = ['{' + x.split('/')[0] + ' ,' + x.split('/')[1] + '}' for x in order_states]
for psf in psfnames:
    fig, ax = plt.subplots(figsize=(7, 7))
    if psf == 'base':
        ax.bar([i for i in range(6)], df_tables_means[psf + 'oe0'].to_numpy()[14:20],
               yerr=df_tables_sds[psf + 'oe0'].to_numpy()[14:20], color=gr1[1], ecolor='black', capsize=4,
               error_kw={'elinewidth': 2, 'capthick': 2})
        ax.set_xticks([i for i in range(6)])
        ax.set_xticklabels(ylabels, rotation=45)
        ax.set_yticks([round(x, 2) for x in np.linspace(0, np.max(df_tables_means[psf + 'oe0'].to_numpy()[14:20]), 4)])
    else:
        staggered_barplot([df_tables_means[psf + x].to_list()[14:20] for x in oe],
                          [df_tables_sds[psf + x].to_list()[14:20] for x in oe], ax, oe, ylabels, gr1)
        plt.xticks(rotation=45)
    ax.set_xlabel('bistable states')
    ax.set_ylabel('proportion')
    fig.savefig('figures/plots/' + psf + '_bistable.png', dpi=500, bbox_inches='tight')
    plt.close()

# tristable
order_states = ['e/he/m', 'e/hm/m', 'e/he/hm', 'he/hm/m']
ylabels = ['{' + x.split('/')[0] + ', ' + x.split('/')[1] + ', ' + x.split('/')[2] + '}' for x in order_states]
for psf in psfnames:
    fig, ax = plt.subplots(figsize=(7, 7))
    if psf == 'base':
        ax.bar([i for i in range(4)], df_tables_means[psf + 'oe0'].to_numpy()[20:24],
               yerr=df_tables_sds[psf + 'oe0'].to_numpy()[20:24], color=gr1[1], ecolor='black', capsize=4,
               error_kw={'elinewidth': 2, 'capthick': 2})
        ax.set_xticks([i for i in range(4)])
        ax.set_xticklabels(ylabels, rotation=30)
        ax.set_yticks([round(x, 2) for x in np.linspace(0, np.max(df_tables_means[psf + 'oe0'].to_numpy()[20:24]), 4)])
    else:
        staggered_barplot([df_tables_means[psf + x].to_list()[20:24] for x in oe],
                          [df_tables_sds[psf + x].to_list()[20:24] for x in oe], ax, oe, ylabels, gr1)
        plt.xticks(rotation=30)
    ax.set_xlabel('tristable states')
    ax.set_ylabel('proportion')
    fig.savefig('figures/plots/' + psf + '_tristable.png', dpi=500, bbox_inches='tight')
    plt.close()


# combined monostable enrichment plot
def plot_fill(ax, position, width, arr, colors, ers):  # custom function to show all circuits' enrichment together
    ax.fill_between([position - width / 2, position + width / 2], [0, 0], [arr[0], arr[0]], color=colors[0])
    ax.fill_between([position - width / 2, position + width / 2], [arr[0], arr[0]], [arr[0] + arr[1], arr[0] + arr[1]],
                    color=colors[1])
    ax.fill_between([position - width / 2, position + width / 2], [arr[0] + arr[1], arr[0] + arr[1]],
                    [arr[0] + arr[1] + arr[2], arr[0] + arr[1] + arr[2]], color=colors[2])
    ax.fill_between([position - width / 2, position + width / 2], [arr[0] + arr[1] + arr[2], arr[0] + arr[1] + arr[2]],
                    [arr[0] + arr[1] + arr[2] + arr[3], arr[0] + arr[1] + arr[2] + arr[3]], color=colors[3])
    ax.errorbar([position] * 4, [arr[0], arr[0] + arr[1], arr[0] + arr[1] + arr[2], arr[0] + arr[1] + arr[2] + arr[3]],
                yerr=ers, color='black', capsize=4, capthick=2, elinecolor='black', fmt='none', elinewidth=2)


fig, ax = plt.subplots(figsize=(7, 7))
count = 0
for psf in psfnames:
    if psf in ['base', 'grhl2partial']:  # skip base and GRHL2-KD
        continue
    plot_fill(ax, count + 0.1, 0.2, df_tables_means[psf + 'de10'].to_list()[10:14], sc10,
              df_tables_sds[psf + 'de10'].to_list()[10:14])
    plot_fill(ax, count + 0.35, 0.2, df_tables_means[psf + 'oe0'].to_list()[10:14], sc10,
              df_tables_sds[psf + 'oe0'].to_list()[10:14])
    plot_fill(ax, count + 0.6, 0.2, df_tables_means[psf + 'oe10'].to_list()[10:14], sc10,
              df_tables_sds[psf + 'oe10'].to_list()[10:14])
    count += 1
ticks = [[i + 0.1, i + 0.35, i + 0.6] for i in range(3)]
ticklabels = [['d', 'r', 'o'] for i in range(3)]
ax.set_xticks(np.array(ticks).flatten())
ax.set_xticklabels(np.array(ticklabels).flatten())
ax.set_yticks(np.linspace(0, 1, 5))
fig.savefig('figures/plots/combined_monostable_enrichment.png', dpi=500, bbox_inches='tight')
plt.close()
