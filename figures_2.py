import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

helper = __import__('helpers')

psfnames = ['base', 'grhl2full', 'grhl2partial', 'ovol', 'nrf2']
oe = ['de10', 'oe0', 'oe10']

df_tables_means = pd.read_csv('figures/tables_numbers/stemness_means.csv')
df_tables_sds = pd.read_csv('figures/tables_numbers/stemness_sds.csv')
df_significance = pd.read_csv('figures/tables_numbers/stat_significance.csv')

sc1 = ['#7a5195', '#003f5c', '#ef5675', '#ffa600']
sc10 = ['#011627', '#ff9f1c', '#e71d36', '#2ec486']
gr1 = ['#1e3d58', '#057dcd', '#43b0f1']

plt.rc('axes', labelsize=15)
plt.rc('xtick', labelsize=20)
plt.rc('ytick', labelsize=20)
plt.rc('legend', fontsize=18)


def staggered_barplot(bars, ebs, ax, labels, ticklabels, colors):
    n = len(bars)
    assert n == len(ebs)
    w = 0.8 / n
    start = [i - n * w / 2 + w / 2 for i in range(1, len(bars[0]) + 1)]
    for i in range(n):
        ax.bar(start, bars[i], color=colors[i], width=w, yerr=ebs[i], ecolor='black', label=labels[i], capsize=6 - n,
               error_kw={'elinewidth': 2, 'capthick': 2})
        start = [i + w for i in start]
    ax.set_yticks([round(x, 2) for x in np.linspace(0, np.max(np.array(bars).flatten()), 5)])
    ax.set_xticks([i for i in range(1, len(bars[0]) + 1)])
    ax.set_xticklabels(ticklabels)


# p1 plots for all circuits
for psf in psfnames:
    fig, ax = plt.subplots(figsize=(7, 7))
    if psf == 'base':
        ax.bar([i for i in range(4)], df_tables_means[psf + 'oe0'].to_numpy()[0:4],
               yerr=df_tables_sds[psf + 'oe0'].to_numpy()[0:4], color=gr1[1], ecolor='black', capsize=4,
               error_kw={'elinewidth': 2, 'capthick': 2})
        ax.set_xticks([i for i in range(4)])
        ax.set_xticklabels(['e', 'he', 'hm', 'm'])
        ax.set_yticks([round(x, 2) for x in np.linspace(0, np.max(df_tables_means[psf + 'oe0'].to_numpy()[0:4]), 5)])
        helper.significance_annotate((0, 0.34), (1, 0.34),
                                     df_significance['p'].to_numpy()[df_significance['params'] == 'base_p1_e_he'],
                                     ax, h=0.015, h2=0.05, h3=0.01)
        helper.significance_annotate((2, 0.47), (3, 0.47),
                                     df_significance['p'].to_numpy()[df_significance['params'] == 'base_p1_hm_m'],
                                     ax, h=0.015, h2=0.05, h3=0.01)
        helper.significance_annotate((0, 0.54), (3, 0.54),
                                     df_significance['p'].to_numpy()[df_significance['params'] == 'base_p1_e_m'],
                                     ax, h=0.015, h2=0.05, h3=0.01)
        helper.significance_annotate((1, 0.5), (2, 0.5),
                                     df_significance['p'].to_numpy()[df_significance['params'] == 'base_p1_he_hm'],
                                     ax, h=0.015, h2=0.05, h3=0.01)
    else:
        staggered_barplot([df_tables_means[psf + x].to_list()[0:4] for x in oe],
                          [df_tables_sds[psf + x].to_list()[0:4] for x in oe], ax, oe, ['e', 'he', 'hm', 'm'], gr1)
    ax.set_xlabel('phenotypes')
    ax.set_ylabel('p1')
    fig.savefig('figures/plots/' + psf + '_stemness_p1.png', dpi=500, bbox_inches='tight')
    plt.close()

# p2 plots
for psf in psfnames:
    fig, ax = plt.subplots(figsize=(7, 7))
    if psf == 'base':
        ax.bar([i for i in range(4)], df_tables_means[psf + 'oe0'].to_numpy()[4:8],
               yerr=df_tables_sds[psf + 'oe0'].to_numpy()[4:8], color=gr1[1], ecolor='black', capsize=4,
               error_kw={'elinewidth': 2, 'capthick': 2})
        ax.set_xticks([i for i in range(4)])
        ax.set_xticklabels(['e', 'he', 'hm', 'm'])
        ax.set_yticks([round(x, 2) for x in np.linspace(0, np.max(df_tables_means[psf + 'oe0'].to_numpy()[4:8]), 5)])
        helper.significance_annotate((0, 0.26), (1, 0.26),
                                     df_significance['p'].to_numpy()[df_significance['params'] == 'base_p2_e_he'],
                                     ax, h=0.015, h2=0.05, h3=0.01)
        helper.significance_annotate((2, 0.35), (3, 0.35),
                                     df_significance['p'].to_numpy()[df_significance['params'] == 'base_p2_hm_m'],
                                     ax, h=0.015, h2=0.05, h3=0.01)
        helper.significance_annotate((0, 0.41), (3, 0.41),
                                     df_significance['p'].to_numpy()[df_significance['params'] == 'base_p2_e_m'],
                                     ax, h=0.015, h2=0.05, h3=0.01)
        helper.significance_annotate((1, 0.38), (2, 0.38),
                                     df_significance['p'].to_numpy()[df_significance['params'] == 'base_p2_he_hm'],
                                     ax, h=0.015, h2=0.05, h3=0.01)

    else:
        staggered_barplot([df_tables_means[psf + x].to_list()[4:8] for x in oe],
                          [df_tables_sds[psf + x].to_list()[4:8] for x in oe], ax, oe, ['e', 'he', 'hm', 'm'], gr1)
    ax.set_xlabel('phenotypes')
    ax.set_ylabel('p2')
    fig.savefig('figures/plots/' + psf + '_stemness_p2.png', dpi=500, bbox_inches='tight')
    plt.close()

# p1/p2 comparison across circuits (only oe0) hybrid vs non-hybrids
fig, ax = plt.subplots(figsize=(5, 7))
vals = [df_tables_means[x + 'oe0'].to_list()[8] for x in psfnames if not (x == 'grhl2partial')]
vals2 = [df_tables_sds[x + 'oe0'].to_list()[8] for x in psfnames if not (x == 'grhl2partial')]
vals_2 = [df_tables_means[x + 'oe0'].to_list()[9] for x in psfnames if not (x == 'grhl2partial')]
vals2_2 = [df_tables_sds[x + 'oe0'].to_list()[9] for x in psfnames if not (x == 'grhl2partial')]
staggered_barplot([vals, vals_2], [vals2, vals2_2], ax, ['hybrid', 'non-hybrid'], ['base', 'GRHL2', 'OVOL', 'NRF2'],
                  [sc1[1], sc1[-1]])
fig.savefig('figures/plots/combined_p1_hyb_nonhyb.png', dpi=500, bbox_inches='tight')
plt.close()

fig, ax = plt.subplots(figsize=(5, 7))
vals = [df_tables_means[x + 'oe0'].to_list()[10] for x in psfnames if not (x == 'grhl2partial')]
vals2 = [df_tables_sds[x + 'oe0'].to_list()[10] for x in psfnames if not (x == 'grhl2partial')]
vals_2 = [df_tables_means[x + 'oe0'].to_list()[11] for x in psfnames if not (x == 'grhl2partial')]
vals2_2 = [df_tables_sds[x + 'oe0'].to_list()[11] for x in psfnames if not (x == 'grhl2partial')]
ax.fill_between([0, 0.5], [0, 0], [vals[0], vals[0]], color=sc1[1])
ax.fill_between([0.75, 1.25], [0, 0], [vals[1], vals[1]], color=sc1[1])
ax.fill_between([1.5, 2], [0, 0], [vals[2], vals[2]], color=sc1[1])
ax.fill_between([2.25, 2.75], [0, 0], [vals[3], vals[3]], color=sc1[1])
ax.fill_between([0, 0.5], [vals[0], vals[0]], [vals[0] + vals_2[0], vals[0] + vals_2[0]], color=sc1[-1])
ax.fill_between([0.75, 1.25], [vals[1], vals[1]], [vals[1] + vals_2[1], vals[1] + vals_2[1]], color=sc1[-1])
ax.fill_between([1.5, 2], [vals[2], vals[2]], [vals[2] + vals_2[2], vals[2] + vals_2[2]], color=sc1[-1])
ax.fill_between([2.25, 2.75], [vals[3], vals[3]], [vals[3] + vals_2[3], vals[3] + vals_2[3]], color=sc1[-1])
ax.errorbar([0.25, 1, 1.75, 2.5], vals, yerr=vals2, elinewidth=2, capsize=5, capthick=2, fmt='none', ecolor='black')
ax.errorbar([0.25, 1, 1.75, 2.5], [vals[i] + vals_2[i] for i in range(4)], yerr=vals2_2, elinewidth=2, capsize=5,
            capthick=2, fmt='none', ecolor='black')
ax.set_xticks([0.25, 1, 1.75, 2.5])
ax.set_xticklabels(['base', 'GRHL2', 'OVOL', 'NRF2'])
fig.savefig('figures/plots/combined_p2_hyb_nonhyb.png', dpi=500, bbox_inches='tight')
plt.close()

# sliding stemness window plot for base
df_list = [pd.read_csv('collated_data/base/oe0/base' + str(x) + '_clusData.csv') for x in range(1, 6)]
df_list_2 = [pd.read_csv('collated_data/base/oe0/base' + str(x) + '_collated.csv') for x in range(1, 6)]
# import the unnormalized as well as the normalized base DataFrames
norm_values = [(np.mean(x['LIN28']), np.std(x['LIN28'], ddof=1)) for x in df_list_2]
# get the Z-normalization constants for all the runs in the base circuit

unnormalized_lin28_range = (-11.81, 12.08)  # run tables_and_numbers_2.py to print this value
ulr = unnormalized_lin28_range  # abbreviation for easy typing (purpose defeated by this long comment)
middle_percentages = np.linspace(10, 90, 20)
mn = np.mean(ulr)
rng = ulr[1] - ulr[0]
candidate_stem_windows = [(mn - x / 2 / 100 * rng, mn + x / 2 / 100 * rng) for x in middle_percentages]
combined_probs = []
for p in ['e', 'he', 'hm', 'm']:
    phenotype_probs = []
    for sw in candidate_stem_windows:
        prob_vec = []
        for i in range(5):
            norm_sw = ((sw[0] - norm_values[i][0]) / norm_values[i][1], (sw[1] - norm_values[i][0]) / norm_values[i][1])
            df = df_list[i]
            hs_boolmask = np.logical_and(df['LIN28'] >= norm_sw[0], df['LIN28'] <= norm_sw[1])
            prob_vec.append(np.sum(np.logical_and(df['phenotype'] == p, hs_boolmask)) / np.sum(hs_boolmask))
        phenotype_probs.append(np.mean(prob_vec))
    combined_probs.append(np.array(phenotype_probs))
fig, ax = plt.subplots(figsize=(7, 7))
ax.fill_between(middle_percentages, [0] * 20, combined_probs[0], label='e', color=sc1[0])
ax.fill_between(middle_percentages, combined_probs[0], combined_probs[0] + combined_probs[1], label='he', color=sc1[1])
ax.fill_between(middle_percentages, combined_probs[0] + combined_probs[1],
                combined_probs[0] + combined_probs[1] + combined_probs[2], label='hm', color=sc1[2])
ax.fill_between(middle_percentages, combined_probs[0] + combined_probs[1] + combined_probs[2],
                combined_probs[0] + combined_probs[1] + combined_probs[2] + combined_probs[3], label='m', color=sc1[3])
xtick_vec = [middle_percentages[i] for i in range(0, 20, 3)]
ax.set_xticks(xtick_vec)
ax.set_xticklabels([str(round(x, 1)) for x in xtick_vec])
ax.set_yticks([round(x, 2) for x in np.linspace(0, 1, 5)])
ax.axvline(30, 0, 1, linestyle='--', color='black')
ax.set_xlabel('middle % of LIN28 range')
ax.set_ylabel('p2 for each phenotype')
ax.legend()
fig.savefig('figures/plots/sliding_stemness.png', dpi=500, bbox_inches='tight')
plt.close()
