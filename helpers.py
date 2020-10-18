import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score, silhouette_samples, calinski_harabasz_score
from sklearn.metrics import davies_bouldin_score


# Helper Functions
def shifted_hill_foo(y, b0, b, n):  # Shifted Hill Function
    # also representable as y + (1 - y)(1/(1 + (b/b0) ** n))
    return (b0 ** n) / ((b0 ** n) + (b ** n)) + y * (b ** n) / ((b0 ** n) + (b ** n))


def column_normalize_foo(dframe, reference_dframe):
    # To normalize the gene expression values (log2-scale) using the mean and sd from the reference_run
    allgenes = ['u200', 'ZEB', 'SNAIL', 'LIN28', 'let7', 'NFkB', 'GRHL2', 'KEAP1', 'NRF2', 'ECad',
                'OVOL']  # all possible genes in all the runs
    newdframe = pd.DataFrame()
    for i in dframe.columns:
        x = dframe[i].to_numpy()
        if i in allgenes:
            if i in reference_dframe.columns:
                m = np.mean(reference_dframe[i])
                s = np.std(reference_dframe[i], ddof=1)
                newdframe[i] = (x - m) / s
            else:  # if the gene not present in the reference DataFrame (only for completeness)
                m = np.mean(dframe[i])
                s = np.std(dframe[i], ddof=1)
                newdframe[i] = (x - m) / s
        else:  # if the column is not a gene-column
            newdframe[i] = x
    return newdframe


def count_foo(model):
    # To count children at different levels of a HierarchicalClustering tree
    # modified from https://scikit-learn.org/stable/auto_examples/cluster/plot_agglomerative_dendrogram.html
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count
    return counts


def custom_dendrogram_foo(hc, counts, n):
    # To return the linkage_matrix expected by scipy.cluster.hierarchy.dendrogram having data of only the last n merges
    # Explanation:
    # nodes merged at ith iteration have index n + i, with original observations having i < n
    # for z = linkage matrix, z[i,0] and z[i,1] are merged at ith iteration to form a cluster of index n+i, and these
    # nodes have a distance z[i,2] and the number of original observations (i<n) in the new cluster is z[i,3]
    # Also, number of merges = number of observations - 1
    children = hc.children_
    distances = hc.distances_
    original_n = len(counts) + 1
    seen_clusters = dict()  # the cluster numbers we have already observed
    linkage_matrix = np.zeros((n, 4))
    iteration = 0
    leafnodes = -1  # begin at -1 so that the first assigned value is 0 (leafnodes + 1)
    for c1, c2 in children[-n:, :]:
        if c1 not in seen_clusters:
            seen_clusters[c1] = leafnodes + 1
            leafnodes += 1
        if c2 not in seen_clusters:
            seen_clusters[c2] = leafnodes + 1
            leafnodes += 1
        c1 = seen_clusters[c1]
        c2 = seen_clusters[c2]
        seen_clusters[original_n + len(counts) - n + iteration] = n + 1 + iteration
        linkage_matrix[iteration, :] = [c1, c2, distances[-n + iteration], counts[-n + iteration]]
        iteration += 1
    return linkage_matrix


def hierarchical_cluster_labels(hc, n):
    # To cut the HierarchicalClustering tree to obtain n clusters and return the labels
    # starts joining clusters from bottom until the requisite number is reached
    labs = np.array([x for x in range(len(hc.labels_))])
    cluster_groups = dict([(x, [x]) for x in range(len(hc.labels_))])
    curr_clusters = len(hc.labels_)
    n_observations = len(hc.labels_)
    i = 0
    while curr_clusters > n:
        c1, c2 = hc.children_[i, :]
        gr = cluster_groups.pop(c1) + cluster_groups.pop(c2)
        labs[np.array(gr)] = labs[gr[0]]  # choose any one common label to apply to all
        cluster_groups[n_observations + i] = gr
        curr_clusters -= 1
        i += 1
    uq = np.unique(labs).tolist()
    labs = [uq.index(x) for x in labs]  # assign cluster numbers to observations based on the unique labels found
    return labs


def overall_scores_foo(data, labs):
    # To compute the cluster-quality metrics for the given labels and data using the euclidean distance norm
    ss = silhouette_score(data, metric='euclidean', labels=labs)
    ch = calinski_harabasz_score(data, labels=labs)
    db = davies_bouldin_score(data, labels=labs)
    return ss, ch, db


def silhouette_plot(data, labs, ax):  # not used. Kept for the benefit of a possible (though improbable) reader :D
    # To plot a Silhouette plot (individual sample silhouette widths arranged in descending fashion)
    ss = silhouette_samples(data, labels=labs, metric='euclidean')
    ss.sort()
    ax.fill_between([i for i in range(len(ss))], ss, [0 for i in range(len(ss))])
    ax.hlines(np.mean(ss), 0, len(ss), linestyles='dashed')  # line for the average silhouette width
    ax.set_xlabel('Sorted Observations')
    ax.set_ylabel('Silhouette widths for each observation')


def custom_boxplot(ax, arr, position, col, median):
    # Prettier custom boxplot-equivalent with colors
    # arr contains the requisite percentiles, position is the center of the box (0.5 is the width of the box)
    ax.fill_between([position - 0.25, position + 0.25], [arr[2], arr[2]], [arr[3], arr[3]], color=col)
    ax.fill_between([position - 0.25, position + 0.25], [arr[1], arr[1]], [arr[2], arr[2]], color=col, alpha=0.6)
    ax.fill_between([position - 0.25, position + 0.25], [arr[3], arr[3]], [arr[4], arr[4]], color=col, alpha=0.6)
    ax.fill_between([position - 0.25, position + 0.25], [arr[0], arr[0]], [arr[1], arr[1]], color=col, alpha=0.3)
    ax.fill_between([position - 0.25, position + 0.25], [arr[4], arr[4]], [arr[5], arr[5]], color=col, alpha=0.3)
    ax.plot([position - 0.4, position + 0.4], [median, median], color=col, linewidth=2)


def bond_strength(y, mu, bigG, k, interaction_type='Act'):
    if interaction_type == 'Act':
        return y * bigG / mu / k  # bond strength
    elif interaction_type == 'Inh':
        return bigG / mu / k / y  # inverted lambda for inhibitory relations
    else:
        assert False, 'Incorrect Type.'


def asymmetry_bond_strength(b1, b2):
    # returns asymmetry parameter for two bonds (used for double negative feedback loops)
    log_ratio = np.log2(b1 / b2)
    return log_ratio


def significance_annotate(xy1, xy2, p, ax, h=0.5, h2=0.3, h3=0.1):
    # allows adding the significance bars to denote p-value in the plots
    if p > 0.01:
        s = '+'
    elif p > 0.001:
        s = '*'
    elif p > 0.0001:
        s = '**'
    else:
        s = '***'
    m = max(xy1[1], xy2[1])  # located at the higher of the two y-values
    ax.plot([xy1[0], xy1[0]], [m + h2, m + h2 + h], color='black')
    ax.plot([xy1[0], xy2[0]], [m + h2 + h, m + h2 + h], color='black')
    ax.plot([xy2[0], xy2[0]], [m + h2, m + h2 + h], color='black')
    ax.text((xy1[0] + xy2[0]) / 2, m + h2 + h + h3, s, horizontalalignment='center', fontsize=14)


def significance_annotate_inverted(xy1, xy2, p, ax, h=0.5, h2=0.3, h3=0.1):
    # allows adding the significance bars to denote p-value in the plots
    if p > 0.01:
        s = '+'
    elif p > 0.001:
        s = '*'
    elif p > 0.0001:
        s = '**'
    else:
        s = '***'
    m = max(xy1[1], xy2[1])  # located at the higher of the two y-values
    ax.plot([xy1[0], xy1[0]], [m - h2, m - h2 - h], color='black')
    ax.plot([xy1[0], xy2[0]], [m - h2 - h, m - h2 - h], color='black')
    ax.plot([xy2[0], xy2[0]], [m - h2, m - h2 - h], color='black')
    ax.text((xy1[0] + xy2[0]) / 2, m - h2 - h - h3, s, horizontalalignment='center', verticalalignment='top',
            fontsize=14)
