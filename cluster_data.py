import numpy as np
import pandas as pd
import os
from sklearn.cluster import KMeans
import time  # To know how much time it took (for time management)

helpers = __import__('helpers')

# To compute the clustering of core genes and plot/tabulate basic clustering information
t_start = time.time()

np.random.seed(12345)

paths = ['GRHL2/full/oe0', 'GRHL2/full/oe10', 'GRHL2/full/de10',
         'GRHL2/partial/oe0', 'GRHL2/partial/oe10', 'GRHL2/partial/de10',
         'NRF2/oe0', 'NRF2/oe10', 'NRF2/de10',
         'OVOL/oe0', 'OVOL/oe10', 'OVOL/de10',
         'base/oe0']

raw_data_paths = ['collated_data/' + x for x in paths]  # where the csv files are stored
save_paths = ['figures/clustering/'] * 13  # where the plots will be stored (kept as a list for uniformity)

df_cluster_indices = pd.DataFrame()  # to store all the cluster indices to later use for plotting and for easy reference

for path_to_folder, save_path in zip(raw_data_paths, save_paths):
    t_per_folder = time.time()
    print('\nStarting ' + path_to_folder)
    files = os.listdir(path_to_folder)
    circuit_names = [f.split('_')[0] for f in files if 'oe0Normalized.csv' in f]  # names of the runs in the folder
    oe_level = path_to_folder.split('/')[-1]  # the over/down-expression level (from among oe0, oe10, de10)
    # assuming no name has '_' in it
    df_list = []  # stores the DataFrames for all the runs/replicates
    df_param_list = []
    for circuit_name in circuit_names:
        df_list.append(pd.read_csv(path_to_folder + '/' + circuit_name + '_oe0Normalized.csv'))
        df_param_list.append(pd.read_csv(path_to_folder + '/' + circuit_name + '_parameters_formatted.csv'))
    psfname = circuit_names[0][:-1]  # from among ovol, grhl2full, grhl2partial, base, nrf2

    tier1_genes = ['u200', 'ZEB', 'LIN28', 'let7']  # coreGenes
    # main clustering and saving the DataFrames along with phenotype calling
    final_labels = []  # to store the final labels of all replicates (corresponding to cluster number = 4)
    main_ss_vec = []  # to store the cluster quality indices for different runs
    main_db_vec = []
    main_ch_vec = []
    main_inertia_vec = []
    for run_index in range(5):
        main_data = df_list[run_index].loc[:, tier1_genes].to_numpy()
        major_labs_base = [np.array([]), np.array([])]  # to allow indexing directly by c (number of clusters)
        inertia_vec = []  # store KMeans inertia (one of the cluster metrics)
        for c in range(2, 9):
            kmeans_obj = KMeans(c, init='k-means++', n_init=50, random_state=1, algorithm='full', max_iter=500,
                                tol=1e-5)
            kmeans_obj.fit(main_data)
            major_labs_base.append(np.array(kmeans_obj.labels_, copy=True))
            inertia_vec.append(kmeans_obj.inertia_)
            if c == 4:  # save the 4 cluster solution
                final_labels.append(np.array(kmeans_obj.labels_, copy=True))
        ssvec, chvec, dbvec = [], [], []  # Silhouette score, Calinski-Harabasz criteria, Davies-Bouldin score
        for c in range(2, 9):
            ss, ch, db = helpers.overall_scores_foo(main_data, major_labs_base[c])
            ssvec.append(ss)
            chvec.append(ch)
            dbvec.append(db)
        main_ss_vec.append(ssvec)
        main_db_vec.append(dbvec)
        main_ch_vec.append(chvec)
        main_inertia_vec.append(inertia_vec)

    # add the across-replicate means and sds for each of the indices for each of the cluster numbers from 2 - 8
    column_prefix = psfname + oe_level
    df_cluster_indices[column_prefix + 'sils_means'] = [np.mean([x[i] for x in main_ss_vec]) for i in range(7)]
    df_cluster_indices[column_prefix + 'sils_sd'] = [np.std([x[i] for x in main_ss_vec], ddof=1) for i in range(7)]
    df_cluster_indices[column_prefix + 'ch_means'] = [np.mean([x[i] for x in main_ch_vec]) for i in range(7)]
    df_cluster_indices[column_prefix + 'ch_sd'] = [np.std([x[i] for x in main_ch_vec], ddof=1) for i in range(7)]
    df_cluster_indices[column_prefix + 'db_means'] = [np.mean([x[i] for x in main_db_vec]) for i in range(7)]
    df_cluster_indices[column_prefix + 'db_sd'] = [np.std([x[i] for x in main_db_vec], ddof=1) for i in range(7)]
    df_cluster_indices[column_prefix + 'inert_means'] = [np.mean([x[i] for x in main_inertia_vec]) for i in range(7)]
    df_cluster_indices[column_prefix + 'inert_sd'] = [np.std([x[i] for x in main_inertia_vec], ddof=1) for i in range(7)]

    # phenotype calling
    phen = ['e', 'he', 'hm', 'm']  # increasing order of ZEB
    for run_index in range(5):
        df_list[run_index]['clus_labels'] = final_labels[run_index]  # add the original cluster labels (not needed)
        df_list[run_index]['phenotype'] = 'Not Computed'  # placeholder until the phenotypes are assigned
        df = df_list[run_index]
        zeb_medians = []  # calculate the ZEB medians of the 4 clusters
        for c in range(4):
            zeb_medians.append(np.median(df.loc[df['clus_labels'].to_numpy() == c, 'ZEB']))
        clus_labels = [0, 1, 2, 3]
        sorted_zeb_medians = sorted(zeb_medians)
        for i in range(4):
            zeb_value = sorted_zeb_medians[i]  # in increasing order
            clus_label_value = clus_labels[zeb_medians.index(zeb_value)]
            df.loc[df['clus_labels'].to_numpy() == clus_label_value, 'phenotype'] = phen[i]  # assign phenotype
        df.to_csv(path_to_folder + '/' + psfname + str(run_index + 1) + '_clusData.csv', index=False)
    print('Time for last folder: ' + str(round(time.time() - t_per_folder, 2)) + ' seconds.')

df_cluster_indices.to_csv(save_paths[0] + '/' + 'cluster_indices.csv', index=False)
print('\nAll finished. It took ' + str(round((time.time() - t_start) / 3600, 2)) + ' hours.')
