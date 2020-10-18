import numpy as np
import pandas as pd
import os
from sklearn.decomposition import PCA

helpers = __import__('helpers')
#  To compute proportions (and the respective standard deviations) for PCA and the link strength analysis in the four
#  phenotypes

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
parameters = ['PCA_component_1', 'PCA_component_2', 'PCA_explained_var_1', 'PCA_explained_var_2',
              'zeb_u200_asym_e', 'zeb_u200_asym_he', 'zeb_u200_asym_hm', 'zeb_u200_asym_m',
              'lin28_let7_asym_e', 'lin28_let7_asym_he', 'lin28_let7_asym_hm', 'lin28_let7_asym_m',
              'total_coupling_e', 'total_coupling_he', 'total_coupling_hm', 'total_coupling_m']
# each of the PCA component is a space separated (in order of the columns) coefficients for each of the genes

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

    # PCA
    allgenes = ['u200', 'ZEB', 'SNAIL', 'LIN28', 'let7', 'NFkB', 'GRHL2', 'KEAP1', 'NRF2',
                'ECad', 'OVOL']  # all genes available in all of the analysis
    # doing for only 1 replicate (run 1)
    df_of_interest = df_list[0].loc[:, [x for x in df_list[0].columns if x in allgenes]]
    # df with only gene columns (log2 Z-Normalized (oe0 reference))
    pca_object = PCA(svd_solver='full')
    pca_object.fit(df_of_interest.to_numpy())
    pca1 = pca_object.components_[0, :]
    pca2 = pca_object.components_[1, :]
    pca_var1 = pca_object.explained_variance_ratio_[0]
    pca_var2 = pca_object.explained_variance_ratio_[1]

    append_vector_means.append(' '.join([str(x) for x in pca1]))
    append_vector_means.append(' '.join([str(x) for x in pca2]))
    append_vector_means.append(pca_var1)
    append_vector_means.append(pca_var2)

    append_vector_sds += [np.nan, np.nan, np.nan, np.nan]

    # link strength analysis
    # zeb_u200_asym
    run_wise_asyms = []
    for run_index in range(5):
        # getting the b1, b2, asym indices
        df = df_param_list[run_index].loc[df_param_list[run_index]['num_states'] == 1, :]  # only monostable systems
        df2 = df_list[run_index].loc[df_list[run_index]['num_states'] == 1, :]
        y1 = df['Inh_of_u200ToZEB'].to_numpy()
        mu1 = df['Trd_of_u200ToZEB'].to_numpy()
        bigG1 = df['Prod_of_u200'].to_numpy()
        k1 = df['Deg_of_u200'].to_numpy()
        b1 = helpers.bond_strength(y1, mu1, bigG1, k1, 'Inh')
        y2 = df['Inh_of_ZEBTou200'].to_numpy()
        mu2 = df['Trd_of_ZEBTou200'].to_numpy()
        bigG2 = df['Prod_of_ZEB'].to_numpy()
        k2 = df['Deg_of_ZEB'].to_numpy()
        b2 = helpers.bond_strength(y2, mu2, bigG2, k2, 'Inh')
        asym = helpers.asymmetry_bond_strength(b2, b1)  # ZEB_u200/u200_ZEB, i.e. higher asym means ZEB suppresses u200
        phen_wise_asym = []
        for p in ['e', 'he', 'hm', 'm']:
            models = df2.loc[df2['phenotype'] == p, 'model_id'].to_numpy()  # all monostable model ids of p
            bool_mask = (df['model_id'].isin(models))  # boolean mask for all param sets corresponding to the models
            phen_wise_asym.append(np.mean(asym[bool_mask]))
        run_wise_asyms.append(phen_wise_asym)

    append_vector_means += [np.mean([x[i] for x in run_wise_asyms]) for i in range(4)]
    append_vector_sds += [np.std([x[i] for x in run_wise_asyms], ddof=1) for i in range(4)]



    # lin28_let7_asym
    run_wise_asyms = []
    for run_index in range(5):
        # getting the b1, b2, asym indices
        df = df_param_list[run_index].loc[df_param_list[run_index]['num_states'] == 1, :]
        df2 = df_list[run_index].loc[df_list[run_index]['num_states'] == 1, :]
        y1 = df['Inh_of_LIN28Tolet7'].to_numpy()
        mu1 = df['Trd_of_LIN28Tolet7'].to_numpy()
        bigG1 = df['Prod_of_LIN28'].to_numpy()
        k1 = df['Deg_of_LIN28'].to_numpy()
        b1 = helpers.bond_strength(y1, mu1, bigG1, k1, 'Inh')
        y2 = df['Inh_of_let7ToLIN28'].to_numpy()
        mu2 = df['Trd_of_let7ToLIN28'].to_numpy()
        bigG2 = df['Prod_of_let7'].to_numpy()
        k2 = df['Deg_of_let7'].to_numpy()
        b2 = helpers.bond_strength(y2, mu2, bigG2, k2, 'Inh')
        asym = helpers.asymmetry_bond_strength(b1, b2)  # LIN28_let7/let7_LIN28
        phen_wise_asym = []
        for p in ['e', 'he', 'hm', 'm']:
            models = df2.loc[df2['phenotype'] == p, 'model_id'].to_numpy()
            bool_mask = (df['model_id'].isin(models))
            phen_wise_asym.append(np.mean(asym[bool_mask]))
        run_wise_asyms.append(phen_wise_asym)

    append_vector_means += [np.mean([x[i] for x in run_wise_asyms]) for i in range(4)]
    append_vector_sds += [np.std([x[i] for x in run_wise_asyms], ddof=1) for i in range(4)]

    # total coupling
    run_wise_combs = []
    for run_index in range(5):
        # getting the b1, b2, comb_strength indices
        df = df_param_list[run_index]
        df2 = df_list[run_index]
        y1 = df['Inh_of_u200ToLIN28'].to_numpy()
        mu1 = df['Trd_of_u200ToLIN28'].to_numpy()
        bigG1 = df['Prod_of_u200'].to_numpy()
        k1 = df['Deg_of_u200'].to_numpy()
        b1 = helpers.bond_strength(y1, mu1, bigG1, k1, 'Inh')
        y2 = df['Inh_of_let7ToZEB'].to_numpy()
        mu2 = df['Trd_of_let7ToZEB'].to_numpy()
        bigG2 = df['Prod_of_let7'].to_numpy()
        k2 = df['Deg_of_let7'].to_numpy()
        b2 = helpers.bond_strength(y2, mu2, bigG2, k2, 'Inh')
        comb_strength = np.log2(b1 * b2)
        phen_wise_comb = []
        for p in ['e', 'he', 'hm', 'm']:
            models = df2.loc[df2['phenotype'] == p, 'model_id'].to_numpy()
            bool_mask = (df['model_id'].isin(models))
            phen_wise_comb.append(np.mean(comb_strength[bool_mask]))
        run_wise_combs.append(phen_wise_comb)
    append_vector_means += [np.mean([x[i] for x in run_wise_combs]) for i in range(4)]
    append_vector_sds += [np.std([x[i] for x in run_wise_combs], ddof=1) for i in range(4)]

    df_tables_numbers_means[psfname + oe_level] = append_vector_means
    df_tables_numbers_sds[psfname + oe_level] = append_vector_sds

df_tables_numbers_means.to_csv(save_paths[0] + '/pca_link_strength_means.csv', index=False)
df_tables_numbers_sds.to_csv(save_paths[0] + '/pca_link_strength_sds.csv', index=False)
