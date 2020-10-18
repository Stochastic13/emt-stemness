import os
from shutil import copyfile

# Collates the raw_data into a single file (with each stable state a separate row with num_states and model_id)
# Collates all the parameter sets into a single file (with each set in a separate row with num_states and model_id)

paths = ['GRHL2/full/oe0', 'GRHL2/full/oe10', 'GRHL2/full/de10',
         'GRHL2/partial/oe0', 'GRHL2/partial/oe10', 'GRHL2/partial/de10',
         'NRF2/oe0', 'NRF2/oe10', 'NRF2/de10',
         'OVOL/oe0', 'OVOL/oe10', 'OVOL/de10',
         'base/oe0', 'uncoupled']
raw_data_paths = ['raw_data/' + x for x in paths]  # where the raw output (dat files) is stored
save_paths = ['collated_data/' + x for x in paths]  # path of the folder to save in

for path_to_folder, save_path in zip(raw_data_paths, save_paths):
    print('\nProcessing Data in ' + path_to_folder + '; Saving files in ' + save_path)
    # list all files in the folder
    files = [f for f in os.listdir(path_to_folder) if os.path.isfile(path_to_folder + '/' + f)]
    # extract circuit names based on existing .topo files
    name_circuit = [x.split('.')[0] for x in files if '.topo' in x]
    print('Found topo files for: ' + ', '.join(name_circuit))
    # filter files of relevance
    relevant_files = dict([(x, []) for x in name_circuit])
    for f in files:
        for name in name_circuit:
            if name in f:  # assuming no name is a subset of another name
                relevant_files[name].append(f)
                break

    # parse data for each circuit separately and output a csv
    for name, files in relevant_files.items():
        # extract gene_names
        prefix_filename = [x.split('.')[0] for x in files if x.split('.')[1] == 'cfg']
        # Above line needed because the prefix of filename changes across unperturbed and perturbed
        assert len(prefix_filename) == 1, 'Zero/More than one cfg file(s) associated with one name'
        with open(path_to_folder + '/' + prefix_filename[0] + '.cfg', 'r') as f:
            rd = f.read()
            rd = rd.split('\n')
            rd = [x.split('\t') for x in rd if len(x) > 0]  # skip empty rows
            ngenes = int([x[1] for x in rd if x[0] == 'NumberOfGenes'][0])  # value corresponding to NumberOfGenes
            gene_names = [x[1] for x in rd[31:(31 + ngenes)]]  # all gene_names begin from line 31 in cfg (0-indexed)
        print('Genes for ' + name + ' : ' + ', '.join(gene_names))
        with open(path_to_folder + '/' + prefix_filename[0] + '.prs', 'r') as f:  # load parameter names
            rd = f.read()
            rd = rd.split('\n')
            colnames = ['model_id', 'num_states'] + [x.split('\t')[0] for x in rd[1:] if len(x) > 0]
        copyfile(path_to_folder + '/' + prefix_filename[0] + '.prs', save_path + '/' + name + '.prs')  # copy prs file
        output_prm = ','.join(colnames) + '\n'  # load column names for the parameter csv file
        with open(path_to_folder + '/' + prefix_filename[0] + '_parameters.dat', 'r') as f:  # load the parameters
            rd = f.read()
            rd = rd.split('\n')
            rd = [','.join(x.split('\t')) for x in rd if len(x) > 0]  # make each row into a csv row
            output_prm += '\n'.join(rd)
        with open(save_path + '/' + name + '_parameters_formatted.csv', 'w', newline='') as f:  # save the file
            f.write(output_prm)
        output = ','.join(['model_id', 'state_id', 'num_states'] + gene_names) + '\n'  # add column names
        for filename in files:
            if 'solution' in filename:
                num_states = filename.split('_')[-1].split('.')[0]  # extract solution number
                with open(path_to_folder + '/' + filename, 'r') as f:  # load file as nested list
                    rd = f.read()
                    rd = rd.split('\n')
                    rd = [x.split('\t') for x in rd if len(x) > 0]
                for row in rd:
                    for state_id in range(int(num_states)):
                        # state_id: 0-indexed identifier of the diff states in a multi-state solution
                        output += row[0] + ','  # add model_id
                        output += str(state_id) + ',' + row[1] + ','  # add state_id and num_states
                        output += ','.join(row[(2 + ngenes * state_id): (2 + ngenes * state_id + ngenes)]) + '\n'
                        # add solution vectors (log2 scale)
        with open(save_path + '/' + name + '_collated.csv', 'w', newline='') as f:  # save the file
            f.write(output)
    print('Done!')
