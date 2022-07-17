This code is in adjunct to the paper linked below.

**Satwik Pasani; Sarthak Sahoo; Mohit Kumar Jolly; *Hybrid E/M phenotype(s) and stemness:
a mechanistic connection embedded in network topology*, 2020**

[paper link](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7794989/)

The `raw_data` zip file must be extracted first. It contains a nested folder organization according to the circuit and the level of expression (i.e. overexpression (oe10), downexpression (de10) or unperturbed (oe0)). Each final folder has the 5 `.topo` files for the 5 replicates of the circuit (except the uncoupled circuit), with one `scriptfile_all.sh` file per folder containing the bash script to run RACIPE on all the `.topo` files in the folder. The terminal outputs of the 5 runs, as evident from the bash script, is stored in the `termoutput` files. The actual outputted files are not included due to their large size.

The remaining files are the python scripts to reproduce the tables/figures created in the analysis. Before running, two empty folders titled `collated_data` and `figures` must be created in the working directory. `figures` must contain the subfolders `plots` and `tables_numbers` and `clustering`. The order in which the scripts must ideally be run is as follows: `collate_data` -> `normalize_data` -> `cluster_data` -> `tables_and_numbers_1/2/3` -> `statistical_significance_1/2` -> `figures_1/2/3/misc` -> `additional_analysis/2`

