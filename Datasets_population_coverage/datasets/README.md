# Pan-allele predictors - get and preprocess the datasets

This folder contains all the scripts neccessary to download and prepare datasets used by MHCFlurry2.0 and NetMHCpan4.1. 

### NetMHCpan4.1 datasets

Download and untar the dataset used for training NetMHCpan4.1 from [here](https://services.healthtech.dtu.dk/suppl/immunology/NAR_NetMHCpan_NetMHCIIpan/). Place the content inside the folder `./NetMHCpan_train`. Follow the instructions in `./NetMHCpan4.1_prepare_data.ipynb` notebook to prepare and clean the data by:

1. Reading in the list of all alleles given by `allelelist` file.
2. Loading and extracting the loci for all the binding affinity datasets in `c00*_ba` files.
3. Loading and extracting the loci for all the eluted ligand datasets in `c00*_el` files. Here we load only the monoallelic portions of the data.
4. Deconvoluting the multiallelic portion of the dataset.
5. Merging the monoalleling and multiallelic portions of the datasets. 
6. Visualizing the resulting data.


### MHCFlurry2.0 datasets

Download Data_S3, Data_S5 and Data_S6 files from the supplementary material found [here](https://data.mendeley.com/datasets/zx3kjzc3yx/3). Place downloaded files inside the `./mhcflurry` folder. Follow the steps in `./MHCFlurry_prepare_data.ipynb`to prepare the data:

1. Read in the binding affinity data given by `Data_S3.csv`.
2. Read in the monoallelic mass spec data given by `Data_S3.csv` and `Data_S5.csv`.
3. Deconvolute the multiallelic data from `Data_S6.csv`.
4. Combine the multiallelic and monoallelic data.
5. Visualize the resulting data.

### Visualize dataset contents

Check the interactive dashboard with the dataset contents here: `./Visualize_datasets.ipynb`