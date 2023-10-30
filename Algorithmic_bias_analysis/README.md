# Algorithmic Bias analysis

The following README explains a bit how to reproduce the results in the algorithmic bias section of the paper. 

## Code

The files are the following:
- `Prepare_dataset.R`: Filters and prepares the dataset from Pyke et al. used in the analysis. Namely, peptides with any chemical modifications are removed, and Negatives are sampled from the human proteome. 
- `Motif_generation.R`: Motif generation for calculating FOOP scores (see paper for more details)
- `Data_analysis.R`: PPV/FOOP analysis and plot generation 

The data for the PPV/FOOP analysis and for generating the numbers and the plots can be found in the following link:
https://rice.box.com/s/4m8gvvooy67kcmrsakyyqhsmi3s2u3nj

P.S.: In order to generated NetMHCpan4.1 and MHCFlurry scores by yourself, you will need to install these two tools first.
P.S.2: Reach out to the authors or post an issue in case you find something wrong in the analysis!
