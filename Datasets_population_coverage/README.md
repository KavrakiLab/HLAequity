# Pan-allele predictors - dataset bias 

Here we calculate the population coverage of different training datasets (NetMHCpan4.1 and MHCFlurry2.0 - binding affinity and mass spectrometry data). 

### 1. Download and prepare the datasets

Find the detailed instructions and notebooks for downloading the datasets inside the `./datasets/` folder.

    - Prepare the NetMHCpan4.1 datasets 
    - Prepare the MHCFlurry2.0 datastes
    - Visualize the dataset content 
    
### 2. Calculate the population coverage

Use the `./CalculateCoverage_IEDB.ipynb` notebook to calculate the population coverage of the datasets. This notebook will guide you through the following steps:

    - Downloading and setting up the IEDB population coverage tool from [here](http://tools.iedb.org/population/download/)
    - Running the population coverage for a single population and a single dataset
    - Running the population coverage for all datasets and all populations
    
    
### 3. Visualize the population coverage of the datasets

Run the `./Visualize_results.ipynb` notebook to see how different datasets cover the given populations and how this coverage relates to different income levels and ancestries. 