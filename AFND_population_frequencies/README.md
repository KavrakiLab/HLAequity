# AFND database - getting the population frequencies

This folder contains all the scripts neccessary to download and prepare the AFND population frequencies used in the downstream analysis.

### 1. Download the AFND allele frequencies per locus

You can use the `AFND_get_data_per_locus.ipynb` notebook to downlaod the data from  http://www.allelefrequencies.net/ and clean this data per locus. This notebook contains the following steps:

1. We query the AFND: http://www.allelefrequencies.net/ using the basic AFND API given [here](http://www.allelefrequencies.net/extaccess.asp)
2. We format the allele codes according to the standard nomenclature (bellow). We keep the alleles which contain the _allele group_ and the _specific HLA protein_ fields.
3. We map the populations to countries.
4. We remove the duplicate entries for "population"-"allele" pairs.
5. We perform sanity checks on the data (we make sure that the allele frequencies per population add up to 1)

**Output:** As the output of this notebook you will find the csv files `AFND_data_locus_*.csv` containing the population frequencies of different alleles across the chosen loci. 
 

### 2. Combine the AFND frequencies across different loci

Use the `AFND_combine_loci_data.ipynb` notebook to combine the AFND data for different loci (i.e., A/B/C). 


1. We load the per-locus data
2. We clean this data keeping only the populations that have information for each of the selected loci

**Output:** The final allele-frequency data is saved to the `AFND_data_locus_all.csv` file. 


### 3. Visualize and analyze the AFND data (optional)


Use the `AFND_visualize_populations.ipynb` notebook to visualize the data related to the allele frequencies across the populations.

1. View the raw table with the AFND data
2. Select and compare populations and their allele frequencies 
3. View the dendogram and a heatmap to see which populations have similar allele frequencies
4. Inspect the 2D umap embeddings on the populations and select populations in this embedding to plot their allele frequencies

