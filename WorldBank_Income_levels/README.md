# Map the populations and their countries to the income levels

1. We take download the country income levels from the World Bank [here](https://datahelpdesk.worldbank.org/knowledgebase/articles/906519-world-bank-country-and-lending-groups) (current classification by income in XLSX format).
2. We convert the xlsx format to the csv format with comma as separator (using excel or other tools).
3. We take the populations found in the AFND and their mapping to countries and map the populations to the income levels. Related code can be found in the `Map_populations.ipynb` notebook. 
4. The result is the `population_income_map.csv` that contains the resulting mapping.