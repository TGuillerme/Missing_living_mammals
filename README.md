# Assessment of available anatomical characters for phylogenetic analysis among living mammals
[Thomas Guillerme](http://tguillerme.github.io) and [Natalie Cooper](http://nhcooper123.github.io/).

This repository contains all the code and data used in the manuscript.
###### Manuscript in review for Biology Letters.
A preprint is available [here](http://biorxiv.org/content/early/2015/07/28/022970).

-------

## Data <a href="http://figshare.com/articles/Assessment_of_cladistic_data_availability_for_living_mammals/1575729"><img src="http://tguillerme.github.io/images/logo-FS.png" height="26" widht="26"/></a> 
All the raw data (cladistic matrices) and transformed data (cladistic matrices with corrected taxonomy entries) are abailable in the [data folder](https://github.com/TGuillerme/Missing_living_mammals/tree/master/Data) on this repository. The raw data are also avaialble on [figshare](http://figshare.com/articles/Assessment_of_cladistic_data_availability_for_living_mammals/1575729). 

## Analysis
All the individual functions for this analysis (and their testing) are avaiable in the [function folder](https://github.com/TGuillerme/Missing_living_mammals/tree/master/Functions).

The analysis is divided into three steps
* 1.Extracting the living taxa from the matrices
* 2.Analysing data availability and distribution per order
* 3.Summarising the results

Note that the first steps take some time (several hours on your usual desktop machine) so the results from this steps is made directly available [here](https://github.com/TGuillerme/Missing_living_mammals/tree/master/Data/List_of_matching_taxa).

#### 1-Extracting the living taxa from the matrices
This step sorts all the OTUS from all the matrices into either `living` or `fossil` taxa in order to have the list of all the available mammals with morphological (cladistic) data. The code for this step is available [here](https://github.com/TGuillerme/Missing_living_mammals/blob/master/Analysis/Extracting_living_taxa.R).

#### 2-Analysing data availability and distribution per order
The second step analysis the data availability within each mammalian order. It calculates the proportion of living mammals with avaiable cladistic data as well as the distribution of this data accross the phylogeny. The code for this step is available [here](https://github.com/TGuillerme/Missing_living_mammals/blob/master/Analysis/Analysis_per_order.R).

#### 3-Summarising the results
Finally, the third step summarises the results of the analysis in step 2 and generates the manuscript figure/table (fable!) containing all the information of the availability and distribution of living mammals cladistic data. The code for this step is available [here](https://github.com/TGuillerme/Missing_living_mammals/blob/master/Analysis/Results_summary.R).
Note that the output is a table in `LaTeX` format.

#### Supplementary analysis
* [Google Scholar search rarefaction](https://github.com/TGuillerme/Missing_living_mammals/blob/master/Analysis/Rarefaction.R)
* [Sampling effort](https://github.com/TGuillerme/Missing_living_mammals/blob/master/Analysis/Sampling_effort.R)
