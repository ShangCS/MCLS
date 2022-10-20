# Multi-omics clustering for cancer subtyping based on latent subspace learning
## Table of contents
* [Descriptions](#descriptions)
* [Preparations](#preparations)
## Descriptions
### 1.Realworld datasets
All realworld datasets are saved in ``` "data/..." ``` including the following seven multi-omics datasets:
* BIC
* COAD
* GBM
* KRCCC
* LSCC
* Bladder
* Brain
### 2.Code files descriptions
#### 1.Complete_data.R
Complete multi-omics datasets clusteirng analysis. The default test data is set as the ```BIC```.
#### 2.Incomplete_data.R
Incomplete multi-omics datasets clusteirng analysis. The default test data is set as the ```BIC```.
#### 3.run_data_col.R
Complete multi-omics data integration and get the embedding result.
#### 4.specClust.R
The spetcral clusteing method applied in the latent subspace to get the results.
#### 5.enrichment.R
Using clinical data to do enrichment analysis.


## Preparations
(1) Download and import the datasets.

(2) Run ```Complete_data.R``` or ```Incomplete_data``` file.


