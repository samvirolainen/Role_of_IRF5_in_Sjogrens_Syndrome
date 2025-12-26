# Replication of a Sjögren's Syndrome Genome-Wide Association Study (GWAS) 

## 🧬 Project Overview 
Primary Sjögren's Syndrome (pSS) is a complex autoimmune disease that is characterized by inflammation and tissue destruction of the salivary and lacrimal glands, though many patients also experience extraglandular symptoms including synovitis and interstitial lung disease (PMID: 36914790). The etiology of the disease is complex, with a combination of genetic and environmental factors playing a role in disease risk and pahtogenesis. Twin studies show a monozygotic twin concordance rate of ~30%, indicating a modest but important role of genetics in disease risk. Prior Genome-Wide Association Studies (GWAS) have identified several key risk loci for pSS, including  _IRF5-TNPO3_, _TYK2_, and _STAT4_. This goal of this project is to replicate a published GWAS (PMID: 28076899) for pSS. 

## Project Goals 
- Quality Control (QC) of GWAS summary statistics and replication of significant associations
- Replicate informative plots: Quantile-Quantile (QQ) plot, Manhattan Plot 
- Identify genome-wide signficant loci for pSS 

## Why is this Approach Useful for Drug Development? 
Population genetics approaches nominate specific genes and pathways that are important in the genetic risk for complex diseases such as pSS. For drug development, this information can be used to inform the development of pre-clinical platforms (in vitro/in vivo/ex vivo models) to validate potential disease targets and pathways. Since autoimmune diseases share common pathways, identifying risk loci for one autoimmune disease is often informative for others. For example, the _IRF5_ locus is also highly associated with risk for Systemic Lupus Erythematosus.  

## Data Availability
- Sjögren's Syndrome GWAS data is pulled from: Taylor KE, et al. Genome-Wide Association Analysis Reveals Genetic Heterogeneity of Sjögren's Syndrome According to Ancestry. Arthritis Rheumatol. 2017 Jun;69(6):1294-1305. doi: 10.1002/art.40040.
- Genetic Coordinates for the IRF5-TNPO3 Locus is from: Kottyan LC, et al. The IRF5–TNPO3 association with systemic lupus erythematosus has two components that other autoimmune disorders variably share. Hum Mol Genet. 2015 Jan 15;24(2):582-96. doi: 10.1093/hmg/ddu455.

Raw and processed data from the publication are **not distributed** in this repository due to licensing and size restrictions.  
Users must obtain data directly from the original publication and associated repositories.

## Dataset Overview
The GWAS data in this project are from Taylor KE, et al and are availble under accession GCST012796 in the GWAS Catalog. The dataset contains 585 pSS cases and 1546 controls of European Acestry. The study design analytical strategy are available in the methods section of the publication: https://acrjournals.onlinelibrary.wiley.com/doi/10.1002/art.40040. 

## Project Organization (Notebooks written in R Markdown) 
- Notebook 1: Download and Prepare GWAS
- Nobteook 2: Harmonize and Clean GWAS Data
- Notebook 3: GWAS Visualization
- Notebook 4: _IRF5_ Locus Analysis

## Main Findings 
https://github.com/samvirolainen/Role_of_IRF5_in_Sjogrens_Syndrome/blob/main/results/qq_sjogrens.png

## Environment
This project uses a reproducible conda environment (see `environment.yml`):

```bash
conda env create -f environment.yml
conda activate coloc_r_env
```
## License
This project’s code is released under the MIT License (see `LICENSE`).
