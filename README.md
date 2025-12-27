# Replication of a Sjögren's Syndrome Genome-Wide Association Study (GWAS) 

## 🧬 Project Overview 
Primary Sjögren's Syndrome (pSS) is a complex autoimmune disease that is characterized by inflammation and tissue destruction of the salivary and lacrimal glands, though many patients also experience extraglandular symptoms including synovitis and interstitial lung disease (PMID: 36914790). The etiology of the disease is complex, with a combination of genetic and environmental factors playing a role in disease risk and pathogenesis. Twin studies show a monozygotic twin concordance rate of ~30%, indicating a modest but important role of genetics in disease risk. Prior Genome-Wide Association Studies (GWAS) have identified several key risk loci for pSS, including  _IRF5-TNPO3_, _TYK2_, and _STAT4_. 

This project reproduces a published GWAS of pSS (PMID: 28076899) and demonstrates how population genetics results can be translated into hypotheses for target validation and pathway prioritization in autoimmune disease.

## Project Goals 
- Quality Control (QC) of GWAS summary statistics and replication of significant associations
- Replicate informative plots: Quantile-Quantile (QQ) plot, Manhattan Plot 
- Identify genome-wide significant loci for pSS
- Propose pre-clinical strageies to validate the identified targets.
  
## Why is this Approach Useful for Drug Development? 
Population genetics approaches nominate specific genes and pathways that are important in the genetic risk for complex diseases such as pSS. For drug development, this information can be used to inform the development of pre-clinical platforms (_in vitro/in vivo/ex vivo_ models) to validate potential disease targets and dissect immune pathways. Since autoimmune diseases share common pathways, identifying risk loci for one autoimmune disease is often informative for others. For example, the _IRF5_ locus is also highly associated with risk for Systemic Lupus Erythematosus, so validation of this target could also guide  indication selection decisions. 

## Data Availability
- Sjögren's Syndrome GWAS data is pulled from: Taylor KE, et al. Genome-Wide Association Analysis Reveals Genetic Heterogeneity of Sjögren's Syndrome According to Ancestry. Arthritis Rheumatol. 2017 Jun;69(6):1294-1305. doi: 10.1002/art.40040.
- Genetic Coordinates for the IRF5-TNPO3 Locus is from: Kottyan LC, et al. The IRF5–TNPO3 association with systemic lupus erythematosus has two components that other autoimmune disorders variably share. Hum Mol Genet. 2015 Jan 15;24(2):582-96. doi: 10.1093/hmg/ddu455.

Raw and processed data from the publication are **not distributed** in this repository due to licensing and size restrictions.  
Users must obtain data directly from the original publication and associated repositories.

## Dataset Overview
The GWAS data in this project are from Taylor KE, et al and are availble under accession GCST012796 in the GWAS Catalog. The dataset contains 585 pSS cases and 1546 controls (both of European Ancestry). The study design analytical strategy are available in the methods section of the publication: https://acrjournals.onlinelibrary.wiley.com/doi/10.1002/art.40040. 

## Project Organization (Notebooks written in R Markdown) 
- Notebook 1: Download and Prepare GWAS Data 
- Nobteook 2: Harmonize and Clean GWAS Data
- Notebook 3: GWAS Visualization
- Notebook 4: _IRF5_ Locus Analysis

## Main Findings (also under Results) 
- QQ Plot: https://github.com/samvirolainen/Role_of_IRF5_in_Sjogrens_Syndrome/blob/main/results/qq_sjogrens.png
- Manhattan Plot: https://github.com/samvirolainen/Role_of_IRF5_in_Sjogrens_Syndrome/blob/main/results/manhattan_sjogrens.png

## Interpretation of GWAS Plots 
- QQ plot: most points lie on the diagonal, but at the tail there is a strong deviation from the diagonal to indicate true genetic associations.
- Manhattan plot: Coherent linkage disequilibrium (LD) blocks form with little random noise, indicating strong agreement among the genetic signals.

## Replicated pSS Risk Loci 

|Locus           | Chromosome | Biological Interpretation | Possible Next Steps for Target Validation 
| ------------- | ------------- | ----------------------- | ----------------- | 
| _HLA_  | chr 6  | Encodes Major Histocompatibility Complex (MHC), a key component of antigen presentation in the immune system, highly associated with many autoimmune diseases | Examine how _HLA_ risk haplotypes alter other molecular phenotypes, such as protein abundance (through protein Quantitative Trait Loci, pQTLs). Other studies have associated HLA risk haplotypes with protein abundance: https://www.nature.com/articles/s41467-024-50583-8 , though colocalization methods for the _HLA_ region remain challenging to implement. 
| _IRF5-TNPO3_ | chr 7  | Encodes IRF5, a transcription factor and key regulator of the Type I Interferon (IFN) pathway, a clinically validated innate immune pathway that drives several autoimmune diseases including pSS | Generate _in vitro_ assays to assess impact of _IRF5_ haplotypes on IFN pathway. A variant-by-variant approach could utilize a luciferase reporter system, while a high-throughput approach (such as a Massively Parallel Reporter Assay) could screen several variants at a time. 
| _STAT4_ | chr 2  | Encodes STAT4, a transcription factor that regulates gene expression programs downstream of several immune cytokine signaling pathways | Obtain primary immune cells from healthy donors with and without _STAT4_ risk haplotypes and stimulate the cells with cytokines that activate _STAT4_ and readout activation of these pathways (such as at the RNA and protein level). Readouts impacted by  _STAT4_ pSS risk haplotypes could serve as more specific drug targets and be more beneficial to modulate than STAT4 itself. 

## Environment
This project uses a reproducible conda environment (see `environment.yml`):

```bash
conda env create -f environment.yml
conda activate coloc_r_env
```
## License
This project’s code is released under the MIT License (see `LICENSE`).
