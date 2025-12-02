# Colocalization of Autoimmune Risk Loci with Quantitative Trait Loci (QTLs) 

## 🧬 Scientific Rationale 
Autoimmune diseases such as Lupus (SLE), Juvenile Idiopathic Arthritis (JIA), and Sjögren's Syndrome result from a complex combination of genetic and environmental facotrs. Because of this complexity, identifying drug targets and targetable pathways for autoimmune diseases is often difficult. Genetic studies can provide clues of the specific genes, pathways, and cell-specific mediators of a specific disease. This project leverages human genetic data to 1) Identify risk genes and pathways driving risk for autoimmune diseases 
and 2) Propose how these risk genes contribute to disease and in which cell types they are important 

## Project Goals 
- Identify and quality control Genome-Wide Association Study (GWAS) data for specific autoimmune diseases and perform fine mapping to identify causal risk variants 
- Identify Quantitative Trait Loci (QTL) data (matched by genetic ancestry and genome bulild to the GWAS data) from relevant tissues and perform fine mapping to identify the variants driving the QTL signal 
- Perform colocalization analysis and identify risk variants driving genotype-dependent gene expression (eQTL) and protien abundance (pQTL)

## Why is this Approach Useful for Target Identification and Drug Development? 
This approach nominates specific immune pathways that are important for autoimmmune disease pathogenesis and prioritizes relevant cell-types in these pathways. In early discovery, this information can be used to inform the development of pre-clinical platforms (in vitro/in vivo/ex vivo models) to validate potential disease targets. 

## Data Availability
-Sjögren's Syndrome GWAS data is pulled from: Taylor KE, et al. Genome-Wide Association Analysis Reveals Genetic Heterogeneity of Sjögren's Syndrome According to Ancestry. Arthritis Rheumatol. 2017 Jun;69(6):1294-1305. doi: 10.1002/art.40040.
-Genetic Coordinates for the IRF5-TNPO3 Locus is from: Kottyan LC, et al. The IRF5–TNPO3 association with systemic lupus erythematosus has two components that other autoimmune disorders variably share. Hum Mol Genet. 2015 Jan 15;24(2):582-96. doi: 10.1093/hmg/ddu455.

Raw and processed data from the publication are **not distributed** in this repository due to licensing and size restrictions.  
Users must obtain data directly from the original publication or associated repositories.

## Environment
This project uses a reproducible conda environment (see `environment.yml`):

```bash
conda env create -f environment.yml
conda activate coloc_r_env
```
## License
This project’s code is released under the MIT License (see `LICENSE`).
