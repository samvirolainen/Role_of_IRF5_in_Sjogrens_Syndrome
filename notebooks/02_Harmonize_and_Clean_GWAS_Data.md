02 — Harmonize & Clean GWAS Data
================
Sam Virolainen
2025-11-26

## Load GWAS Data

``` r
gwas <- readRDS("../data/processed/harmonized/gwas_sjogrens_raw.rds")
head(gwas)
```

    ##    chromosome base_pair_location effect_allele odds_ratio standard_error
    ##         <int>              <int>        <char>      <num>          <num>
    ## 1:          1             752566             C     0.9576        0.09064
    ## 2:          1             752721             T     0.9825        0.09870
    ## 3:          1             754063             T     1.0010        0.15770
    ## 4:          1             759036             A     0.9477        0.17140
    ## 5:          1             771967             A     0.9775        0.11330
    ## 6:          1             774047             A     1.1410        0.91840
    ##    ci_lower ci_upper p_value          beta
    ##       <num>    <num>   <num>         <num>
    ## 1:   0.8017    1.144  0.6326 -0.0433251247
    ## 2:   0.8097    1.192  0.8581 -0.0176549352
    ## 3:   0.7350    1.364  0.9941  0.0009995003
    ## 4:   0.6774    1.326  0.7541 -0.0537172825
    ## 5:   0.7829    1.221  0.8409 -0.0227569871
    ## 6:   0.1887    6.905  0.8854  0.1319050709

## Normalize Chromosome Formatting

``` r
gwas$chromosome <- as.character(gwas$chromosome)
gwas$chromosome[gwas$chromosome == "X"] <- "23"
gwas$chromosome[gwas$chromosome == "Y"] <- "24"
gwas$chromosome <- as.integer(gwas$chromosome)
```

## Normalize Allele Format

``` r
gwas$effect_allele <- toupper(gwas$effect_allele)
gwas$other_allele <- toupper(gwas$other_allele)
```

## Remove non-SNPs and Invalid Genotype Calls

``` r
valid <- c("A", "C", "G", "T")
gwas <- gwas %>%
filter(effect_allele %in% valid)
```

## Create Standardized Variant ID

``` r
gwas$variant_id <- paste0(
"chr", gwas$chromosome, ":",
gwas$base_pair_location, "_",
gwas$effect_allele
)
```

## Summarize Statistics After Filtering

``` r
summary(gwas)
```

    ##    chromosome     base_pair_location  effect_allele        odds_ratio       
    ##  Min.   : 1.000   Min.   :     2220   Length:1442662     Min.   :0.000e+00  
    ##  1st Qu.: 4.000   1st Qu.: 30312648   Class :character   1st Qu.:1.000e+00  
    ##  Median : 8.000   Median : 68236469   Mode  :character   Median :1.000e+00  
    ##  Mean   : 9.158   Mean   : 77750434                      Mean   :2.791e+06  
    ##  3rd Qu.:14.000   3rd Qu.:114732824                      3rd Qu.:1.000e+00  
    ##  Max.   :23.000   Max.   :249212725                      Max.   :1.265e+10  
    ##  standard_error         ci_lower         ci_upper         p_value      
    ##  Min.   :    0.062   Min.   :0.0000   Min.   :0.4533   Min.   :0.0000  
    ##  1st Qu.:    0.069   1st Qu.:0.7674   1st Qu.:1.1180   1st Qu.:0.2470  
    ##  Median :    0.084   Median :0.8400   Median :1.1920   Median :0.4992  
    ##  Mean   :   85.938   Mean   :0.8118   Mean   :   Inf   Mean   :0.4997  
    ##  3rd Qu.:    0.121   3rd Qu.:0.8970   3rd Qu.:1.2980   3rd Qu.:0.7514  
    ##  Max.   :24380.000   Max.   :2.8200   Max.   :   Inf   Max.   :1.0000  
    ##       beta           other_allele        variant_id       
    ##  Min.   :-20.92710   Length:1442662     Length:1442662    
    ##  1st Qu.: -0.06443   Class :character   Class :character  
    ##  Median : -0.00070   Mode  :character   Mode  :character  
    ##  Mean   : -0.06035                                        
    ##  3rd Qu.:  0.06203                                        
    ##  Max.   : 23.26092

## Save Cleaned Dataset

``` r
saveRDS(gwas, "../data/processed/harmonized/gwas_sjogrens_clean.rds")
```
