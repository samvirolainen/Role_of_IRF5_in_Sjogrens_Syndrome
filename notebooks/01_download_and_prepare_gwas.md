01 — Download & Prepare GWAS Data
================
Sam Virolainen
2025-11-22

\#Set up Working Directory

## Set up Directory Structure

``` r
dir.create("data/raw/gwas", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed/harmonized", recursive = TRUE, showWarnings = FALSE)
```

## Load the GWAS Dataset

``` r
gwas_path <- "../data/raw/gwas/GCST012796_buildGRCh37.csv"

if (!file.exists(gwas_path)) {
  stop("GWAS file not found. Make sure it is placed in data/raw/gwas/")
}

gwas <- fread(gwas_path, sep = ",")
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
    ##    ci_lower ci_upper p_value
    ##       <num>    <num>   <num>
    ## 1:   0.8017    1.144  0.6326
    ## 2:   0.8097    1.192  0.8581
    ## 3:   0.7350    1.364  0.9941
    ## 4:   0.6774    1.326  0.7541
    ## 5:   0.7829    1.221  0.8409
    ## 6:   0.1887    6.905  0.8854

## Convert OR to Beta

``` r
gwas$beta <- log(gwas$odds_ratio)
```

## Inspect the Data Structure

``` r
glimpse(gwas)
```

    ## Rows: 1,442,662
    ## Columns: 9
    ## $ chromosome         <int> 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, …
    ## $ base_pair_location <int> 752566, 752721, 754063, 759036, 771967, 774047, 779…
    ## $ effect_allele      <chr> "C", "T", "T", "A", "A", "A", "G", "A", "G", "A", "…
    ## $ odds_ratio         <dbl> 0.9576, 0.9825, 1.0010, 0.9477, 0.9775, 1.1410, 0.9…
    ## $ standard_error     <dbl> 0.09064, 0.09870, 0.15770, 0.17140, 0.11330, 0.9184…
    ## $ ci_lower           <dbl> 0.8017, 0.8097, 0.7350, 0.6774, 0.7829, 0.1887, 0.8…
    ## $ ci_upper           <dbl> 1.144, 1.192, 1.364, 1.326, 1.221, 6.905, 1.170, 1.…
    ## $ p_value            <dbl> 0.6326, 0.8581, 0.9941, 0.7541, 0.8409, 0.8854, 0.7…
    ## $ beta               <dbl> -0.0433251247, -0.0176549352, 0.0009995003, -0.0537…

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
    ##       beta          
    ##  Min.   :-20.92710  
    ##  1st Qu.: -0.06443  
    ##  Median : -0.00070  
    ##  Mean   : -0.06035  
    ##  3rd Qu.:  0.06203  
    ##  Max.   : 23.26092

## Verify Required Columns for Downstream Analysis

``` r
required_cols <- c(
  "chromosome",
  "base_pair_location",
  "effect_allele",
  "odds_ratio",
  "standard_error",
  "p_value"
)

required_cols %in% names(gwas)
```

    ## [1] TRUE TRUE TRUE TRUE TRUE TRUE

## Check for Missing Values

``` r
colSums(is.na(gwas[, ..required_cols]))
```

    ##         chromosome base_pair_location      effect_allele         odds_ratio 
    ##                  0                  0                  0                  0 
    ##     standard_error            p_value 
    ##                  0                  0

## Save the raw dataset for use in Notebook 2

``` r
saveRDS(gwas, "../data/processed/harmonized/gwas_sjogrens_raw.rds")
```

## Notebook 02 will:

# Align alleles between GWAS and eQTL

# Ensure consistent genome build

# Filter down to overlapping variants

# Prepare data structures for coloc

# Select loci or genes for testing
