<!-- README.md is generated from README.Rmd. Please edit that file -->

# Model-based pup production estimation (mbpp) for northern fur seals

This package implements 2 forms of model-based pup production estimation
method: (1) A fully Bayesian MCMC via JAGS and (2) An asymptotically
approximate model via TMB.

## Installation

The package can be installed using the `devtools` package:

``` r
### Download model fitting code and compile
devtools::install_github("dsjohnson/mbpp")
#> Skipping install of 'mbpp' from a github remote, the SHA1 (ba057c9c) has not changed since last install.
#>   Use `force = TRUE` to force installation
library(mbpp)
mbpp::compile_mbpp_tmb()
#> Note: Using Makevars in /Users/djohnson/.R/Makevars
#> [1] 0
```

The last line is necessary to compile the `TMB` source code. This only
needs to be done once.

## Example analysis

### Load data

``` r
library(tidyverse)

data(nfs_mark)
data(nfs_resight)
# collect into islands for parallel fitting
nfs_mark <- nfs_mark %>% group_by(icode) %>% nest() %>% 
  rename(mark_data = data)
nfs_resight <- nfs_resight %>% group_by(icode) %>% nest() %>% 
  rename(resight_data = data)

nfs_data <- left_join(nfs_mark, nfs_resight)
#> Joining, by = "icode"

head(nfs_data)
#> # A tibble: 2 x 3
#>   icode mark_data         resight_data     
#>   <chr> <list>            <list>           
#> 1 SNG   <tibble [6 × 5]>  <tibble [24 × 5]>
#> 2 SNP   <tibble [13 × 5]> <tibble [52 × 5]>
```

### Parallel processing

``` r
plan(multisession, workers=2)

## JAGS model fitting
set.seed(123)
nfs_data <- nfs_data %>% 
  mutate(
    mcmc_samp = future_map2(mark_data, resight_data,
                            ~jags_mbpp(mark_data=.x, resight_data = .y,
                                       n.iter=110000, n.burn=10000,
                                       n.chains=1
                            )
    )
  )

## Asymptotic model fitting
st <- Sys.time()
nfs_data <- nfs_data %>% 
  mutate(
    asym = future_map2(mark_data, resight_data,
                       ~asym_mbpp(
                         mark_data=.x, resight_data = .y,
                         control = list(eval.max=10000, iter.max=5000),
                       )
    )
  )
plan("default")
```

### Create parametric bootstrap sample

``` r
nfs_data <- nfs_data %>% 
  mutate(
    asym_boot = map(asym, ~boot_asym(.x,size=1000))
  )
```

## Summarize \(N\) estimation

``` r
N_data <- nfs_data %>% select(icode, mcmc_samp) %>% 
  mutate(
    N_sum = map(nfs_data$mcmc_samp, 
                ~{
                  ci <- coda::HPDinterval(coda::mcmc(.x$fitting$mcmc_samples$N))
                  cbind(.x$N, ci) 
                }
    )
  ) %>% select(-mcmc_samp) %>% unnest() %>% mutate(method="MCMC")


N_data <- nfs_data %>% select(icode, asym) %>% 
  mutate(
    N_sum = map(asym, 
                ~{
                  cbind(.x$N, lower=.x$N$N-1.96*.x$N$se_N, upper=.x$N$N+1.96*.x$N$se_N) 
                }
    )
  ) %>% select(-asym) %>% unnest() %>% mutate(method="Asymp. approx.") %>% 
  bind_rows(N_data,.)

N_data <- nfs_data %>% select(icode, asym_boot) %>% 
  mutate(
    N_sum = map(asym_boot, 
                ~{
                  ci <- coda::HPDinterval(coda::mcmc(.x$N_boot))
                  cbind(.x$N, ci) 
                }
    )
  ) %>% select(-asym_boot) %>% unnest() %>% mutate(method="Bootstrap") %>% 
  bind_rows(N_data,.)
```

### Figure of reaults

![](Nplot.png)<!-- -->

# Disclaimer

*This software package is developed and maintained by scientists at the
NOAA Fisheries Alaska Fisheries Science Center and should be considered
a fundamental research communication. The reccomendations and
conclusions presented here are those of the authors and this software
should not be construed as official communication by NMFS, NOAA, or the
U.S. Dept. of Commerce. In addition, reference to trade names does not
imply endorsement by the National Marine Fisheries Service, NOAA. While
the best efforts have been made to insure the highest quality, tools
such as this are under constant development and are subject to change.*
