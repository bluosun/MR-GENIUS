# MR-GENIUS

The "genius" R package implements mendelian randomization G-estimation under no interaction with unmeasured selection 
(MR GENIUS), based on work by Tchetgen Tchetgen, E., Sun, B. and Walter, S. (2017).

# How to Install

To install this package, use devtools:

```r
devtools::install_github("bluosun/MR-GENIUS")
```
(requires R version >= 3.4.1) 
# Overview
The package currently offers one function "genius" which returns the MR GENIUS estimate, estimated variance as well as corresponding confidence interval.

```r
#Y: A numeric vector of outcomes
#A: A numeric vector of exposures (binary values should be coded in 1/0)
#G: A numeric matrix of IVs; each column stores values for one IV
#alpha: Significance level for confidence interval (default value=0.05)
#lower: The lower end point of the causal effect interval to be searched (default value=-10) 
#upper: The upper end point of the causal effect interval to be searched (default value=-10) 

genius(Y, A, G, alpha = 0.05, lower=-10, upper=10)
```

## References 
Tchetgen Tchetgen, E., Sun, B. and Walter, S. (2017). <a href="https://arxiv.org/abs/1709.07779"> The GENIUS Approach to Robust Mendelian Randomization Inference.</a> arXiv e-prints.


