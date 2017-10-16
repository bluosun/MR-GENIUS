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
The package currently implements two MR GENIUS estimators, *genius_addY* and *genius_mulY*, for additive and multiplicative outcome models, respectively. The functions report the MR GENIUS estimate, estimated variance, confidence interval at user-specified significance level as well as the p-value corresponding to the two-sided Wald test of null causal effect of the exposure on the outcome. 

## MR GENIUS under an additive outcome model

*genius_addY* implements the estimators given in equations (6) and (12) of Tchetgen Tchetgen et al (2017) for single and multiple instruments, respectively. The term E[*A*|*G*] is modelled under the logit and identity links for binary and continuous exposure respectively, with a default linear predictor consisting of the main effects of all available instruments.  

```r
#Y      : A numeric vector of outcomes
#A      : A numeric vector of exposures (binary values should be coded in 1/0)
#G      : A numeric matrix of IVs; each column stores values for one IV (a numeric vector if only a single IV is available).
#formula: An object of class "formula" describing the linear predictor of the model for E[A|G] (default A~G, main effects of all available instruments).
#alpha  : Significance level for confidence interval (default value=0.05)
#lower  : The lower end point of the causal effect interval to be searched (default value=-10) 
#upper  : The upper end point of the causal effect interval to be searched (default value=-10) 

genius_addY(Y,A,G,formula=A~G,alpha=0.05,lower=-10,upper=10) 
```
## MR GENIUS under a multiplicative outcome model

*genius_mulY* implements MR GENIUS as the solution to the empirical version of equation (14) in Tchetgen Tchetgen et al (2017). The term E[*A*|*G*] is modelled under the logit and identity links for binary and continuous exposure respectively, with a default linear predictor consisting of the main effects of all available instruments.  

```r
#Y      : A numeric vector of outcomes
#A      : A numeric vector of exposures (binary values should be coded in 1/0)
#G      : A numeric matrix of IVs; each column stores values for one IV (a numeric vector if only a single IV is available).
#formula: An object of class "formula" describing the linear predictor of the model for E[A|G] (default A~G, main effects of all available instruments).
#alpha  : Significance level for confidence interval (default value=0.05)
#lower  : The lower end point of the causal effect interval to be searched (default value=-10) 
#upper  : The upper end point of the causal effect interval to be searched (default value=-10) 

genius_mulY(Y,A,G,formula=A~G,alpha=0.05,lower=-10,upper=10) 
```
## MR GENIUS under a multiplicative exposure model

*genius_mulA* implements the estimator given in Lemma 3 of Tchetgen Tchetgen et al (2017), under a multiplicative exposure model. By default, the log ratio term in equation (9) is modelled as a linear combination of the main effects of all available instruments.  

```r
#Y      : A numeric vector of outcomes
#A      : A numeric vector of exposures (binary values should be coded in 1/0)
#G      : A numeric matrix of IVs; each column stores values for one IV (a numeric vector if only a single IV is available).
#alpha  : Significance level for confidence interval (default value=0.05)
#lower  : The lower end point of the causal effect interval to be searched (default value=-10) 
#upper  : The upper end point of the causal effect interval to be searched (default value=-10) 

genius_mulA(Y,A,G,alpha=0.05,lower=-10,upper=10) 
```
# Examples

```r
# the following packages are needed to simulate data; they are not required to run "genius" package
library("msm")
library("MASS")

expit <- function(x) {
    exp(x)/(1+exp(x))
}

#################################################
######## Additive outcome model examples ########
#################################################

### example with binary exposure, all instruments invalid ###
# true causal effect, beta = 1.0
# Number of instruments, nIV = 10
# Y: vector of outcomes
# A: vector of exposures
# G: matrix of instruments, one column per instrument

nIV=10; N=5000; beta=1.5;
phi=rep(-0.02,nIV); gamma=rep(-0.15,nIV); alpha=rep(-0.5,nIV);
Gn = mvrnorm(N,rep(0,nIV),diag(rep(1,nIV)))
G  = (Gn>0)*1;
U= as.vector(phi%*%t(G))+ rtnorm(n=N,mean=0.35,lower=0.2,upper=0.5);
A = rbinom(N,1,expit(as.vector(gamma%*%t(G)))+U-0.35-as.vector(phi%*%t(G)));
Y = as.vector(alpha%*%t(G)) + beta*A + U + rnorm(N);

genius_addY(Y,A,G);

### specify a more richly parameterized linear predictor for the model 
### of E[A|G] containing all main effects and pairwise interactions of 
### instruments

colnames(G)=paste("g",1:10,sep="")

genius_addY(Y,A,G,A~(g1+g2+g3+g4+g5+g6+g7+g8+g9+g10)^2);

### example with continous exposure, all instruments invalid ###
# true causal effect, beta = 1.5
# Number of instruments, nIV = 10
# Y: vector of outcomes
# A: vector of exposures
# G: matrix of instruments, one column per instrument

nIV=10; N=500; beta=1.5;
phi=rep(-0.5,nIV); gamma=rep(-2,nIV); alpha=rep(-0.5,nIV);
lambda0=1; lambda1=rep(0.5,nIV);
Gn = mvrnorm(N,rep(0,nIV),diag(rep(1,nIV)))
G  = (Gn>0)*1;
U = as.vector(phi%*%t(G))+rnorm(N);
A = as.vector(gamma%*%t(G)) +U + rnorm(N,mean=0,sd=abs(lambda0+as.vector(lambda1%*%t(G))));
Y = as.vector(alpha%*%t(G)) + beta*A + U + rnorm(N);

genius_addY(Y,A,G);

#######################################################
######## Multiplicative outcome model examples ########
#######################################################


### example with binary exposure, all instruments invalid ###
# true causal effect, beta = 1.5
# Number of instruments, nIV = 10
# Y: vector of outcomes
# A: vector of exposures
# G: matrix of instruments, one column per instrument

nIV=10; N=2000; beta=1.5;
phi=rep(-0.02,nIV); gamma=rep(-0.15,nIV); alpha=rep(-0.5,nIV);
Gn = mvrnorm(N,rep(0,nIV),diag(rep(1,nIV)))
G  = (Gn>0)*1;
U= as.vector(phi%*%t(G))+ rtnorm(n=N,mean=0.35,lower=0.2,upper=0.5);
A = rbinom(N,1,expit(as.vector(gamma%*%t(G)))+U-0.35-as.vector(phi%*%t(G)));
Y = exp(beta*A)*(as.vector(alpha%*%t(G)) + U) + rnorm(N);

genius_mulY(Y,A,G);

### specify a more richly parameterized linear predictor for the model of E[A|G] 
### containing all main effects and pairwise interactions of instruments                                                       

colnames(G)=paste("g",1:10,sep="")

genius_mulY(Y,A,G,A~(g1+g2+g3+g4+g5+g6+g7+g8+g9+g10)^2);

### continuous exposure
nIV=10; N=2000; beta=0.25; 
phi=rep(0.2,nIV); gamma=rep(0.5,nIV); alpha=rep(0.5,nIV);

lambda0=0.5; lambda1=rep(0.5,nIV);
Gn = mvrnorm(N,rep(0,nIV),diag(rep(1,nIV)))
G  = (Gn>0)*1;
U = as.vector(phi%*%t(G))+rnorm(N);
A = as.vector(gamma%*%t(G)) +U + rnorm(N,mean=0,sd=abs(lambda0+as.vector(lambda1%*%t(G))));
Y = exp(beta*A)*(as.vector(alpha%*%t(G)) + U) + rnorm(N);

genius_mulY(Y,A,G);

########################################################
######## Multiplicative exposure model examples ########
########################################################

nIV=10; N=2000; beta=0.5; 
gamma=rep(0.5,nIV); alpha=rep(0.5,nIV);phi=rep(0.05,nIV);
Gn = mvrnorm(N,rep(0,nIV),diag(rep(1,nIV)))
G  = (Gn>0)*1;
U = as.vector(phi%*%t(G))+rnorm(N);
#exposure generated from negative binomial distribution
A = rnbinom(N,size=10,mu = exp(as.vector(gamma%*%t(G)) +0.1*U)) 
Y = as.vector(alpha%*%t(G)) + beta*A + U + rnorm(N);

genius.mulA(Y,A,G);
```

# References 
Tchetgen Tchetgen, E., Sun, B. and Walter, S. (2017). <a href="https://arxiv.org/abs/1709.07779"> The GENIUS Approach to Robust Mendelian Randomization Inference.</a> arXiv e-prints.


