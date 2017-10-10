#' @title G-Estimation under No-Interaction with Unmeasured Selection
#' 
#' @description
#' Implements G-Estimation under No-Interaction with Unmeasured Selection.
#'
#' @details
#' This function implements estimation of causal effect under an additive outcome model. The estimator is given in equations (6) and (12) of Tchetgen Tchetgen et al (2017) for single and multiple instruments, respectively. The term
#' \eqn{E(A|G)} is modelled under the logit and identity links for binary and continuous exposure respectively, with a linear predictor consisting of the main effects of all available instruments.  
#'
#' @references
#' Tchetgen Tchetgen, E., Sun, B. and Walter, S. (2017). \href{https://arxiv.org/abs/1709.07779}{The GENIUS Approach to Robust Mendelian Randomization Inference}. arXiv e-prints.
#' 
#' @param Y A numeric vector of outcomes.
#' @param A A numeric vector of exposures (binary values should be coded in 0/1).
#' @param G A numeric matrix of instruments; each column stores values for one instrument (a numeric vector if only a single instrument is available).
#' @param alpha Significance level for confidence interval (default value=0.05).
#' @param lower The lower end point of the causal effect interval to be searched (default value=-10).
#' @param upper The upper end point of the causal effect interval to be searched (default value=10).
#'
#' @return A "genius" object containing the following items:
#' \item{beta.est}{The point estimate of causal effect of the exposure on the outcome.}
#' \item{beta.var}{The corresponding estimated variance.}
#' \item{ci}{The corresponding Wald-type confidence interval at specified significance level.}
#' \item{pval}{The p-value for two-sided Wald test of null causal effect of the exposure on the outcome.}
#'
#' @examples
#'# the following packages are needed to simulate data
#'library("msm")
#'library("MASS")
#'
#'### example with binary exposure, all instruments invalid ###
#'# true causal effect, beta = 1.0
#'# Number of instruments, nIV = 10
#'# Y: vector of outcomes
#'# A: vector of exposures
#'# G: matrix of instruments, one column per instrument
#'
#'nIV=10; N=5000; beta=1;
#'phi=rep(-0.02,nIV); gamma=rep(-0.15,nIV); alpha=rep(-0.5,nIV);
#'Gn = mvrnorm(N,rep(0,nIV),diag(rep(1,nIV)))
#'
#'G  = (Gn>0)*1;
#'U= as.vector(phi%*%t(G))+ rtnorm(n=N,mean=0.35,lower=0.2,upper=0.5);
#'A = rbinom(N,1,expit(as.vector(gamma%*%t(G)))+U-0.35-as.vector(phi%*%t(G)));
#'Y = as.vector(alpha%*%t(G)) + beta*A + U + rnorm(N);
#'
#'genius(Y,A,G);
#'
#'### example with continous exposure, all instruments invalid ###
#'
#'nIV=10; N=500; beta=1;
#'phi=rep(-0.5,nIV); gamma=rep(-2,nIV); alpha=rep(-0.5,nIV);
#'lambda0=1; lambda1=rep(0.5,nIV);
#'Gn = mvrnorm(N,rep(0,nIV),diag(rep(1,nIV)))
#'
#'G  = (Gn>0)*1;
#'U = as.vector(phi%*%t(G))+rnorm(N);
#'A = as.vector(gamma%*%t(G)) +U + rnorm(N,mean=0,sd=abs(lambda0+as.vector(lambda1%*%t(G))));
#'Y = as.vector(alpha%*%t(G)) + beta*A + U + rnorm(N);
#'
#'genius(Y,A,G);
#'@export
genius <- function(Y,A,G,alpha=0.05,lower=-10,upper=10) {

	if (is.data.frame(Y)) {
		Y=as.vector(data.matrix(Y));
	}
	if (is.data.frame(A)) {
		A=as.vector(data.matrix(A));
	}
	if (is.data.frame(G)) {
		G=data.matrix(G)
	}
	if (class(G) == "matrix") {
		#number of IVs
  		nIV =dim(G)[2];
		#sample size
		N = dim(G)[1];
		if (nIV==1) {G=as.vector(G);}
	} else {
  		nIV =1;
		N = length(G);
	}
	A.binary = all(A %in% 0:1);

	if (nIV>1) {
		if (A.binary) {
      		glm.AG = stats::glm(A~G, family=binomial, x=T);
		} else {
			glm.AG = stats::lm(A~G, x=T);
		}

		#GMM 
		gmm.fun <- function(beta, X) {
			sweep(X,2,apply(X,2,mean),"-")*(A-glm.AG$fit)*(Y-beta*A)
		}

		gmm.out = gmm::gmm(gmm.fun,x=G,t0=c(lower,upper),optfct="optimize",
			type="iterative",wmatrix="optimal",vcov="iid",centeredVcov=TRUE);

		beta.est=gmm.out$coef;

		if (A.binary) {
			mm= mm.binary(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N); 
			M = M.binary(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N);
			B = B.binary(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N);
		} else {
			mm= mm.linear(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N); 
			M = M.linear(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N);
			B = B.linear(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N);
		}

		beta.var= diag((1/N)*(solve(B)%*%M)%*%mm%*%t((solve(B)%*%M)))[2*nIV+2];
	}
	else {
		if (A.binary) {
      		glm.AG = stats::glm(A~G, family=binomial, x=T);
		} else {
			glm.AG = stats::lm(A~G, x=T);
		}
		#single-IV estimator 
		beta.est=mean((G-mean(G))* (A-glm.AG$fit) * Y)/mean((G-mean(G))* (A-glm.AG$fit) * A);

		if (A.binary) {
			mm= mm.binary(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N); 
			B = B.binary(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N);
		} else {
			mm= mm.linear(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N); 
			B = B.linear(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N);
		}
		beta.var=diag((1/N)*(solve(B))%*%mm%*%t((solve(B))))[2*nIV+2];
	}
	ci = beta.est + c(-1,1)*stats::qnorm(1-alpha/2)*sqrt(beta.var);
	pval = 2*pnorm(-abs(beta.est/sqrt(beta.var)));
 	object <- list(beta.est=unname(beta.est), beta.var = beta.var, ci=ci, pval=unname(pval))
  	class(object) <- "genius"
  	return(object)
}

expit <- function(x) {
    exp(x)/(1+exp(x))
}

mm.binary <- function (b, n_iv, design.mat, exposure, iv, outcome, n_sample) {
	
	if (n_iv>1) {
		#G-/mu
		mu.e = sweep(iv,2,apply(iv,2,mean),"-");
		#normal equations
		nr.e = design.mat*as.vector(exposure-expit(design.mat%*%b[1:(n_iv+1)]));
		#IV equations
		iv.e = sweep(iv,2,apply(iv,2,mean),"-")*as.vector(exposure-expit(design.mat%*%b[1:(n_iv+1)]))*(outcome-b[n_iv+2]*exposure);
		#centering the IV moment conditions
		iv.e = sweep(iv.e,2, apply(iv.e,2,mean),"-");
		tildem = cbind(mu.e, nr.e, iv.e);
	}
	else {
		#G-/mu
		mu.e = iv-mean(iv);
		#normal equations
		nr.e = design.mat*as.vector(exposure-expit(design.mat%*%b[1:(n_iv+1)]));
		#IV equations
		iv.e = (iv-mean(iv))*as.vector(exposure-expit(design.mat%*%b[1:(n_iv+1)]))*(outcome-b[n_iv+2]*exposure);
		#centering the IV moment conditions
		iv.e = iv.e-mean(iv.e);
		tildem = cbind(mu.e, nr.e, iv.e);
	}
	(1/n_sample)*(t(tildem)%*%tildem);
}

mm.linear <- function (b, n_iv, design.mat, exposure, iv, outcome, n_sample) {

	if (n_iv>1) {
		#G-/mu
		mu.e = sweep(iv,2,apply(iv,2,mean),"-");
		#normal equations
		nr.e = design.mat*as.vector(exposure-(design.mat%*%b[1:(n_iv+1)]));
		#IV equations
		iv.e = sweep(iv,2,apply(iv,2,mean),"-")*as.vector(exposure-(design.mat%*%b[1:(n_iv+1)]))*(outcome-b[n_iv+2]*exposure);
		#centering the IV moment conditions
		iv.e = sweep(iv.e,2, apply(iv.e,2,mean),"-");
		tildem =cbind(mu.e, nr.e, iv.e);
	}
	else {
		#G-/mu
		mu.e = iv-mean(iv);
		#normal equations
		nr.e = design.mat*as.vector(exposure-(design.mat%*%b[1:(n_iv+1)]));
		#IV equations
		iv.e = (iv-mean(iv))*as.vector(exposure-(design.mat%*%b[1:(n_iv+1)]))*(outcome-b[n_iv+2]*exposure);
		#centering the IV moment conditions
		iv.e = iv.e-mean(iv.e);
		tildem =cbind(mu.e, nr.e, iv.e);
	}
	(1/n_sample)*(t(tildem)%*%tildem);
}

W.opt.linear <- function (b, n_iv, design.mat, exposure, iv, outcome, n_sample) {
	iv.e = sweep(iv,2,apply(iv,2,mean),"-")*as.vector(exposure-(design.mat%*%b[1:(n_iv+1)]))*(outcome-b[n_iv+2]*exposure);
	(1/n_sample)*(t(iv.e)%*%iv.e);
}

W.opt.binary <- function (b, n_iv, design.mat, exposure, iv, outcome, n_sample) {
	iv.e = sweep(iv,2,apply(iv,2,mean),"-")*as.vector(exposure-expit(design.mat%*%b[1:(n_iv+1)]))*(outcome-b[n_iv+2]*exposure);
	(1/n_sample)*(t(iv.e)%*%iv.e);
}

dudbeta.linear <- function (b, n_iv, design.mat, exposure, iv, outcome, n_sample) {
	db = -sweep(iv,2,apply(iv,2,mean),"-")*as.vector(exposure-(design.mat%*%b[1:(n_iv+1)]))*exposure;
	t(apply(db,2,mean))
}

dudbeta.binary <- function (b, n_iv, design.mat, exposure, iv, outcome, n_sample) {
	db = -sweep(iv,2,apply(iv,2,mean),"-")*as.vector(exposure-expit(design.mat%*%b[1:(n_iv+1)]))*exposure;
	t(apply(db,2,mean))
}

M.binary <- function (b, n_iv, design.mat, exposure, iv, outcome, n_sample) {
	M = matrix(0,(2*n_iv +2),(3*n_iv+1));
	M[1:(2*n_iv +1),1:(2*n_iv +1)] = diag((2*n_iv +1));
	A_prime = dudbeta.binary(b, n_iv, design.mat, exposure, iv, outcome, n_sample);
	Omega_inv = W.opt.binary(b, n_iv, design.mat, exposure, iv, outcome, n_sample); 
	M[(2*n_iv +2),(2*n_iv +2):(3*n_iv+1)] = A_prime%*%Omega_inv;
	M
}

M.linear <- function (b, n_iv, design.mat, exposure, iv, outcome, n_sample) {
	M = matrix(0,(2*n_iv +2),(3*n_iv+1));
	M[1:(2*n_iv +1),1:(2*n_iv +1)] = diag((2*n_iv +1));
	A_prime = dudbeta.linear(b, n_iv, design.mat, exposure, iv, outcome, n_sample);
	Omega_inv = W.opt.linear(b, n_iv, design.mat, exposure, iv, outcome, n_sample); 
	M[(2*n_iv +2),(2*n_iv +2):(3*n_iv+1)] = A_prime%*%Omega_inv;
	M
}

B.binary <- function (b, n_iv, design.mat, exposure, iv, outcome, n_sample) {

	B = matrix(0,(2*n_iv +2),(2*n_iv+2));
	B[1:n_iv,1:n_iv] = -diag(n_iv);
	#p(1-p)
	pvar = as.vector(expit(design.mat%*%b[1:(n_iv+1)])*(1-expit(design.mat%*%b[1:(n_iv+1)])));
	B[(n_iv+1):(2*n_iv+1),(n_iv+1):(2*n_iv+1)] = -(1/n_sample)*(t(design.mat)%*%(pvar*design.mat));

	if (n_iv>1) {
		A_prime = dudbeta.binary(b, n_iv, design.mat, exposure, iv, outcome, n_sample);
		Omega_inv = W.opt.binary(b, n_iv, design.mat, exposure, iv, outcome, n_sample); 
		dudmu = -diag(n_iv)*mean(as.vector(exposure-expit(design.mat%*%b[1:(n_iv+1)]))*(outcome-b[n_iv+2]*exposure));
		dudpsi = -(1/n_sample)*t(sweep(iv,2,apply(iv,2,mean),"-")*(outcome-b[n_iv+2]*exposure))%*%(pvar*design.mat);
		B[(2*n_iv+2),] = (A_prime%*%Omega_inv)%*%cbind(dudmu, dudpsi,t(A_prime));
	}
	else {
		dudmu = -1*mean(as.vector(exposure-expit(design.mat%*%b[1:(n_iv+1)]))*(outcome-b[n_iv+2]*exposure));
		dudpsi= -1*apply( (iv-mean(iv))*(outcome-b[n_iv+2]*exposure)*(pvar*design.mat) ,2, mean);
		dudb  = -mean((iv-mean(iv))*as.vector(exposure-expit(design.mat%*%b[1:(n_iv+1)]))*exposure);
		B[(2*n_iv+2),] = c(dudmu, dudpsi, dudb);
	}
	B
}

B.linear <- function (b, n_iv, design.mat, exposure, iv, outcome, n_sample) {

	B = matrix(0,(2*n_iv +2),(2*n_iv+2));
	B[1:n_iv,1:n_iv] = -diag(n_iv);
	B[(n_iv+1):(2*n_iv+1),(n_iv+1):(2*n_iv+1)] = -(1/n_sample)*(t(design.mat)%*%(design.mat));
	if (n_iv>1) {
		A_prime = dudbeta.linear(b, n_iv, design.mat, exposure, iv, outcome, n_sample);
		Omega_inv = W.opt.linear(b, n_iv, design.mat, exposure, iv, outcome, n_sample); 
		dudmu = -diag(n_iv)*mean(as.vector(exposure-(design.mat%*%b[1:(n_iv+1)]))*(outcome-b[n_iv+2]*exposure));
		dudpsi = -(1/n_sample)*t(sweep(iv,2,apply(iv,2,mean),"-")*(outcome-b[n_iv+2]*exposure))%*%design.mat;
		B[(2*n_iv+2),] = (A_prime%*%Omega_inv)%*%cbind(dudmu, dudpsi,t(A_prime));
	}
	else {
		dudmu = -1*mean(as.vector(exposure-(design.mat%*%b[1:(n_iv+1)]))*(outcome-b[n_iv+2]*exposure));
		dudpsi= -1*apply( (iv-mean(iv))*(outcome-b[n_iv+2]*exposure)*(design.mat) ,2, mean);
		dudb  = -mean((iv-mean(iv))*as.vector(exposure-(design.mat%*%b[1:(n_iv+1)]))*exposure);
		B[(2*n_iv+2),] = c(dudmu, dudpsi, dudb);
	}
	B
}