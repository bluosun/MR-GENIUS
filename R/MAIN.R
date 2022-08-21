#' @title MR GENIUS under additive outcome model
#' 
#' @description
#' Implements MR GENIUS under an additive outcome model.
#'
#' @details
#' This function implements the estimators given in equations (6) and (12) of Tchetgen Tchetgen et al (2017) for single and multiple instruments, respectively. The term
#' \eqn{E(A|G)} is modelled under the logit and identity links for binary and continuous exposure respectively, with a default linear predictor consisting of the main effects of all available instruments.  
#'
#' @references
#' Tchetgen Tchetgen, E., Sun, B. and Walter, S. (2017). \href{https://arxiv.org/abs/1709.07779}{The GENIUS Approach to Robust Mendelian Randomization Inference}. arXiv e-prints.
#' 
#' @param Y A numeric vector of outcomes.
#' @param A A numeric vector of exposures (binary values should be coded in 0/1).
#' @param G A numeric matrix of instruments; each column stores values for one instrument (a numeric vector if only a single instrument is available).
#' @param formula An object of class "formula" describing the linear predictor of the model for \eqn{E(A|G)} (default \eqn{A~G}, main effects of all available instruments).
#' @param alpha Significance level for confidence interval (default value=0.05).
#' @param lower The lower end point of the causal effect interval to be searched (default value=-10).
#' @param upper The upper end point of the causal effect interval to be searched (default value=10).
#'
#' @return A "genius" object containing the following items:
#' \item{beta.est}{The point estimate of the causal effect (on the additive scale) of the exposure on the outcome.}
#' \item{beta.var}{The corresponding estimated variance.}
#' \item{ci}{The corresponding Wald-type confidence interval at specified significance level.}
#' \item{pval}{The p-value for two-sided Wald test of null causal effect (on the additive scale) of the exposure on the outcome.}
#'
#' @examples
#'# the following packages are needed to simulate data
#'library("msm")
#'library("MASS")
#'expit <- function(x) {
#'    exp(x)/(1+exp(x))
#'}
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
#'genius_addY(Y,A,G);
#'
#'### specify a more richly parameterized linear predictor for the model 
#'### of E[A|G] containing all main effects and pairwise interactions of 
#'### instruments
#'
#'colnames(G)=paste("g",1:10,sep="")
#'
#'genius_addY(Y,A,G,A~(g1+g2+g3+g4+g5+g6+g7+g8+g9+g10)^2);
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
#'genius_addY(Y,A,G);
#'
#'@export

genius_addY <- function(Y,A,G,formula=A~G,alpha=0.05,lower=-10,upper=10) {

	if (is.data.frame(Y)) {
		Y=as.vector(data.matrix(Y));
	}
	if (is.data.frame(A)) {
		A=as.vector(data.matrix(A));
	}
	if (is.data.frame(G)) {
		G=data.matrix(G)
	}
	if (inherits(G, "matrix")) {
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
      		glm.AG = stats::glm(formula(formula), family=binomial, x=T, data=cbind(A,as.data.frame(G)));
		} else {
			glm.AG = stats::lm(formula(formula), x=T, data=cbind(A,as.data.frame(G)));
		}
		#Multiple-IV GMM
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

		beta.var= diag((1/N)*(solve(B)%*%M)%*%mm%*%t((solve(B)%*%M)))[nIV+length(glm.AG$coef)+1];
	}
	else {
		if (A.binary) {
      		glm.AG = stats::glm(formula(formula), family=binomial, x=T, data=cbind(A,as.data.frame(G)));
		} else {
			glm.AG = stats::lm(formula(formula), x=T, data=cbind(A,as.data.frame(G)));
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
		beta.var=diag((1/N)*(solve(B))%*%mm%*%t((solve(B))))[nIV+length(glm.AG$coef)+1];
	}
	ci = beta.est + c(-1,1)*stats::qnorm(1-alpha/2)*sqrt(beta.var);
	pval = 2*pnorm(-abs(beta.est/sqrt(beta.var)));
 	object <- list(beta.est=unname(beta.est), beta.var = beta.var, ci=ci, pval=unname(pval))
  	class(object) <- "genius"
  	return(object)
}






#' @title MR GENIUS under multiplicative outcome model
#' 
#' @description
#' Implements MR GENIUS under a multiplicative outcome model.
#'
#' @details
#' This function implements MR GENIUS as the solution to the empirical version of equation (14) in Tchetgen Tchetgen et al (2017). The term
#' \eqn{E(A|G)} is modelled under the logit and identity links for binary and continuous exposure respectively, with a default linear predictor consisting of the main effects of all available instruments.  
#'
#' @references
#' Tchetgen Tchetgen, E., Sun, B. and Walter, S. (2017). \href{https://arxiv.org/abs/1709.07779}{The GENIUS Approach to Robust Mendelian Randomization Inference}. arXiv e-prints.
#' 
#' @param Y A numeric vector of outcomes.
#' @param A A numeric vector of exposures (binary values should be coded in 0/1).
#' @param G A numeric matrix of instruments; each column stores values for one instrument (a numeric vector if only a single instrument is available).
#' @param formula An object of class "formula" describing the linear predictor of the model for \eqn{E(A|G)} (default \eqn{A~G}, main effects of all available instruments).
#' @param alpha Significance level for confidence interval (default value=0.05).
#' @param lower The lower end point of the causal effect interval to be searched (default value=-10).
#' @param upper The upper end point of the causal effect interval to be searched (default value=10).
#'
#' @return A "genius" object containing the following items:
#' \item{beta.est}{The point estimate of the causal effect (on the multiplicative scale) of the exposure on the outcome.}
#' \item{beta.var}{The corresponding estimated variance.}
#' \item{ci}{The corresponding Wald-type confidence interval at specified significance level.}
#' \item{pval}{The p-value for two-sided Wald test of null causal effect (on the multiplicative scale) of the exposure on the outcome.}
#'
#' @examples
#'#the following packages are needed to simulate data
#'library("msm")
#'library("MASS")
#'
#'### examples under multiplicative outcome model, all instruments invalid ###
#'# true causal effect, beta = 1.5
#'# Number of instruments, nIV = 10
#'# Y: vector of outcomes
#'# A: vector of exposures
#'# G: matrix of instruments, one column per instrument
#'
#'### binary exposure
#'nIV=10; N=2000; beta=1.5;
#'phi=rep(-0.02,nIV); gamma=rep(-0.15,nIV); alpha=rep(-0.5,nIV);
#'Gn = mvrnorm(N,rep(0,nIV),diag(rep(1,nIV)))
#'G  = (Gn>0)*1;
#'U= as.vector(phi%*%t(G))+ rtnorm(n=N,mean=0.35,lower=0.2,upper=0.5);
#'A = rbinom(N,1,expit(as.vector(gamma%*%t(G)))+U-0.35-as.vector(phi%*%t(G)));
#'Y = exp(beta*A)*(as.vector(alpha%*%t(G)) + U) + rnorm(N);
#'
#'genius_mulY(Y,A,G);
#'
#'### specify a more richly parameterized linear predictor for the model of E[A|G] 
#'### containing all main effects and pairwise interactions of instruments                                                       
#'
#'colnames(G)=paste("g",1:10,sep="")
#'
#'genius_mulY(Y,A,G,A~(g1+g2+g3+g4+g5+g6+g7+g8+g9+g10)^2);
#'
#'### continuous exposure
#'nIV=10; N=2000; beta=0.25; 
#'phi=rep(0.2,nIV); gamma=rep(0.5,nIV); alpha=rep(0.5,nIV);
#'lambda0=0.5; lambda1=rep(0.5,nIV);
#'Gn = mvrnorm(N,rep(0,nIV),diag(rep(1,nIV)))
#'G  = (Gn>0)*1;
#'U = as.vector(phi%*%t(G))+rnorm(N);
#'A = as.vector(gamma%*%t(G)) +U + rnorm(N,mean=0,sd=abs(lambda0+as.vector(lambda1%*%t(G))));
#'Y = exp(beta*A)*(as.vector(alpha%*%t(G)) + U) + rnorm(N);
#'
#'genius_mulY(Y,A,G);
#'
#'@export

genius_mulY <- function(Y,A,G,formula=A~G,alpha=0.05,lower=-10,upper=10) {

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
      		glm.AG = stats::glm(formula(formula), family=binomial, x=T, data=cbind(A,as.data.frame(G)));
		} else {
			glm.AG = stats::lm(formula(formula), x=T, data=cbind(A,as.data.frame(G)));
		}
		#Multiple-IV GMM
		gmm.fun <- function(beta, X) {
			sweep(X,2,apply(X,2,mean),"-")*(A-glm.AG$fit)*(Y*exp(-beta*A))
		}
		gmm.out = gmm::gmm(gmm.fun,x=G,t0=c(lower,upper),optfct="optimize",
			type="iterative",wmatrix="optimal",vcov="iid",centeredVcov=TRUE);

		beta.est=gmm.out$coef;

		if (A.binary) {
			mm= mm.binary.mulY(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N); 
			M = M.binary.mulY(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N);
			B = B.binary.mulY(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N);
		} else {
			mm= mm.linear.mulY(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N); 
			M = M.linear.mulY(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N);
			B = B.linear.mulY(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N);
		}

		beta.var= diag((1/N)*(solve(B)%*%M)%*%mm%*%t((solve(B)%*%M)))[nIV+length(glm.AG$coef)+1];
	}
	else {
		if (A.binary) {
      		glm.AG = stats::glm(formula(formula), family=binomial, x=T, data=cbind(A,as.data.frame(G)));
		} else {
			glm.AG = stats::lm(formula(formula), x=T, data=cbind(A,as.data.frame(G)));
		}
		#single-IV estimator 
		est.fun <- function(x) {
			mean((G-mean(G))* (A-glm.AG$fit) * (Y*exp(-x*A)));
		}

		beta.est=stats::uniroot(est.fun, c(lower,upper))$root;

		if (A.binary) {
			mm= mm.binary.mulY(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N); 
			B = B.binary.mulY(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N);
		} else {
			mm= mm.linear.mulY(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N); 
			B = B.linear.mulY(c(glm.AG$coef,beta.est), nIV, glm.AG$x, A, G, Y, N);
		}
		beta.var=diag((1/N)*(solve(B))%*%mm%*%t((solve(B))))[nIV+length(glm.AG$coef)+1];
	}
	ci = beta.est + c(-1,1)*stats::qnorm(1-alpha/2)*sqrt(beta.var);
	pval = 2*pnorm(-abs(beta.est/sqrt(beta.var)));
 	object <- list(beta.est=unname(beta.est), beta.var = beta.var, ci=ci, pval=unname(pval))
  	class(object) <- "genius"
  	return(object)
}


#' @title MR GENIUS under multiplicative exposure model
#' 
#' @description
#' Implements MR GENIUS under a multiplicative exposure model.
#'
#' @details
#' This function implements the estimator given in Lemma 3 of Tchetgen Tchetgen et al (2017), under a multiplicative exposure model. By default, the log ratio term in equation (9) is modelled as a linear combination of the main effects of all available instruments.  
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
#' \item{beta.est}{The point estimate of the causal effect (on the additive scale) of the exposure on the outcome.}
#' \item{beta.var}{The corresponding estimated variance.}
#' \item{ci}{The corresponding Wald-type confidence interval at specified significance level.}
#' \item{pval}{The p-value for two-sided Wald test of null causal effect (on the additive scale) of the exposure on the outcome.}
#'
#' @examples
#'#the following package is needed to simulate data
#'library("MASS")
#'nIV=10; N=2000; beta=0.5; 
#'gamma=rep(0.5,nIV); alpha=rep(0.5,nIV);phi=rep(0.05,nIV);
#'Gn = mvrnorm(N,rep(0,nIV),diag(rep(1,nIV)))
#'G  = (Gn>0)*1;
#'U = as.vector(phi%*%t(G))+rnorm(N);
#'#exposure generated from negative binomial distribution
#'A = rnbinom(N,size=10,mu = exp(as.vector(gamma%*%t(G)) +0.1*U)) 
#'Y = as.vector(alpha%*%t(G)) + beta*A + U + rnorm(N);
#'
#'genius_mulA(Y,A,G);
#'
#'
#' @export

genius_mulA <- function(Y,A,G,alpha=0.05,lower=-10,upper=10) {


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

	      #Estimation for pi (vector of parameters in equation 9)
		pi.func <- function(pi,A,G) {
			apply(as.vector(A*exp(-G%*%pi))*sweep(G,2,apply(G,2,mean),"-"),2,mean)
		}
		pi.fit <- BB::dfsane(par = rep(0,nIV), fn = pi.func,
    			    control = list(NM = TRUE, M = 100, noimp = 500, trace = FALSE),
    			    A=A, G = G)

		#GMM for multiple-IV, empirical version of equation in Lemma 3
		gmm.fun <- function(beta, X) {
			sweep(X,2,apply(X,2,mean),"-")*(as.vector(A*exp(-X%*%pi.fit$par))-mean((A*exp(-X%*%pi.fit$par))))*(Y-beta*A)
		}
		gmm.out = gmm::gmm(gmm.fun,x=G,t0=c(lower,upper),optfct="optimize",
			type="iterative",wmatrix="optimal",vcov="iid",centeredVcov=TRUE);

		beta.est=gmm.out$coef;

		#estimated asymptotic variance
		mm= mm.mulA(beta.est, pi.fit$par, nIV, A, G, Y, N); 
		M = M.mulA(beta.est, pi.fit$par, nIV,  A, G, Y, N);
		B = B.mulA(beta.est, pi.fit$par, nIV,  A, G, Y, N);

		beta.var= diag((1/N)*(solve(B)%*%M)%*%mm%*%t((solve(B)%*%M)))[2*nIV+2];

	} else {

		pi.func <- function(pi,A,G) {
			mean(A*exp(-G*pi)*(G-mean(G)));
		}
		pi.fit <- BB::BBoptim(par = rep(0,nIV), fn = pi.func, A=A, G = G, quiet=TRUE);
		#single-IV estimator 
		beta.est = mean((G-mean(G))*(as.vector(A*exp(-G*pi.fit$par))-mean((A*exp(-G*pi.fit$par))))*Y) /
			     mean((G-mean(G))*(as.vector(A*exp(-G*pi.fit$par))-mean((A*exp(-G*pi.fit$par))))*A);

		#estimated asymptotic variance
		mm= mm.mulA(beta.est, pi.fit$par, nIV, A, G, Y, N); 
		B = B.mulA(beta.est, pi.fit$par, nIV, A, G, Y, N);

		beta.var=diag((1/N)*(solve(B))%*%mm%*%t((solve(B))))[2*nIV+2];
	}

	ci = beta.est + c(-1,1)*stats::qnorm(1-alpha/2)*sqrt(beta.var);
	pval = 2*pnorm(-abs(beta.est/sqrt(beta.var)));
 	object <- list(beta.est=unname(beta.est), beta.var = beta.var, ci=ci, pval=unname(pval))
  	class(object) <- "genius"
  	return(object)
}

#' @title Breusch-Pagan Test
#' 
#' @description
#' Performs the Breusch-Pagan test against heteroskedasticity. This function is exported from the "lmtest" package (Achim Zeileis & Torsten Hothorn, 2002).
#' This function provides a way to test the heteroskedasticity identification condition (5) in Lemma 1 of Tchetgen Tchetgen et al. (2017).
#' @usage
#' bptest(formula, varformula = NULL, studentize = TRUE, data = list())
#' @details
#' The Breusch-Pagan test fits a linear regression model to the residuals
#' of a linear regression model (by default the same explanatory variables are taken as in the main regression
#' model) and rejects if too much of the variance is explained by the additional explanatory variables. Under H0 the test statistic of the Breusch-Pagan test follows a
#' chi-squared distribution with \code{parameter} (the number of regressors without the constant in the model) degrees of freedom. 
#'
#' @param formula a symbolic description for the model to be tested (or a fitted \code{"lm"} object).
#' @param varformula a formula describing only the potential explanatory variables
#'   for the variance (no dependent variable needed). By default the same
#'   explanatory variables are taken as in the main regression model.
#' @param studentize logical. If set to \code{TRUE} Koenker's studentized version of the test statistic will be used.
#' @param data an optional data frame containing the variables in the model. By default the variables are taken from the environment which \code{bptest} is called from.
#'
#' @return A list with class \code{"htest"} containing the following components:
#'  \item{statistic}{the value of the test statistic.}
#'  \item{p.value}{the p-value of the test.}
#'  \item{parameter}{degrees of freedom.}
#'  \item{method}{a character string indicating what type of test was performed.}
#'  \item{data.name}{a character string giving the name(s) of the data.}
#'
#' @examples
#'# the following packages are needed to simulate data
#'library("msm")
#'library("MASS")
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
#'bptest(A~G);
#
#' @references
#' T.S. Breusch & A.R. Pagan (1979),
#' A Simple Test for Heteroscedasticity and Random Coefficient Variation.
#' \emph{Econometrica} \bold{47}, 1287--1294
#'
#' R. Koenker (1981), A Note on Studentizing a Test for Heteroscedasticity.
#' \emph{Journal of Econometrics} \bold{17}, 107--112.
#'
#' W. Kramer & H. Sonnberger (1986),
#' \emph{The Linear Regression Model under Test}. Heidelberg: Physica
#'
#' Achim Zeileis & Torsten Hothorn (2002), Diagnostic Checking in Regression Relationships. \emph{R News} \bold{2(3)}, 7-10. \url{https://CRAN.R-project.org/doc/Rnews/}
#'
#' Tchetgen Tchetgen, E., Sun, B. and Walter, S. (2017). \href{https://arxiv.org/abs/1709.07779}{The GENIUS Approach to Robust Mendelian Randomization Inference}. arXiv e-prints.
#'
#' @importFrom lmtest bptest
#' @name bptest
#' @export
NULL

