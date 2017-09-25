#' G-Estimation under No-Interaction with Unmeasured Selection
#'
#' Implements G-Estimation under No-Interaction with Unmeasured Selection.
#'
#' @param Y A numeric vector of outcomes.
#' @param A A numeric vector of exposures (binary values should be coded in 0/1).
#' @param G A numeric matrix of instruments; each column stores values for one instrument.
#' @param alpha Significance level for confidence interval (default value=0.05).
#' @param lower The lower end point of the causal effect interval to be searched (default value=-10).
#' @param upper The upper end point of the causal effect interval to be searched (default value=10).
#'
#' @return A "genius" object containing the following items:
#' \item{beta.est}{The point estimate of causal effect of the exposure on the outcome.}
#' \item{beta.var}{The corresponding estimated variance.}
#' \item{ci}{The corresponding Wald-type confidence interval at specified significance level.}
#' @export
genius <- function(Y,A,G,alpha=0.05,lower=-10,upper=10) {

	#sample size
	N = dim(G)[1];
	#number of IVs
	nIV = dim(G)[2];
	A.binary = all(A %in% 0:1);

	if (A.binary) {
      	glm.AG = glm(A~G, family=binomial, x=T);
	} else {
		glm.AG = lm(A~G, x=T);
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
	ci = beta.est + c(-1,1)*qnorm(1-alpha/2)*sqrt(beta.var);
	#resid.glm = glm(residuals(glm.AG)^2~G, family = Gamma(link = log));

 	object <- list(beta.est=beta.est, beta.var = beta.var, ci=ci)
  	class(object) <- "genius"
  	return(object)
}

expit <- function(x) {
    exp(x)/(1+exp(x))
}

mm.binary <- function (b, n_iv, design.mat, exposure, iv, outcome, n_sample) {
	#G-/mu
	mu.e = sweep(iv,2,apply(iv,2,mean),"-");
	#normal equations
	nr.e = design.mat*as.vector(exposure-expit(design.mat%*%b[1:(n_iv+1)]));
	#IV equations
	iv.e = sweep(iv,2,apply(iv,2,mean),"-")*as.vector(exposure-expit(design.mat%*%b[1:(n_iv+1)]))*(outcome-b[n_iv+2]*exposure);
	#centering the IV moment conditions
	iv.e = sweep(iv.e,2, apply(iv.e,2,mean),"-");
	tildem = cbind(mu.e, nr.e, iv.e);
	(1/n_sample)*(t(tildem)%*%tildem);
}

mm.linear <- function (b, n_iv, design.mat, exposure, iv, outcome, n_sample) {
	#G-/mu
	mu.e = sweep(iv,2,apply(iv,2,mean),"-");
	#normal equations
	nr.e = design.mat*as.vector(exposure-(design.mat%*%b[1:(n_iv+1)]));
	#IV equations
	iv.e = sweep(iv,2,apply(iv,2,mean),"-")*as.vector(exposure-(design.mat%*%b[1:(n_iv+1)]))*(outcome-b[n_iv+2]*exposure);
	#centering the IV moment conditions
	iv.e = sweep(iv.e,2, apply(iv.e,2,mean),"-");
	tildem =cbind(mu.e, nr.e, iv.e);
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
	A_prime = dudbeta.binary(b, n_iv, design.mat, exposure, iv, outcome, n_sample);
	Omega_inv = W.opt.binary(b, n_iv, design.mat, exposure, iv, outcome, n_sample); 
	dudmu = -diag(n_iv)*mean(as.vector(exposure-expit(design.mat%*%b[1:(n_iv+1)]))*(outcome-b[n_iv+2]*exposure));
	dudpsi = -(1/n_sample)*t(sweep(iv,2,apply(iv,2,mean),"-")*(outcome-b[n_iv+2]*exposure))%*%(pvar*design.mat);
	B[(2*n_iv+2),] = (A_prime%*%Omega_inv)%*%cbind(dudmu, dudpsi,t(A_prime));
	B
}

B.linear <- function (b, n_iv, design.mat, exposure, iv, outcome, n_sample) {
	B = matrix(0,(2*n_iv +2),(2*n_iv+2));
	B[1:n_iv,1:n_iv] = -diag(n_iv);
	B[(n_iv+1):(2*n_iv+1),(n_iv+1):(2*n_iv+1)] = -(1/n_sample)*(t(design.mat)%*%(design.mat));
	A_prime = dudbeta.linear(b, n_iv, design.mat, exposure, iv, outcome, n_sample);
	Omega_inv = W.opt.linear(b, n_iv, design.mat, exposure, iv, outcome, n_sample); 
	dudmu = -diag(n_iv)*mean(as.vector(exposure-(design.mat%*%b[1:(n_iv+1)]))*(outcome-b[n_iv+2]*exposure));
	dudpsi = -(1/n_sample)*t(sweep(iv,2,apply(iv,2,mean),"-")*(outcome-b[n_iv+2]*exposure))%*%design.mat;
	B[(2*n_iv+2),] = (A_prime%*%Omega_inv)%*%cbind(dudmu, dudpsi,t(A_prime));
	B
}