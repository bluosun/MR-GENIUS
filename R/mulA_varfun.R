B.mulA <- function (beta, fit.coef, n_iv, exposure, iv, outcome, n_sample) {
	B = matrix(0,(2*n_iv+2),(2*n_iv+2));
	B[1:(n_iv+1),1:(n_iv+1)] = -diag((n_iv+1));
	if (n_iv>1) {
		#extra term differentiating wrt pi (Aexp(-wG)-mu2)
		B[n_iv+1, (n_iv+2):(2*n_iv+1)] = apply(-as.vector(exposure*exp(-iv%*%fit.coef))*iv,2,mean);
		#Aexp(-wG)(G-mu1)
		B[(n_iv+2):(2*n_iv+1),1:n_iv] = -diag(n_iv)*mean(as.vector(exposure*exp(-iv%*%fit.coef)));
		B[(n_iv+2):(2*n_iv+1),(n_iv+2):(2*n_iv+1)] = -(1/n_sample)*(t(iv)%*%(iv*as.vector(exposure*exp(-iv%*%fit.coef))));
		#(G-mu1)(Aexp(-wG)-mu2)(Y-beta*A)
		A_prime = dudbeta.mulA(beta, fit.coef, n_iv, exposure, iv, outcome, n_sample);
		Omega_inv = W.opt.mulA(beta, fit.coef, n_iv, exposure, iv, outcome, n_sample); 

		dudmu = -diag(n_iv)*mean((as.vector(exposure*exp(-iv%*%fit.coef))-mean((exposure*exp(-iv%*%fit.coef))))*(outcome-beta*exposure));
		dudmu2 = -apply(sweep(iv,2,apply(iv,2,mean),"-")*(outcome-beta*exposure),2,mean);
		dudw = -(1/n_sample)*(t(sweep(iv,2,apply(iv,2,mean),"-"))%*%(iv*as.vector(exposure*exp(-iv%*%fit.coef))*(outcome-beta*exposure)));
		B[(2*n_iv+2),] = (A_prime%*%Omega_inv)%*%cbind(dudmu, dudmu2, dudw,t(A_prime));

	} else {
		#extra term differentiating wrt pi (Aexp(-wG)-mu2)
		B[n_iv+1, (n_iv+2):(2*n_iv+1)] = mean(as.vector(exposure*exp(-iv*fit.coef)));
		#Aexp(-wG)(G-mu1)
		B[(n_iv+2):(2*n_iv+1),1:n_iv] = -1*mean(as.vector(exposure*exp(-iv*fit.coef)));
		B[(n_iv+2):(2*n_iv+1),(n_iv+2):(2*n_iv+1)] = -mean(iv*iv*as.vector(exposure*exp(-iv*fit.coef)));
		#(G-mu1)(Aexp(-wG)-mu2)(Y-beta*A)
		dudmu = -1*mean((as.vector(exposure*exp(-iv*fit.coef))-mean((exposure*exp(-iv*fit.coef))))*(outcome-beta*exposure));
		dudmu2 = -mean((iv-mean(iv))*(outcome-beta*exposure));
		dudw = -mean((iv-mean(iv))*(iv*as.vector(exposure*exp(-iv*fit.coef))*(outcome-beta*exposure)));
		dudbeta = -mean((iv-mean(iv))*(as.vector(exposure*exp(-iv*fit.coef))-mean((exposure*exp(-iv*fit.coef))))*exposure)
		B[(2*n_iv+2),] = c(dudmu,dudmu2,dudw,dudbeta);
	}
	B
}

M.mulA <- function (beta, fit.coef, n_iv,  exposure, iv, outcome, n_sample) {
	M = matrix(0,(2*n_iv+2),(3*n_iv+1));
	M[1:(2*n_iv+1),1:(2*n_iv+1)] = diag((2*n_iv+1));
	A_prime = dudbeta.mulA(beta, fit.coef, n_iv,  exposure, iv, outcome, n_sample);
	Omega_inv = W.opt.mulA(beta, fit.coef, n_iv,  exposure, iv, outcome, n_sample); 
	M[(2*n_iv+2),(2*n_iv+2):(3*n_iv+1)] = A_prime%*%Omega_inv;
	M
}

dudbeta.mulA <- function (beta, fit.coef, n_iv, exposure, iv, outcome, n_sample) {
	db = -sweep(iv,2,apply(iv,2,mean),"-")*(as.vector(exposure*exp(-iv%*%fit.coef))-mean((exposure*exp(-iv%*%fit.coef))))*exposure;
	t(apply(db,2,mean))
}

W.opt.mulA <- function (beta, fit.coef, n_iv, exposure, iv, outcome, n_sample) {
	iv.e = sweep(iv,2,apply(iv,2,mean),"-")*(as.vector(exposure*exp(-iv%*%fit.coef))-mean((exposure*exp(-iv%*%fit.coef))))*(outcome-beta*exposure);
	(1/n_sample)*(t(iv.e)%*%iv.e);
}

mm.mulA <- function (beta, fit.coef, n_iv, exposure, iv, outcome, n_sample) {
	if (n_iv>1) {
		#G-\mu
		mu.e = sweep(iv,2,apply(iv,2,mean),"-");
		#Aexp(-wG)-\mu2
		mu2.e = as.vector(exposure*exp(-iv%*%fit.coef))-mean((exposure*exp(-iv%*%fit.coef)));
		#normal equations (eq 9)
		nr.e = as.vector(exposure*exp(-iv%*%fit.coef))*sweep(iv,2,apply(iv,2,mean),"-")
		#IV equations
		iv.e = sweep(iv,2,apply(iv,2,mean),"-")*(as.vector(exposure*exp(-iv%*%fit.coef))-mean((exposure*exp(-iv%*%fit.coef))))*(outcome-beta*exposure);
		#centering the IV moment conditions
		iv.e = sweep(iv.e,2, apply(iv.e,2,mean),"-");
		tildem = cbind(mu.e, mu2.e, nr.e, iv.e);
	}
	else {
		#G-\mu
		mu.e = iv-mean(iv)
		#Aexp(-wG)-\mu2
		mu2.e = exposure*exp(-iv*fit.coef)-mean(exposure*exp(-iv*fit.coef));
		#normal equations (eq 9)
		nr.e = (exposure*exp(-iv*fit.coef))*(iv-mean(iv));
		#IV equations
		iv.e = (iv-mean(iv))*(exposure*exp(-iv*fit.coef)-mean(exposure*exp(-iv*fit.coef)))*(outcome-beta*exposure);
		#centering the IV moment conditions
		iv.e = iv.e-mean(iv.e);
		tildem = cbind(mu.e, mu2.e, nr.e, iv.e);
	}
	(1/n_sample)*(t(tildem)%*%tildem);
}

