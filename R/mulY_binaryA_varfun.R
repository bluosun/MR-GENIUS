mm.binary.mulY <- function (b, n_iv, design.mat, exposure, iv, outcome, n_sample) {
	
	dim_a = length(b)-1;
	if (n_iv>1) {
		#G-/mu
		mu.e = sweep(iv,2,apply(iv,2,mean),"-");
		#normal equations
		nr.e = design.mat*as.vector(exposure-expit(design.mat%*%b[1:dim_a]));
		#IV equations
		iv.e = sweep(iv,2,apply(iv,2,mean),"-")*as.vector(exposure-expit(design.mat%*%b[1:dim_a]))*(outcome*exp(-b[dim_a+1]*exposure));
		#centering the IV moment conditions
		iv.e = sweep(iv.e,2, apply(iv.e,2,mean),"-");
		tildem = cbind(mu.e, nr.e, iv.e);
	}
	else {
		#G-/mu
		mu.e = iv-mean(iv);
		#normal equations
		nr.e = design.mat*as.vector(exposure-expit(design.mat%*%b[1:dim_a]));
		#IV equations
		iv.e = (iv-mean(iv))*as.vector(exposure-expit(design.mat%*%b[1:dim_a]))*(outcome*exp(-b[dim_a+1]*exposure));
		#centering the IV moment conditions
		iv.e = iv.e-mean(iv.e);
		tildem = cbind(mu.e, nr.e, iv.e);
	}
	(1/n_sample)*(t(tildem)%*%tildem);
}


W.opt.binary.mulY <- function (b, n_iv, design.mat, exposure, iv, outcome, n_sample) {
	dim_a=length(b)-1;
	iv.e = sweep(iv,2,apply(iv,2,mean),"-")*as.vector(exposure-expit(design.mat%*%b[1:dim_a]))*(outcome*exp(-b[dim_a+1]*exposure));
	(1/n_sample)*(t(iv.e)%*%iv.e);
}

dudbeta.binary.mulY <- function (b, n_iv, design.mat, exposure, iv, outcome, n_sample) {
	dim_a=length(b)-1;
	db = -sweep(iv,2,apply(iv,2,mean),"-")*as.vector(exposure-expit(design.mat%*%b[1:dim_a]))*exposure*outcome*exp(-b[dim_a+1]*exposure);
	t(apply(db,2,mean))
}

M.binary.mulY <- function (b, n_iv, design.mat, exposure, iv, outcome, n_sample) {
	dim_a=length(b)-1;
	M = matrix(0,(n_iv+dim_a+1),(dim_a+2*n_iv));
	M[1:(n_iv+dim_a),1:(n_iv+dim_a)] = diag((n_iv+dim_a));
	A_prime = dudbeta.binary.mulY(b, n_iv, design.mat, exposure, iv, outcome, n_sample);
	Omega_inv = W.opt.binary.mulY(b, n_iv, design.mat, exposure, iv, outcome, n_sample); 
	M[(n_iv+dim_a+1),(n_iv+dim_a+1):(dim_a+2*n_iv)] = A_prime%*%Omega_inv;
	M
}

B.binary.mulY <- function (b, n_iv, design.mat, exposure, iv, outcome, n_sample) {

	dim_a = length(b)-1;
	B = matrix(0,(dim_a+n_iv+1),(dim_a+n_iv+1));
	B[1:n_iv,1:n_iv] = -diag(n_iv);
	#p(1-p)
	pvar = as.vector(expit(design.mat%*%b[1:dim_a])*(1-expit(design.mat%*%b[1:dim_a])));
	B[(n_iv+1):(dim_a+n_iv),(n_iv+1):(dim_a+n_iv)]= -(1/n_sample)*(t(design.mat)%*%(pvar*design.mat));

	if (n_iv>1) {
		A_prime = dudbeta.binary.mulY(b, n_iv, design.mat, exposure, iv, outcome, n_sample);
		Omega_inv = W.opt.binary.mulY(b, n_iv, design.mat, exposure, iv, outcome, n_sample); 
		dudmu = -diag(n_iv)*mean(as.vector(exposure-expit(design.mat%*%b[1:dim_a]))*(outcome*exp(-b[dim_a+1]*exposure)));
		dudpsi = -(1/n_sample)*t(sweep(iv,2,apply(iv,2,mean),"-")*(outcome*exp(-b[dim_a+1]*exposure)))%*%(pvar*design.mat);
		B[(dim_a+n_iv+1),] = (A_prime%*%Omega_inv)%*%cbind(dudmu, dudpsi,t(A_prime));
	}
	else {
		dudmu = -1*mean(as.vector(exposure-expit(design.mat%*%b[1:dim_a]))*(outcome*exp(-b[dim_a+1]*exposure)));
		dudpsi= -1*apply( (iv-mean(iv))*(outcome*exp(-b[dim_a+1]*exposure))*(pvar*design.mat) ,2, mean);
		dudb  = -mean((iv-mean(iv))*as.vector(exposure-expit(design.mat%*%b[1:dim_a]))*exposure*outcome*exp(-b[dim_a+1]*exposure));
		B[(dim_a+n_iv+1),] = c(dudmu, dudpsi, dudb);
	}
	B
}






