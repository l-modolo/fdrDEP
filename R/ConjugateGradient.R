ComputeCG = function(covariates, distances.included, dgammA, gammA, trans.par, iter.CG, ptol, v)
{
	#if(v) print("CG step")
	# cat("Saving input of ComputeCG\n")
	#save(covariates, distances.included, dgammA, gammA, trans.par, iter.CG, ptol, v, file=paste("/home/vmiele/Collaborations/Modolo/fdrDEP/tests/databenchs/inputComputeCG",dim(covariates)[1],".RData",sep=""))
	#stop("Done")
	gradient.old = ComputeGradient.C(covariates, distances.included, dgammA, gammA, trans.par, v = v)
	if(length(gradient.old) == 1)
	{
		if(v) cat("Error in ComputeGradient\n")
		return(-1)
	}
	
	phi = rbind(rep(0, length(gradient.old)),-gradient.old)
	
	difference = 1
	niter = 0
	while(difference > ptol & niter < iter.CG)
	{
		trans.par.old = trans.par
		niter = niter + 1
		
		tmp = LineSearch.C(covariates, distances.included, dgammA, gammA, trans.par, phi, iter.CG, ptol, v = v)
		if(length(tmp) == 1)
		{
			if(v) cat("Error in LineSearch\n")
			return(-1)
		}
		
		trans.par = tmp$trans.par
		
		gradient.new = ComputeGradient.C(covariates, distances.included, dgammA, gammA, trans.par, v = v)
		if(length(gradient.new) == 1)
		{
			if(v) cat("Error in ComputeGradient\n")
			return(-1)
		}
		
		PR = sum((gradient.new - gradient.old)*gradient.new) / sum(gradient.old^2)
		if(is.nan(PR) | PR < 0)
		{
			PR = 0
		}
		
		phi = rbind(rep(0, length(gradient.new)), -gradient.new + PR * phi[2,])
		gradient.old = gradient.new
		if(distances.included & trans.par[2,4]<0 )
		{
			trans.par[2,4] = abs(trans.par[2,4])
		}
		difference = max(abs(trans.par[2,-1] - trans.par.old[2,-1]))
	}
	if(length(tmp) == 1)
	{
		return(-1)
	}
	else
	{
		tmp = pii.A.C(covariates, distances.included, trans.par, v = v)
		return(list(pii = tmp$pii, A = tmp$A, trans.par = trans.par))
	}
}

ComputeGradient.C = function(covariates, distances.included, dgammA, gammA, trans.par, v)
{
	gradient <- rep(0,3+dim(covariates)[2])
	pii <- rep(0,2)
	A <- array(0,c(2,2,dim(covariates)[1]-1))
	res <- tryCatch({
		.C('C_ComputeGradient',
			as.integer(dim(covariates)[1]),
			as.integer(dim(covariates)[2]),
			as.double(covariates),
			as.integer(distances.included),
			as.double(dgammA),
			as.double(gammA),
			as.double(trans.par),
			as.integer(v),
			as.double(pii),
			as.double(A),
			resgradient = as.double(gradient))
	}, warning = function(e) {return(-1)}, simpleError = function(e) {return(-1)}, error = function(e) {return(-1)})
	if(length(res) == 1)
		return(-1)
	return(res$resgradient)
}

pii.A.C = function(covariates, distances.included, trans.par, v)
{
	pii <- rep(0,2)
	A <- array(0,c(2,2,dim(covariates)[1]-1))
	res <- .C('C_piiA',
		as.integer(dim(covariates)[1]),
		as.integer(dim(covariates)[2]),
		as.double(covariates),
		as.integer(distances.included),
		as.double(trans.par),
		as.integer(v),
		respii = as.double(pii), resA = as.double(A))
	return(list(A = array(res$resA, dim(A)), pii = res$respii))
}

ComputeCG.C = function(covariates, distances.included, dgammA, gammA, trans.par, iter.CG, ptol, v)
{
	pii = rep(0,2)
	A = array(0,c(2,2,dim(covariates)[1]-1))
	res = tryCatch({
		.C('C_ComputeCG',
			as.integer(dim(covariates)[1]),
			as.integer(dim(covariates)[2]),
			as.double(covariates),
			as.integer(distances.included),
			as.double(dgammA),
			as.double(gammA),
			as.integer(iter.CG),
			as.double(ptol),
			as.integer(v),
			restrans.par = as.double(trans.par), respii = as.double(pii), resA = as.double(A))
	}, warning = function(e) {return(-1)}, simpleError = function(e) {return(-1)}, error = function(e) {return(-1)})
	
	if(length(res) == 1)
		return(-1)
	return(list(pii = res$respii, A = array(res$resA, dim(A)), trans.par = matrix(res$restrans.par, nrow=2)))
}

LineSearch.C = function(covariates, distances.included, dgammA, gammA, trans.par, phi, iter.CG, ptol, v)
{
	nu = 0
	pii <- rep(0,2)
	A <- array(0,c(2,2,dim(covariates)[1]-1))
	allright = TRUE
	res = tryCatch({
		.C('C_LineSearch',
			as.integer(dim(covariates)[1]),
			as.integer(dim(covariates)[2]),
			as.double(covariates),
			as.integer(distances.included),
			as.double(dgammA), 
			as.double(gammA),
			restrans.par = as.double(trans.par),
			as.double(phi),
			as.integer(iter.CG),
			as.double(ptol),
			as.integer(v),
			as.integer(allright),
			resnu = as.double(nu), as.double(pii), as.double(A))
	}, warning = function(e) {return(-1)}, simpleError = function(e) {return(-1)}, error = function(e) {return(-1)})
	if(length(res) == 1)
		return(-1)
	return(list(nu = res$resnu, trans.par = matrix(res$restrans.par, nrow=2)))
}

LineSearch = function(Z, dist.included=TRUE, dgamma, gamma, trans.par, phi, iter.CG, ptol, v){

	delta1 = phi[1,]
	delta2 = phi[2,]
	trans.par1 = trans.par[1,]
	trans.par2 = trans.par[2,]
	
	p <- dim(Z)[2]
	# delta1 = (lambda1_H, sigma11_H, sigma21_H, rho1_H)
	delta_H <- array(0,c(2,2,c(dim(Z)[1]-1)))
	delta1_H0 <- delta1[1] + sum(delta1[-c(1:3)]*Z[1,])
	delta2_H0 <- delta2[1] + sum(delta2[-c(1:3)]*Z[1,])
    
	tmp11 <- t(delta1[-c(1:3)]*t(Z[-1,]))
	tmp21 <- t(delta1[-c(1:3)]*t(Z[-1,]))
	tmp12 <- t(delta2[-c(1:3)]*t(Z[-1,]))
	tmp22 <- t(delta2[-c(1:3)]*t(Z[-1,]))

        if(dist.included==TRUE){
          delta1_H0 <- delta1[1] + sum(delta1[-c(1:4)]*Z[1,-1])
          delta2_H0 <- delta2[1] + sum(delta2[-c(1:4)]*Z[1,-1])
  	 
          #tmp11 <- t(c(delta1[4],delta1[-c(1:4)])*t(Z[-1,]))
          tmp21 <- t(c(-delta1[4],delta1[-c(1:4)])*t(Z[-1,]))
          #tmp12 <- t(c(delta2[4],delta2[-c(1:4)])*t(Z[-1,]))
          tmp22 <- t(c(-delta2[4],delta2[-c(1:4)])*t(Z[-1,]))       
        }
        
	denom11 = denom21 = denom12 = denom22 <- rep(0,dim(Z)[1]-1)

	for(i in 1:p){
		denom11 <- denom11 + tmp11[,i]
		denom21 <- denom21 + tmp21[,i]
		denom12 <- denom12 + tmp12[,i]
		denom22 <- denom22 + tmp22[,i]
	}
	
	delta_H[1,1,] <- delta1[2] + denom11
	delta_H[2,1,] <- delta1[3] + denom21
	delta_H[1,2,] <- delta2[2] + denom12
	delta_H[2,2,] <- delta2[3] + denom22	
	
	nu.new <- 0
	difference <- 1
	iter <- 0

	while(difference > ptol & iter < iter.CG){	
		iter <- iter + 1
	
		trans.par1.new <- trans.par1 + nu.new*delta1
		trans.par2.new <- trans.par2 + nu.new*delta2

		tmp.trans.prob = pii.A.C(Z, dist.included, rbind(trans.par1.new, trans.par2.new), v = v)
		
		pii <- tmp.trans.prob$pii
		A <- tmp.trans.prob$A
		
		f <- sum(c(delta1_H0,delta2_H0)*(gamma[1,] - pii))
		#cat("dQ_tmp :", f, "\n")
		for(i in 1:2){
			for(j in 1:2){	
				f <- f + sum(delta_H[i,j,]*(dgamma[i,j,] - gamma[-dim(gamma)[1],i]*A[i,j,]))
			}
		}
		#cat("dQ :", f, "\n")

	     fprime <- -sum(c(delta1_H0,delta2_H0)^2*pii*(1 - pii))
	     #cat("dQ2_tmp :", fprime, "\n")
		for(i in 1:2){
			for(j in 1:2){
				fprime <- fprime - sum(delta_H[i,j,]^2*gamma[-dim(gamma)[1],i]*A[i,j,]*(1 - A[i,j,]))
			}
		}
		#cat("dQ2 :", fprime, "\n")
		nu.old <- nu.new
		nu.new <- nu.old - f/fprime
		difference <- abs(f)
		if(is.na(difference)||iter > 100){
			nu.new <- NaN
			break
		}
	}
	trans.par1.new <- trans.par1 + nu.new*delta1
	trans.par2.new <- trans.par2 + nu.new*delta2
	
	return(list(nu=nu.new,trans.par= rbind(trans.par1.new, trans.par2.new)))
}

ComputeGradient = function(covariates, distances.included, dgammA, gammA, trans.par, v)
{
	gradient <- rep(0,3+dim(covariates)[2])
	
	tmp.trans.prob = pii.A.C(covariates, distances.included, trans.par, v = v)
	
	pii <- tmp.trans.prob$pii
	A <- tmp.trans.prob$A
	
	gradient[1] <- gammA[1,2] - pii[2]	# G1
	gradient[2] <- sum(dgammA[1,2,] - gammA[-dim(gammA)[1],1]*A[1,2,])	# G2
	gradient[3] <- sum(dgammA[2,2,] - gammA[-dim(gammA)[1],2]*A[2,2,])	# G3
	
	tmp <- rep(0,dim(dgammA)[3])
	
	for(i in 1:2){
		tmp <- tmp + gammA[-dim(covariates)[1],i]*A[i,2,] # G422
	}
	gradient[-c(1:3)] <- (gammA[1,2] - pii[2])*covariates[1,] +
	apply(matrix((gammA[-1,2] - tmp)*covariates[-1,],ncol=dim(covariates)[2]),2,sum) # G41 + G42
	
	if(distances.included==TRUE){
		gradient[4] <- (gammA[1,2] - pii[2])*covariates[1,1] +
		sum(dgammA[1,2,]*covariates[-1,1]) -
		sum(dgammA[2,2,]*covariates[-1,1]) -
		sum(gammA[-dim(covariates)[1],1]*A[1,2,]*covariates[-1,1]) +
		sum(gammA[-dim(covariates)[1],2]*A[2,2,]*covariates[-1,1])
	}
	
	return(-gradient)
}

pii.A = function(covariates, distances.included, trans.par, v)
{
	# covariates is a matrix of size m x p
	p <- dim(covariates)[2]
	trans.par1 = trans.par[1,]
	trans.par2 = trans.par[2,]
	
	pii <- rep(0,2)
	A <- array(0,c(2,2,dim(covariates)[1]-1))
	
	pii[1] <- sum(exp(trans.par1[1] + sum(trans.par1[-c(1:3)]*covariates[1,])))/
	( sum(exp(trans.par1[1] + sum(trans.par1[-c(1:3)]*covariates[1,])))
	+ sum(exp(trans.par2[1] + sum(trans.par2[-c(1:3)]*covariates[1,]))) )
	pii[2] = 1 - pii[1]
	
	tmp12 <- t(trans.par2[-c(1:3)]*t(covariates[-1,]))
	tmp22 <- t(trans.par2[-c(1:3)]*t(covariates[-1,]))
	
	if(distances.included==TRUE)
	{
		tmp12 <- t(c(trans.par2[4],trans.par2[-c(1:4)])*t(covariates[-1,]))
		tmp22 <- t(c(-trans.par2[4],trans.par2[-c(1:4)])*t(covariates[-1,]))
	}
	
	denom12 = denom22 <- rep(0,dim(covariates)[1]-1)
	for(i in 1:p)
	{
		denom12 <- denom12 + tmp12[,i]
		denom22 <- denom22 + tmp22[,i]
	}
	
	num11 <- exp(trans.par1[2])
	num21 <- exp(trans.par1[3])
	num12 <- exp(trans.par2[2] + denom12)
	num22 <- exp(trans.par2[3] + denom22)
	
	
	A[1,1,] <- num11/(num11 + num12)
	A[1,2,] <- num12/(num11 + num12)
	A[2,1,] <- num21/(num21 + num22)
	A[2,2,] <- num22/(num21 + num22)
	
	return(list(A = A, pii = pii))
}

pii.A.Call = function(covariates, distances.included, trans.par, v)
{
	pii <- rep(0,2)
	A <- array(0,c(2,2,dim(covariates)[1]-1))
	.Call('Call_piiA',A);
	return(list(A = array(res$resA, dim(A)), pii = res$respii))
}


