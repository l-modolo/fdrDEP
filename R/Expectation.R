Expectation = function(parameters, Mvar)
{
	res = list()
	f0x = c()
	f1x = c()
	gammA = matrix(rep(0, parameters[['NUM']]*2), parameters[['NUM']], 2, byrow=TRUE)
	omega = c()

	f0x = (2*dnorm(parameters[['zvalues']][parameters[['zvalues']]!=0], Mvar$f0[1], Mvar$f0[2]))
	f0x[parameters[['zvalues']]==0] = 1
	
	# f1x
	if(parameters[['f1']] == "kernel" | parameters[['f1']] == "kernel.symetric")
	{
		f1x = Mvar$f1
	}
	else
	{
		f1x = rep(0, parameters[['NUM']])
		if(parameters[['f1_compartiments']] == 1)
		{
			f1x = dnorm(parameters[['zvalues']], Mvar$f1[1], Mvar$f1[2])
		}
		else
		{
			for (c in 1:parameters[['f1_compartiments']])
			{
				f1x = f1x+Mvar$pc[c]*dnorm(parameters[['zvalues']], Mvar$f1[c, 1], Mvar$f1[c, 2])
			}
		}
	}
	
	if(parameters[['dependency']] == "none")
	{
		gammA[,1]  =  ifelse( f0x[parameters[['zvalues']]==0], 1, Mvar$ptheta[1]*f0x / (Mvar$ptheta[1]*f0x + Mvar$ptheta[2]*f1x) )
		gammA[,2]  =  1 - gammA[,1]
		if(parameters[['f1']] != "kernel" & parameters[['f1']] != "kernel.symetric" & parameters[['f1_compartiments']] > 1)
		{
			# omega
			omega = matrix(rep(0, parameters[['NUM']]*parameters[['f1_compartiments']]), parameters[['NUM']], parameters[['f1_compartiments']], byrow=TRUE)
			for (c in 1:parameters[['f1_compartiments']])
			{
				f1c = dnorm(parameters[['zvalues']], Mvar$f1[c, 1], Mvar$f1[c, 2])
				omega[, c] = gammA[, 2] * Mvar$pc[c]*f1c/f1x
				if(length(Mvar$f0) >= 3)
				{
					omega[parameters[['zvalues']]==0, c] = 0
				}
			}
		}
		c0  =  f0x*Mvar$ptheta[1] + f1x*Mvar$ptheta[2]
		return(list(gammA = gammA, omega = omega, c0 = c0))
	}
	
	# forward
	alpha = matrix(rep(0, parameters[['NUM']]*2), parameters[['NUM']], 2, byrow=TRUE)
	c0 = rep(0, parameters[['NUM']])
	alpha[1, 1] = Mvar$pii[1]*f0x[1]
	alpha[1, 2] = Mvar$pii[2]*f1x[1]
	c0[1] = 1/sum(alpha[1, ])
	alpha[1, ] = c0[1]*alpha[1, ]
	alpha.tmp  =  tryCatch({
				.C('calAlpha',alpha=as.numeric(alpha),c0=as.numeric(c0),as.numeric(Mvar$A),as.numeric(f0x),as.numeric(f1x),as.integer(parameters[['NUM']]))
			 }, warning = function(e) {return(-1)}, simpleError = function(e) {return(-1)}, error = function(e) {return(-1)})
	if(length(alpha.tmp) == 1)
	{
		if(parameters[['v']]) cat("Error in calAlpha\n")
		return(-1)
	}
	alpha  =  alpha.tmp$alpha
	dim(alpha)  =  c(parameters[['NUM']],2)
	c0  =  alpha.tmp$c0
	
	# backward
	beta = matrix(rep(0, parameters[['NUM']]*2), parameters[['NUM']], 2, byrow=TRUE)
	beta[parameters[['NUM']], 1] = c0[parameters[['NUM']]]
	beta[parameters[['NUM']], 2] = c0[parameters[['NUM']]]
	beta.tmp  =  tryCatch({
				.C('calBeta',beta=as.numeric(beta),as.numeric(c0),as.numeric(Mvar$A),as.numeric(f0x),as.numeric(f1x),as.integer(parameters[['NUM']]))
			 }, warning = function(e) {return(-1)}, simpleError = function(e) {return(-1)}, error = function(e) {return(-1)})
	if(length(beta.tmp) == 1)
	{
		if(parameters[['v']]) cat("Error in calBeta\n")
		return(-1)
	}
	beta  =  beta.tmp$beta
	dim(beta)  =  c(parameters[['NUM']],2)
	
	# lfdr
	lfdr = rep(0, parameters[['NUM']])
	lfdr.tmp  =  tryCatch({
				.C('calLfdr',as.numeric(alpha),as.numeric(beta),lfdr=as.numeric(lfdr),as.integer(parameters[['NUM']]))
			 }, warning = function(e) {return(-1)}, simpleError = function(e) {return(-1)}, error = function(e) {return(-1)})
	if(length(lfdr.tmp) == 1)
	{
		if(parameters[['v']]) cat("Error in calLfdr\n")
		return(-1)
	}
	lfdr  =  lfdr.tmp$lfdr
	
	# gammA & dgammA
	gammA = matrix(rep(0,(parameters[['NUM']]*2)), parameters[['NUM']], 2, byrow=TRUE)
	gammA[parameters[['NUM']], ] = c(lfdr[parameters[['NUM']]], 1-lfdr[parameters[['NUM']]])
	dgammA = array(rep(0, (parameters[['NUM']]-1)*4), c(2, 2, (parameters[['NUM']]-1)))
	gammA.tmp  =  tryCatch({
				.C('calGamma',as.numeric(alpha),as.numeric(beta),as.numeric(Mvar$A),as.numeric(f0x),as.numeric(f1x),gamma=as.numeric(gammA),dgamma=as.numeric(dgammA),as.integer(parameters[['NUM']]))
			 }, warning = function(e) {return(-1)}, simpleError = function(e) {return(-1)}, error = function(e) {return(-1)})
	if(length(gammA.tmp) == 1)
	{
		if(parameters[['v']]) cat("Error in calGamma\n")
		return(-1)
	}
	gammA  =  gammA.tmp$gamma
	dgammA  =  gammA.tmp$dgamma
	dim(gammA)  =  c(parameters[['NUM']],2)
	dim(dgammA)  =  c(2, 2, (parameters[['NUM']]-1))
	
	if(parameters[['f1']] != "kernel" & parameters[['f1']] != "kernel.symetric" & parameters[['f1_compartiments']] > 1)
	{
		# omega
		omega = matrix(rep(0, parameters[['NUM']]*parameters[['f1_compartiments']]), parameters[['NUM']], parameters[['f1_compartiments']], byrow=TRUE)
		for (c in 1:parameters[['f1_compartiments']])
		{ 
			f1c = dnorm(parameters[['zvalues']], Mvar$f1[c, 1], Mvar$f1[c, 2])
			omega[, c] = gammA[, 2] * Mvar$pc[c]*f1c/f1x
			if(length(Mvar$f0) >= 3)
			{
				omega[parameters[['zvalues']]==0, c] = 0
			}
		}
	}
	return(list(gammA = gammA, dgammA = dgammA, omega = omega, c0 = c0, trans.par = Mvar$trans.par))
}


