Expectation = function(zvalues, Mvar, hypothesis, alternativeDistribution, alternativeCompartmentNumber, dependency, v)
{
#	if(v) print("E step")
	res = list()
	NUM = length(zvalues)
	delta = length(zvalues[zvalues==0])/length(zvalues)
	f0x = c()
	f1x = c()
	gammA = matrix(rep(0, NUM*2), NUM, 2, byrow=TRUE)
	omega = c()
	# f1x
	if(alternativeDistribution == "kernel" | alternativeDistribution == "kernel.symetric")
	{
		f1x = Mvar$f1
	}
	else
	{
		f1x = rep(0, NUM)
		if(alternativeCompartmentNumber == 1)
		{
			f1x = dnorm(zvalues, Mvar$f1[1], Mvar$f1[2])
		}
		else
		{
			for (c in 1:alternativeCompartmentNumber)
			{
				f1x = f1x+Mvar$pc[c]*dnorm(zvalues, Mvar$f1[c, 1], Mvar$f1[c, 2])
			}
		}
	}
	
	# f0x
	if(length(Mvar$f0) < 3)
	{
		f0x = dnorm(zvalues, Mvar$f0[1], Mvar$f0[2])
	}
	else
	{
		f0x = delta * (zvalues==0) + (1-delta)*dnorm(zvalues, Mvar$f0[1], Mvar$f0[2])# * (x!=0)
		f1x[zvalues==0] = 0
	}
	
	if(dependency == "none")
	{
		gammA[,1]  =  Mvar$ptheta[1]*f0x / (Mvar$ptheta[1]*f0x + Mvar$ptheta[2]*f1x)
		gammA[,2]  =  1 - gammA[,1]
		
		if(alternativeDistribution != "kernel" & alternativeDistribution != "kernel.symetric" & alternativeCompartmentNumber > 1)
		{
			# omega
			omega = matrix(rep(0, NUM*alternativeCompartmentNumber), NUM, alternativeCompartmentNumber, byrow=TRUE)
			for (c in 1:alternativeCompartmentNumber)
			{
				f1c = dnorm(zvalues, Mvar$f1[c, 1], Mvar$f1[c, 2])
				omega[, c] = gammA[, 2] * Mvar$pc[c]*f1c/f1x
				if(length(Mvar$f0) >= 3)
				{
					omega[zvalues==0, c] = 0
				}
			}
		}
		c0  =  f0x*Mvar$ptheta[1] + f1x*Mvar$ptheta[2]
		return(list(gammA = gammA, omega = omega, c0 = c0))
	}
	
	# forward
	alpha = matrix(rep(0, NUM*2), NUM, 2, byrow=TRUE)
	c0 = rep(0, NUM)
	alpha[1, 1] = Mvar$pii[1]*f0x[1]
	alpha[1, 2] = Mvar$pii[2]*f1x[1]
	c0[1] = 1/sum(alpha[1, ])
	alpha[1, ] = c0[1]*alpha[1, ]
	alpha.tmp  =  .C('calAlpha',alpha=as.numeric(alpha),c0=as.numeric(c0),as.numeric(Mvar$A),as.numeric(f0x),as.numeric(f1x),as.integer(NUM))
	alpha  =  alpha.tmp$alpha
	dim(alpha)  =  c(NUM,2)
	c0  =  alpha.tmp$c0
	
	# backward
	beta = matrix(rep(0, NUM*2), NUM, 2, byrow=TRUE)
	beta[NUM, 1] = c0[NUM]
	beta[NUM, 2] = c0[NUM]
	beta.tmp  =  .C('calBeta',beta=as.numeric(beta),as.numeric(c0),as.numeric(Mvar$A),as.numeric(f0x),as.numeric(f1x),as.integer(NUM))
	beta  =  beta.tmp$beta
	dim(beta)  =  c(NUM,2)
	
	# lfdr
	lfdr = rep(0, NUM)
	lfdr.tmp  =  .C('calLfdr',as.numeric(alpha),as.numeric(beta),lfdr=as.numeric(lfdr),as.integer(NUM))
	lfdr  =  lfdr.tmp$lfdr
	
	# gammA & dgammA
	gammA = matrix(rep(0,(NUM*2)), NUM, 2, byrow=TRUE)
	gammA[NUM, ] = c(lfdr[NUM], 1-lfdr[NUM])
	dgammA = array(rep(0, (NUM-1)*4), c(2, 2, (NUM-1)))
	gammA.tmp  =  .C('calGamma',as.numeric(alpha),as.numeric(beta),as.numeric(Mvar$A),as.numeric(f0x),as.numeric(f1x),gamma=as.numeric(gammA),dgamma=as.numeric(dgammA),as.integer(NUM))
	gammA  =  gammA.tmp$gamma
	dgammA  =  gammA.tmp$dgamma
	dim(gammA)  =  c(NUM,2)
	dim(dgammA)  =  c(2, 2, (NUM-1))
	
	if(alternativeDistribution != "kernel" & alternativeDistribution != "kernel.symetric" & alternativeCompartmentNumber > 1)
	{
		# omega
		omega = matrix(rep(0, NUM*alternativeCompartmentNumber), NUM, alternativeCompartmentNumber, byrow=TRUE)
		for (c in 1:alternativeCompartmentNumber)
		{ 
			f1c = dnorm(zvalues, Mvar$f1[c, 1], Mvar$f1[c, 2])
			omega[, c] = gammA[, 2] * Mvar$pc[c]*f1c/f1x
			if(length(Mvar$f0) >= 3)
			{
				omega[zvalues==0, c] = 0
			}
		}
	}
	
	return(list(gammA = gammA, dgammA = dgammA, omega = omega, c0 = c0, trans.par = Mvar$trans.par))
}


