Maximisation = function(zvalues, covariates, distances.included, Evar, hypothesis, alternativeDistribution, alternativeCompartmentNumber, dependency, ptol, iter.CG, v)
{
#	if(v) print("M step")
	NUM = length(zvalues)
	f0 = c(0, 1)
	pc  =  rep(1, alternativeCompartmentNumber)/alternativeCompartmentNumber
	mus  =  rep(0, alternativeCompartmentNumber)
	sds  =  rep(1, alternativeCompartmentNumber)
	f1  =  cbind(mus, sds)
	ptheta  =  apply(Evar$gammA,2,sum)/NUM
	
	q6 = 0
	sd0 = 1
	
	if(hypothesis == "two.sided")
	{
		f0  =  c(0,1)
	}
	else
	{
		f0  =  c(0,1,-1)
	}
	
	# f1
	if(alternativeDistribution == "kernel")
	{
		kern.f1  =  density(zvalues,weights=Evar$gammA[,2]/sum(Evar$gammA[,2]))
		f1  =  approx(kern.f1$x, kern.f1$y, zvalues, rule = 2, ties="ordered")$y
	}
	else
	{
		if(alternativeDistribution == "kernel.symetric")
		{
			kern.f1  =  density(c(zvalues,2*f0[1]-zvalues),weights=c(Evar$gammA[,2],Evar$gammA[,2])/sum(c(Evar$gammA[,2],Evar$gammA[,2])))
			f1  =  approx(kern.f1$x, kern.f1$y, zvalues, rule = 2, ties="ordered")$y
		}
		else
		{
			if(alternativeCompartmentNumber == 1)
			{
				q1  =  sum(Evar$gammA[, 2])
				q2  =  sum(Evar$gammA[, 2]*zvalues)
				mu1  =  q2/q1
				q3  =  sum(Evar$gammA[, 2]*(zvalues-mu1)*(zvalues-mu1))
				sd1  =  sqrt(q3/q1)
				f1  =  c(mu1, sd1)
			}
			else
			{
				mus  =  1:alternativeCompartmentNumber
				sds  =  1:alternativeCompartmentNumber
				for (c in 1:alternativeCompartmentNumber)
				{
					q1  =  sum(Evar$omega[, c])
					q2  =  sum(Evar$gammA[, 2])
					pc[c]  =  q1/q2
					q3  =  sum(Evar$omega[, c]*zvalues)
					mus[c]  =  q3/q1
					q4  =  sum(Evar$omega[, c]*(zvalues-mus[c])*(zvalues-mus[c]))
					sds[c]  =  sqrt(q4/q1)
				}
				f1  =  cbind(mus, sds)
			}
		}
	}
	
	if(alternativeDistribution == "kernel" | alternativeDistribution == "kernel.symetric")
	{
		if(hypothesis != "two.sided")
		{
			f1[zvalues == 0] = 0
		}
	}
	if(dependency == "none")
	{
		ptheta = apply(Evar$gammA,2,sum)/NUM
		return(list(ptheta = ptheta, pc = pc, f0 = f0, f1 = f1))
	}
	
	if(length(covariates)>0)
	{
		CG = list()
		CG = tryCatch({ComputeCG(covariates, distances.included, Evar$dgammA, Evar$gammA, Evar$trans.par, iter.CG = iter.CG, ptol, v = v) }, warning = function(e) {print(e)}, error = function(e) {print(e)})
		print("ComputeCG")
		print(CG$pii)
		print(CG$A[,,1])
		print(CG$trans.par)
		CG = tryCatch({ComputeCG.C(covariates, distances.included, Evar$dgammA, Evar$gammA, Evar$trans.par, iter.CG = iter.CG, ptol, v = v) }, warning = function(e) {print(e)}, error = function(e) {print(e)})
		print("ComputeCG.C")
		print(CG$pii)
		print(CG$A[,,1])
		print(CG$trans.par)
		if(length(CG) == 1)
		{
			if(v) cat("Error in CG\n")
			return(-1)
		}
		return(list(pii = CG$pii, ptheta = ptheta, pc=pc, A = CG$A , trans.par = CG$trans.par, f0 = f0, f1 = f1))
	}
	else
	{
		pii =  c(0.95, 0.05)
		for (i in 1:2)
		{
			pii[i]  =  Evar$gammA[1, i]
		}
		
		A  =  array(c(0.95, 0.05, 0.05, 0.95),c(2,2, NUM-1))
		for (i in 1:2)
		{
			for (j in 1:2)
			{ 
				q1  =  sum(Evar$dgammA[i, j, ])
				q2  =  sum(Evar$gammA[1:(NUM-1), i])
				A[i, j,]  =  q1/q2
			}
		}
		return(list(pii = pii, ptheta = ptheta, pc = pc, A = A, f0 = f0, f1 = f1))
	}
	
	
}

