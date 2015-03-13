Maximisation = function(parameters, Evar)
{
	pc = 0
	ptheta = apply(Evar$gammA,2,sum)/(parameters[['NUM']])
	# ptheta = apply(Evar$gammA[parameters[[zvalues]]!=0,],2,sum)/(parameters[['NUM']]*(1-parameters[['delta']]))
	mu0 = 0
	q6 = 0
	sd0 = 1
	
	# q6 = (Evar$gammA[, 1]*parameters[['zvalues']]*parameters[['zvalues']])[parameters[['zvalues']] != 0 & Evar$gammA[, 1] != 0]
	# sd0 = sqrt(sum(q6[is.finite(q6)])/sum(Evar$gammA[, 1][parameters[['zvalues']]!=0 & Evar$gammA[, 1] != 0 & is.finite(q6)]))
	f0  = c(mu0,sd0,-1)
	
	# f1
	if(parameters[['f1']] == "kernel")
	{
		if(parameters[['hypothesis']] == "one.sided")
		{
			# weights = (Evar$gammA[parameters[['zvalues']]!=0,2]/sum(Evar$gammA[parameters[['zvalues']]!=0,2]))
			weights = Evar$gammA[, 2]/sum(Evar$gammA[, 2])
			# kernel epanechnikov are similar to gausian but with a smaller
			# computation time. f1 is not defined in zero thus we estimate it
			# starting from the first z-values superior to zero
			kern.f1  =  density(parameters[['zvalues']], weights=weights, from = min(parameters[['zvalues']][parameters[['zvalues']]!=0]), cut=0, na.rm = T, kernel = "epanechnikov")
			# we set the density of f1 in zero to zero
			kern.f1$x = c(0, kern.f1$x)
			kern.f1$y = c(0, kern.f1$y)
		}
		else
		{
			weights = (Evar$gammA[parameters[['zvalues']]!=0,2]/sum(Evar$gammA[parameters[['zvalues']]!=0,2]))
			kern.f1  =  density(parameters[['zvalues']][parameters[['zvalues']]!=0], weights=weights, cut=0, na.rm = T, kernel = "epanechnikov")
		}
		if(length(kern.f1$x) > 1)
		{
			f1  =  approx(kern.f1$x, kern.f1$y, parameters[['zvalues']], rule = 2, ties="ordered")$y
			f1[parameters[['zvalues']] == 0] = 0
		}
		else
		{
			return(-1)
		}
	}
	else
	{
		if(parameters[['f1']] == "kernel.symetric")
		{
			weights = c(Evar$gammA[,2],Evar$gammA[,2])/sum(c(Evar$gammA[,2],Evar$gammA[,2]))
			weights[is.nan(weights)] = 0
			weights[weights == Inf] = 0
			weights[!is.finite(weights)] = 0
			kern.f1  =  density(c(parameters[['zvalues']],2*f0[1]-parameters[['zvalues']]),weights=weights)
			f1  =  approx(kern.f1$x, kern.f1$y, parameters[['zvalues']], rule = 2, ties="ordered")$y
		}
		else
		{
			pc  =  rep(1, parameters[['f1_compartiments']])/parameters[['f1_compartiments']]
			mus  =  rep(0, parameters[['f1_compartiments']])
			sds  =  rep(1, parameters[['f1_compartiments']])
			f1  =  cbind(mus, sds)
			if(parameters[['f1_compartiments']] == 1)
			{
				q1  =  sum(Evar$gammA[, 2])
				q2  =  sum(Evar$gammA[, 2]*parameters[['zvalues']])
				mu1  =  q2/q1
				q3  =  sum(Evar$gammA[, 2]*(parameters[['zvalues']]-mu1)*(parameters[['zvalues']]-mu1))
				sd1  =  sqrt(q3/q1)
				f1  =  c(mu1, sd1)
			}
			else
			{
				mus  =  1:parameters[['f1_compartiments']]
				sds  =  1:parameters[['f1_compartiments']]
				for (c in 1:parameters[['f1_compartiments']])
				{
					q1  =  sum(Evar$omega[, c])
					q2  =  sum(Evar$gammA[, 2])
					pc[c]  =  q1/q2
					q3  =  sum(Evar$omega[, c]*parameters[['zvalues']])
					mus[c]  =  q3/q1
					q4  =  sum(Evar$omega[, c]*(parameters[['zvalues']]-mus[c])*(parameters[['zvalues']]-mus[c]))
					sds[c]  =  sqrt(q4/q1)
				}
				f1  =  cbind(mus, sds)
			}
		}
	}
	# without dependency
	if(parameters[['dependency']] == "none")
	{
		return(list(ptheta = ptheta, pc = pc, f0 = f0, f1 = f1))
	}
	# with NHMM dependency
	if(length(parameters[['covariates']])>0 & parameters[['dependency']] == "NHMM")
	{
		CG = ComputeCG(parameters[['covariates']], parameters[['distances_included']], Evar$dgammA, Evar$gammA, Evar$trans.par, iter.CG = parameters[['GC_max_iteration']], parameters[['EM_threshold']], v = parameters[['v']])
		if(length(CG) == 1)
		{
			if(parameters[['v']]) cat("Error in CG\n")
			return(-1)
		}
		return(list(pii = CG$pii, ptheta = ptheta, pc=pc, A = CG$A , trans.par = CG$trans.par, f0 = f0, f1 = f1))
	}
	else # with HMM dependency
	{
		pii =  c(-1, -1)
		for (i in 1:2)
		{
			pii[i]  =  Evar$gammA[1, i]
		}
		
		A  =  array(c(-1, -1, -1, -1),c(2,2, parameters[['NUM']]-1))
		for (i in 1:2)
		{
			for (j in 1:2)
			{ 
				q1  =  sum(Evar$dgammA[i, j, ])
				q2  =  sum(Evar$gammA[1:(parameters[['NUM']]-1), i])
				A[i, j,]  =  q1/q2
			}
		}
		return(list(pii = pii, ptheta = ptheta, pc = pc, A = A, f0 = f0, f1 = f1))
	}
}

