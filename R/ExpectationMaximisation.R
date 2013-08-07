ExpectationMaximisation = function(zvalues, covariates, distances.included, Mvar, hypothesis, alternativeDistribution, alternativeCompartmentNumber, dependency, ptol,  maxiter, iter.CG, v)
{
	NUM = length(zvalues)
	difference = 1
	converged=TRUE
	niter = 0
	logL.old = -Inf
	logL = 0
	T = 1
	Evar = list()
	Mvar.old = Mvar
	
	while(difference > ptol && niter <= maxiter)
	{
		niter = niter + 1
		Evar     = Expectation(zvalues = zvalues, 
				Mvar = Mvar.old, 
				hypothesis = hypothesis, 
				alternativeDistribution = alternativeDistribution, 
				alternativeCompartmentNumber = alternativeCompartmentNumber, 
				dependency = dependency, 
				v = v)
		if( length(Evar) == 1)
		{
			if(v) cat('Error in E', '\n')
			converged=FALSE;
			break
		}
		logL  = -sum(log(Evar$c0))
		
		if(dependency == "none")
		{
			logL = -logL
		}
		if( logL < logL.old & abs(logL - logL.old) > 0.1 & niter > 2)
		{
			if(v) cat('Error in EM : logL increasing by', logL.old - logL , 'iteration :', niter, '\n')
			converged=FALSE
			break
		}
		
		Mvar     = Maximisation(zvalues = zvalues, 
				covariates = covariates, 
				distances.included = distances.included, 
				Evar = Evar, 
				hypothesis = hypothesis, 
				alternativeDistribution = alternativeDistribution, 
				alternativeCompartmentNumber = alternativeCompartmentNumber, 
				dependency = dependency, 
				iter.CG = iter.CG, 
				ptol = ptol, 
				v = v)
		if( length(Mvar) == 1)
		{
			if(v) cat('Error in M', '\n')
			converged=FALSE
			break
		}
		
		if(dependency == "none")
		{
			df2 <- abs(Mvar.old$f1 - Mvar$f1)
			difference <- max(df2)
		}
		else
		{
			if(length(covariates)>0)
			{
				df1 = abs(Mvar.old$trans.par[2,-1] - Mvar$trans.par[2,-1])
				df2 = abs(Mvar.old$f1 - Mvar$f1)
				df3  =  ifelse(abs(logL - logL.old) > 0.1, abs(logL - logL.old), 0)
			}
			else
			{
				df1 = abs(Mvar.old$A - Mvar$A)
				df2 = abs(Mvar.old$f1 - Mvar$f1)
				df3  =  ifelse(abs(logL - logL.old) > 0.1, abs(logL - logL.old), 0)
			
			}
			difference = max(df1, df2, df3)
		}
		
		if( is.na(difference) )
		{
			if(v) cat('Error in EM : NA result', '\n')
			converged=FALSE
			break
		}
		logL.old = logL
		Mvar.old = Mvar
	}
	if(!converged)
	{
		return(-1)
	}
	if(v) cat('done in ', niter, ' iterations', '\n')
	if(dependency == "none")
	{
		return(list(ptheta = Mvar.old$ptheta, pc=Mvar.old$pc, f0=Mvar.old$f0, f1=Mvar.old$f1, logL = logL, gammA = Evar$gammA))
	}
	if(length(covariates)>0)
	{
		return(list(pii=Mvar.old$pii, ptheta = Mvar.old$ptheta, pc=Mvar.old$pc, A=Mvar.old$A, trans.par = Mvar.old$trans.par, f0=Mvar.old$f0, f1=Mvar.old$f1, logL = logL, gammA = Evar$gammA))
	}
	else
	{
		return(list(pii=Mvar$pii, ptheta = Mvar$ptheta, pc=Mvar$pc, A=Mvar$A, f0=Mvar$f0, f1=Mvar$f1, logL = logL, gammA = Evar$gammA))
	}
}
