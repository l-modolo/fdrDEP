ExpectationMaximisation = function(parameters, Mvar, burn = FALSE)
{
	difference = 1
	converged=TRUE
	niter = 0
	logL.old = -Inf
	logL = 0
	T = 1
	Mvar.old = Mvar
	Evar = list()

	while(converged & difference > parameters[['EM_threshold']] & niter <= ifelse(burn, parameters[['EM_burn_iterations']], parameters[['EM_max_iteration']]))
	{
		niter = niter + 1
		Evar     = Expectation(parameters = parameters, Mvar = Mvar.old)
		if( length(Evar) == 1)
		{
			if(parameters[['v']]) cat('Error in E', '\n')
			converged=FALSE;
		}
		logL  = -sum(log(Evar$c0))
		
		if(parameters[['dependency']] == "none")
		{
			logL = -logL
		}
		if(parameters[['v']]) print(logL)
		if( logL < logL.old & abs(logL - logL.old) > 0.1 & niter > 2)
		{
			if(parameters[['v']]) cat('Error in EM : logL increasing by', logL.old - logL , 'iteration :', niter, '\n')
			converged=FALSE
		}
		
		Mvar     = Maximisation(parameters = parameters, Evar = Evar)
		if( length(Mvar) == 1)
		{
			if(parameters[['v']]) cat('Error in M', '\n')
			converged=FALSE
		}
		
		if(parameters[['dependency']] == "none")
		{
			df2 <- abs(Mvar.old$f1 - Mvar$f1)
			difference <- max(df2)
		}
		else
		{
			if(length(parameters[['covariates']])>0)
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
			if(parameters[['v']]) cat('Error in EM : NA result', '\n')
			converged=FALSE
		}
		logL.old = logL
		Mvar.old = Mvar
	}
	if(!converged)
	{
		print("not converged")
		return(-1)
	}
	if(parameters[['v']]) cat('done in ', niter, ' iterations', '\n')
	if(parameters[['dependency']] == "none")
	{
		return(list(ptheta = Mvar.old$ptheta, pc=Mvar.old$pc, f0=Mvar.old$f0, f1=Mvar.old$f1, logL = logL, gammA = Evar$gammA))
	}
	if(length(parameters[['covariates']])>0)
	{
		return(list(pii=Mvar.old$pii, ptheta = Mvar.old$ptheta, pc=Mvar.old$pc, A=Mvar.old$A, trans.par = Mvar.old$trans.par, f0=Mvar.old$f0, f1=Mvar.old$f1, logL = logL, gammA = Evar$gammA))
	}
	else
	{
		return(list(pii=Mvar$pii, ptheta = Mvar$ptheta, pc=Mvar$pc, A=Mvar$A, f0=Mvar$f0, f1=Mvar$f1, logL = logL, gammA = Evar$gammA))
	}
}
