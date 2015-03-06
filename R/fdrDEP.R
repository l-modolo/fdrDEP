fdrDEP = function(
	pvalues, 
	covariates = NULL, 
	distances= NULL, 
	hypothesis = "two.sided", 
	side = NULL,
	mu = -1,
	pi_0 = -1,
	alternative_distribution = "kernel", 
	alternative_compartment_number = 2, 
	dependency_structure = "none", 
	seed_number = 20, 
	EM_burn_iterations = 20, 
	EM_threshold = 1e-3, 
	EM_max_iteration = 1000, 
	GC_max_iteration = 1000, 
	thread_number = 1,
	working_dir = "fdrDEP_tmp",
	v = FALSE)
{
	parameters = list()
	parameters[['pvalues']] = pvalues
	parameters[['covariates']] = covariates
	parameters[['distances']] = distances
	parameters[['hypothesis']] = hypothesis
	parameters[['side']] = side
	parameters[['f1']] = alternative_distribution
	parameters[['f1_compartiments']] = alternative_compartment_number
	parameters[['dependency']] = dependency_structure
	parameters[['seed_number']] = seed_number
	parameters[['EM_burn_iterations']] = EM_burn_iterations
	parameters[['EM_threshold']] = EM_threshold
	parameters[['EM_max_iteration']] = EM_max_iteration
	parameters[['GC_max_iteration']] = GC_max_iteration
	parameters[['thread_number']] = thread_number
	parameters[['working_dir']] = working_dir
	parameters[['v']] = v
	parameters[['A11']] = 0.95
	parameters[['A22']] = 0.2
	parameters[['NUM']] = length(parameters[['pvalues']])
	parameters[['distances_included']] = FALSE
	parameters[['lambda']] = 0
	parameters[['mu']] = mu
	parameters[['pi_0']] = pi_0

	parameters = parameters_error(parameters)
	parameters = probit_intern(parameters)
	parameters = covariates_normalisation(parameters)

	seedList_init = initialisation(parameters)
	best_EMvar = list()
	if(parameters[['seed_number']] > 1)
	{
		seed_number = 1
		while(length(best_EMvar) <= 1 & seed_number <= length(seedList_init)+1)
		{
			if(seed_number == length(seedList_init) + 1)
			{
				print("error in all seed, recomputation...")
				seedList_init = initialisation(parameters)
				seed_number = 1
			}
			load((seedList_init[[seed_number]])$file)
			if(v) print(paste("best seed : ",(seedList_init[[seed_number]])$logL))
			best_EMvar = ExpectationMaximisation(parameters, Mvar = seedList)
			if(length(best_EMvar) <= 1) if(v) print("error: Trying next best seed")
			seed_number = seed_number + 1
		}
	}
	else
	{
		load((seedList_init[[1]])$file)
		best_EMvar = ExpectationMaximisation(parameters, Mvar = seedList_init)
	}
	lfdr = best_EMvar$gammA[, 1]
	system(paste(" ps -ef | grep \"R --slave\" | awk '{if($3 == ",Sys.getpid(),")system(\"kill \"$2)}'", sep=""))
	if(length(best_EMvar) > 1)
	{
		logL  =  best_EMvar$logL
		if(v) print(logL)
		if(parameters[['f1']] == "kernel") parameters[['f1_compartiments']] = 1
		if(parameters[['f1']] == "kernel.symetric") parameters[['f1_compartiments']] = 1
		if (hypothesis != "two.sided") {
			if(parameters[['dependency']] == "none")
			{
				BIC = logL - (3*parameters[['f1_compartiments']]+1)*log(parameters[['NUM']])/2 
			}
			else
			{
				BIC = logL - (3*parameters[['f1_compartiments']] + dim(covariates)[2] + 1)*log(parameters[['NUM']])/2
			}
		} else {
			if(parameters[['dependency']] == "none")
			{
				BIC = logL - (3*parameters[['f1_compartiments']])*log(parameters[['NUM']])/2 
			}
			else
			{
				BIC = logL - (3*parameters[['f1_compartiments']] + dim(covariates)[2])*log(parameters[['NUM']])/2
			}
		}
		em.var = list(ptheta = best_EMvar$ptheta, pii = best_EMvar$pii, A = best_EMvar$A, pc = best_EMvar$pc, f0 = best_EMvar$f0, f1 = best_EMvar$f1, LIS=lfdr, logL=logL, BIC=BIC, trans.par = best_EMvar$trans.par, gammA = best_EMvar$gammA, dgammA = best_EMvar$dgammA, zvalues = parameters[['zvalues']], mu = parameters[['mu']], lambda = parameters[['lambda']], pi_0 = parameters[['pi_0']])
		return (em.var)
	} else {
		return(-1)
	}
}

probit = function(pvalues, hypothesis = "one.sided", mu = -1, pi_0 = -1, thread_number = 1, side = NULL, parameters = list())
{
	parameters[['pvalues']] = pvalues
	parameters[['hypothesis']] = hypothesis
	parameters[['mu']] = mu
	parameters[['pi_0']] = pi_0
	parameters[['thread_number']] = thread_number
	parameters[['side']] = side
	return(probit_intern(parameters))
}

probit_intern = function(parameters)
{
	if(parameters[['hypothesis']] == "one.sided")
	{
		if(parameters[['mu']] == -1 | parameters[['pi_0']] == -1)
		{
			LPO_results = LPO(parameters[['pvalues']], p = 1, thread_number = parameters[['thread_number']])
			parameters[['lambda']] = LPO_results$lambda
			parameters[['mu']] = LPO_results$mu
			parameters[['pi_0']] = LPO_results$pi_0
		}
		parameters[['zvalues']] = probit_one_sided(parameters)
	}
	else
	{
		parameters[['zvalues']] = probit_two_sided(parameters)
		parameters[['mu']] = 1
		parameters[['pi_0']] = 0.8
		parameters[['lambda']] = 0
	}
	return(parameters)
}

probit_one_sided = function(parameters)
{
	tmp = parameters[['pvalues']]
	tmp[tmp >= parameters[['mu']]] = parameters[['mu']]
	tmp = tmp/parameters[['mu']]
	return(qnorm(1-tmp/2, 0, 1))
}

probit_two_sided = function(parameters)
{
	return(ifelse(parameters[['side']] == "+", qnorm(1-parameters[['pvalues']]/2, 0, 1), qnorm(parameters[['pvalues']]/2, 0, 1)))
}

LIS_graph = function (zval, LIS, k=F, title="", hyp="one.sided")
{
	if(hyp == "one.sided")
	{
		hist(zval[zval != 0], breaks = sqrt(length(zval)), freq=F, xlab="z-values", main="")
		legend("topright", legend=c("f0", "f1", "f"), 
		col=c("blue", "red", "black"), lty=1, lwd=4, cex=1, bty="n")
		x_tmp = order(zval)
		zval_tmp = zval[x_tmp]
		f1 = LIS$ptheta[2] * LIS$f1[x_tmp]
		f0 = LIS$ptheta[1] * dnorm(zval_tmp, LIS$f0[1], LIS$f0[2]) * 2
		lines(zval_tmp, f1 + f0, col="black", lwd=1)
		lines(zval_tmp, f0, col="blue", lwd=1)
		lines(zval_tmp, f1, col="red", lwd=1)
		lines(zval_tmp, LIS$LIS[order(zval)], col="blue", lwd=1)
		lines(zval_tmp, (1-LIS$LIS[order(zval)]), col="red", lwd=1)
		abline(v=0, col="green")
		# plot(zval_tmp, LIS$LIS[x_tmp], col="red", lwd=1, type="l")
	}
	else
	{
		hist(zval[zval != 0], nclass=sqrt(length(zval)), main=title, freq=F, xlim=c(-10,10), ylim=c(0,1), xlab="z-value", cex.main = 3, cex.lab = 2, cex.axis = 2.5)
		
		x_tmp = order(zval)
		zval_tmp = zval[x_tmp]
		f0 = c()
		f1 = c()
		if(length(LIS$f0) == 2)
		{
			f0 = LIS$ptheta[1] * dnorm(zval_tmp, LIS$f0[1], LIS$f0[2])
		}
		if(length(LIS$f0) == 3)
		{
			f0 = LIS$ptheta[1] * dnorm(zval_tmp, LIS$f0[1], LIS$f0[2])
		}
		if(length(LIS$f0) == 4)
		{
			f0 = LIS$ptheta[1] * (LIS$pnu[1] * dnorm(zval_tmp, LIS$f0[1], LIS$f0[2]) + LIS$pnu[2] * dnorm(zval_tmp, LIS$f0[1], LIS$f0[4] ))
		}
		
		if(k)
		{
			f1 = LIS$ptheta[2] * LIS$f1[x_tmp]
		}
		else
		{
			f1 = rep(0, length(zval))
			for(ell in 1:dim(LIS$f1)[1])
			{
				f1l = LIS$ptheta[2] * LIS$pc[ell]  * dnorm(zval_tmp, LIS$f1[ell,1], LIS$f1[ell,2])
				abline(v=LIS$f1[ell,1], col="red")
				f1 = f1 + f1l
			}
		}
		lines(zval_tmp, f1 + f0, col="green", lwd=2)
		lines(zval_tmp, f0, col="blue", lwd=3)
		lines(zval_tmp, f1, col="red", lwd=3)
		lines(zval_tmp, LIS$LIS[order(zval)], col="blue", lwd=1)
		lines(zval_tmp, (1-LIS$LIS[order(zval)]), col="red", lwd=1)
		plot(zval_tmp, LIS$LIS[x_tmp], col="red", lwd=1, type="l")
	}
}

covariates_normalisation = function(parameters)
{
	if(length(parameters[['distances']])>0)
	{
		parameters[['distances_included']] = TRUE
		if(parameters[['distances_included']] & length(parameters[['distances']]) > 0)
		{
			if(parameters[['v']]){cat('Transforming distances','\n')}
			parameters[['distances']] = log2(parameters[['distances']]+2)
		}
	}
	if(length(parameters[['covariates']])>0){
		if(is.vector(parameters[['covariates']])==TRUE)
		{
			parameters[['covariates']] = matrix(parameters[['covariates']],ncol=1)
			parameters[['covariates']] = scale(parameters[['covariates']])
		}
		else
		{
			parameters[['covariates']] = apply(parameters[['covariates']],2,scale)
		}
		parameters[['covariates']] = cbind(parameters[['distances']],parameters[['covariates']])
	}
	else
	{
		parameters[['covariates']] = cbind(parameters[['distances']])
	}
	return(parameters)
}

parameters_error = function(parameters)
{
	if( !(parameters[['hypothesis']] %in% c("two.sided", "one.sided")) )
	{
		cat('Error: hypothesis must be two.sided, one.sided','\n')
		return(0)
	}
	else
		if(parameters[['v']]) {cat("hypothesis ", parameters[['hypothesis']], '\n')}
	
	if( !(parameters[['dependency']] %in% c("NHMM", "HMM", "none")) )
	{
		cat('Error: modeltype must be NHMM, HMM or none','\n')
		return(0)
	}
	else
		if(parameters[['v']]) {cat("dependency ", parameters[['dependency']], '\n')}
	if(length(parameters[['distances']])>0)
	{
		if(length(which(parameters[['distances']]<0))>0)
		{
			cat('Error: distances must be positive!','\n')
			return(0)
		}
		parameters[['distances_included']] = TRUE
	}
	if(length(parameters[['covariates']]) > 0)
	{
		if(dim(parameters[['covariates']])[1] != parameters[['NUM']])
		{
			if(v) cat('Error: pvalues and covariates are not compatible', '\n')
			return(-1)
		}
		if(length(parameters[['distances']]) > 0 & dim(parameters[['covariates']])[1] != length(parameters[['distances']]))
		{
			if(v) cat('Error: distances and covariates are not compatible', '\n')
			return(-1)
		}
	}

	if( !(parameters[['f1']] %in% c('kernel', 'kernel.symetric','mixnormal')) )
	{
		cat('Error: f1 must be kernel, kernel.symetric or mixnormal','\n')
		return(0)
	}
	else
		if(parameters[['v']]) {cat("f1 ", parameters[['f1']], '\n')}

	if( parameters[['f1']] == 'mixnormal')
		if(parameters[['v']]) {cat("f1 number of ompartiments", parameters[['f1_compartiments']], '\n')}

	return(parameters)
}
