fdrDEP = function(pvalues = x, covariates = NULL, distances = NULL, observerdValues = NULL, hypothesis = "two.sided", threshold = 0, alternativeDistribution = 'kernel', alternativeCompartmentNumber = 2, dependency = 'NHMM', seedNumber = 20, burn = 20, ptol = 1e-3, core = 2, maxiter=1000, iter.CG = 1000, v = F, transformed = F)
{
	cat("hypothesis",hypothesis,'\n')
	cat("threshold",threshold,'\n')
	cat("alternativeDistribution",alternativeDistribution,'\n')
	cat("alternativeCompartmentNumber",alternativeCompartmentNumber,'\n')
	cat("dependency",dependency,'\n')
	cat("seedNumber",seedNumber,'\n')
	cat("burn",burn,'\n')
	cat("core",core,'\n')
	cat("maxiter",maxiter,'\n')
	cat("iter.CG",iter.CG,'\n')
	cat("v",v,'\n')
	
	NUM = length(pvalues)
	
	if(dependency != 'NHMM' & dependency != 'HMM' & dependency != 'none'){
		cat('Error: modeltype must be NHMM, HMM or Indep','\n')
		return(0)
	}
	
	if(alternativeDistribution != 'kernel' & alternativeDistribution != 'kernel.symetric' & alternativeDistribution != 'mixnormal'){
		cat('Error: alttype must be kernel, kernel.symetric or mixnormal','\n')
		return(0)
	}
	zvalues = c()
	if(transformed)
		zvalues = pvalues
	else
		zvalues = zvalueTransform(pvalues, hypothesis, threshold, observerdValues)
	
	distances.included = FALSE
	if(length(distances)>0)
	{
		if(length(which(distances<0))>0)
		{
			cat('Error: distances must be positive!','\n')
			return(0)
		}
		if(v) cat('Transforming distances','\n')
		distances <- log2(distances+2)
		distances.included = TRUE
	}
	
	if(length(covariates)>0){
		if(is.vector(covariates)==TRUE)
		{
			covariates = matrix(covariates,ncol=1)
			covariates = scale(covariates)
		}
		else
		{
			covariates = apply(covariates,2,scale)
		}
		if(dim(covariates)[1] != NUM)
		{
			if(v) cat('Error: pvalues and covariates are not compatible', '\n')
			return(-1)
		}
		if(length(distances) > 0 & dim(covariates)[1] != length(distances))
		{
			if(v) cat('Error: distances and covariates are not compatible', '\n')
			return(-1)
		}
		covariates = cbind(distances,covariates)
	}
	else
	{
		covariates = cbind(distances)
	}
	
	working_dir = paste(paste("fdrtmp",format(Sys.time(), "%H_%M_%S"),sep="_"),abs(rnorm(1))*10e15,sep="_")
	dir.create(working_dir)
	ptol = 1e-2
	difference = 1
	niter = 0
	ptol = 1e-2
	difference = 1
	niter = 0
	EMvar = tryCatch({ 
		initialisation(zvalues = zvalues, 
				covariates = covariates, 
				distances.included = distances.included, 
				hypothesis = hypothesis, 
				alternativeDistribution = alternativeDistribution, 
				alternativeCompartmentNumber = alternativeCompartmentNumber, 
				dependency = dependency, 
				seedNumber = seedNumber, 
				burn = burn, 
				ptol = ptol, 
				core = core, 
				maxiter = maxiter, 
				iter.CG = iter.CG, 
				working_dir = working_dir, 
				v = v) }, error = function(e) {print(e);unlink(working_dir,r=T)})
	
	best_EMvar = list()
	while(length(best_EMvar) <= 1)
	{
		niter = niter + 1
		best_EMvar = loadseed(EMvar[[niter]], working_dir)
		if(v) print(paste("best seed : ",best_EMvar$logL))
		best_EMvar = tryCatch({
			ExpectationMaximisation(zvalues = zvalues, 
						covariates = covariates, 
						distances.included = distances.included, 
						Mvar = best_EMvar, 
						hypothesis = hypothesis, 
						alternativeDistribution = alternativeDistribution, 
						alternativeCompartmentNumber = alternativeCompartmentNumber, 
						dependency = dependency, 
						ptol = ptol, 
						maxiter = maxiter, 
						iter.CG = iter.CG, 
						v=v) }, warning = function(e) {print(e)}, error = function(e) {print(e)})
		if(length(EMvar) <= 1)
		{
			if(v) print("error: Trying next best seed")
			best_EMvar = loadseed(EMvar[[niter]], working_dir)
		}
	}
	
	lfdr = best_EMvar$gammA[, 1]
	unlink(working_dir,r=T)
	if(length(best_EMvar) > 1)
	{
		logL  =  best_EMvar$logL
		if(v) print(logL)
		if(alternativeDistribution == "kernel") alternativeCompartmentNumber = 1
		if(alternativeDistribution == "kernel.symetric") alternativeCompartmentNumber = 1
		if (hypothesis != "two.sided") {
			if(dependency == "none")
			{
				BIC = logL - (3*alternativeCompartmentNumber+1)*log(NUM)/2 
			}
			else
			{
				BIC  =  logL - (3*alternativeCompartmentNumber + dim(covariates)[2] + 1)*log(NUM)/2
			}
		} else {
			if(dependency == "none")
			{
				BIC = logL - (3*alternativeCompartmentNumber)*log(NUM)/2 
			}
			else
			{
				BIC  =  logL - (3*alternativeCompartmentNumber + dim(covariates)[2])*log(NUM)/2
			}
		}
		em.var = list(ptheta = best_EMvar$ptheta, pii = best_EMvar$pii, A = best_EMvar$A, pc = best_EMvar$pc, f0 = best_EMvar$f0, f1 = best_EMvar$f1, LIS=lfdr, logL=logL, BIC=BIC, trans.par = best_EMvar$trans.par, gammA = best_EMvar$gammA, dgammA = best_EMvar$dgammA)
		return (em.var)
	} else {
		return(-1)
	}
}

zvalueTransform = function(pvalues = x, hypothesis = "two.sided", threshold = 0, observerdValues = observerdValues)
{
	if( hypothesis == "two.sided")
	{
		zvalues = qnorm( pvalues /2 ,0, 1)
		return( ifelse( observerdValues > threshold, abs(zvalues), -abs(zvalues) ) )
	}
	else
	{
		n = length(pvalues)
		zval = qnorm(pvalues,0,1)
		# FindRoot[Sqrt[2/Pi]/E^(x^2/2) + (-1 - Erf[x/Sqrt[2]])/2 == 0, {x, 0.468021, 0.635139}, WorkingPrecision -> 27]
		root = 0.5179127159921794137280437785
		zvalues = qnorm(pvalues,0,1) + qnorm(root)
		zvalues[zvalues > 0] = 0
		if(n %%2==0)
			return(ifelse(rep(c(TRUE, FALSE), n/2), abs(zvalues), -abs(zvalues)))
		else
			return(ifelse(rep(c(TRUE, FALSE), n/2+1), abs(zvalues), -abs(zvalues))[-(n+1)])
	}
}



