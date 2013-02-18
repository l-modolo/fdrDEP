initialisation = function(zvalues, covariates, distances.included, hypothesis, alternativeDistribution, alternativeCompartmentNumber, dependency, seedNumber, burn, ptol, core, maxiter, iter.CG, working_dir, v)
{
	seedList = mclapply(1:seedNumber, FUN = function(x, zvalues, covariates, distances.included, hypothesis, alternativeDistribution, alternativeCompartmentNumber, dependency, seedNumber, burn, ptol, core, maxiter, iter.CG, working_dir, v){
			runseed(iter = x, 
				zvalues = zvalues, 
				covariates = covariates, 
				distances.included = distances.included, 
				hypothesis = hypothesis, 
				alternativeDistribution = alternativeDistribution, 
				alternativeCompartmentNumber = alternativeCompartmentNumber,
				dependency = dependency, 
				seedNumber = seedNumber, 
				burn = burn, 
				ptol = ptol, 
				maxiter = maxiter, 
				iter.CG = iter.CG, 
				working_dir = working_dir, 
				v = v)}, 
		mc.cores = core, 
		zvalues = zvalues, 
		covariates = covariates, 
		distances.included = distances.included, 
		hypothesis = hypothesis, 
		alternativeDistribution = alternativeDistribution, 
		alternativeCompartmentNumber = alternativeCompartmentNumber,
		dependency = dependency, 
		seedNumber = seedNumber, 
		burn = burn, 
		ptol = ptol, 
		maxiter = maxiter, 
		iter.CG = iter.CG, 
		working_dir = working_dir, 
		v = v)
	
	logL = c()
	for(i in 1:seedNumber)
	{
		logL[i] = ifelse(length(seedList[[i]]) > 1, (seedList[[i]])$logL, -Inf)
	}
	seedList_ordered = list()
	j = 1
	for(i in order(logL, decreasing = T))
	{
		seedList_ordered [[j]] = seedList[[i]]
		j = j + 1
	}
	
	return(seedList_ordered)
}

loadseed = function(seedList, working_dir)
{
	load(paste(paste(working_dir, "seed_",sep="/"),seedList$file,sep=""))
	return(seedList)
}


runseed = function(iter, zvalues, covariates, distances.included, hypothesis, alternativeDistribution, alternativeCompartmentNumber, dependency, seedNumber, burn, ptol, maxiter, iter.CG, working_dir, v)
{
	seedList = list(logL=-Inf)
	while(length(seedList)<=2)
	{
		seed_Evar     = tryCatch({
			seedinit(zvalues = zvalues, 
					covariates = covariates, 
					distances.included = distances.included, 
					A11 = 0.95, A22 = 0.2, 
					alternativeDistribution = alternativeDistribution,
					alternativeCompartmentNumber = alternativeCompartmentNumber, 
					hypothesis = hypothesis, 
					dependency = dependency,
					v = v) }, warning = function(e) {print(e)}, error = function(e) {print(e)})
		seed_Mvar     = tryCatch({
			Maximisation(zvalues = zvalues, 
				covariates = covariates, 
				distances.included = distances.included, 
				Evar = seed_Evar, 
				hypothesis = hypothesis, 
				alternativeDistribution = alternativeDistribution, 
				alternativeCompartmentNumber = alternativeCompartmentNumber, 
				dependency = dependency, 
				iter.CG = iter.CG, 
				ptol = ptol, 
				v = v) }, warning = function(e) {print(e)}, error = function(e) {print(e)})
		if( length(seed_Mvar) != 1)
		{
			rm(seed_Evar)
			seedList = tryCatch({ 
				ExpectationMaximisation(zvalues = zvalues, 
					covariates = covariates, 
					distances.included = distances.included, 
					Mvar = seed_Mvar, 
					hypothesis = hypothesis, 
					alternativeDistribution = alternativeDistribution, 
					alternativeCompartmentNumber = alternativeCompartmentNumber, 
					dependency = dependency, 
					ptol = ptol, 
					maxiter = burn, 
					iter.CG = iter.CG, 
					v=v) }, warning = function(e) {print(e)}, error = function(e) {print(e)})
			rm(seed_Mvar)
		}
		else
		{
			if(v) print("Error in M")
		}
	}
	if(v) cat("seed : ",iter,"/",seedNumber,"   logL :", seedList$logL, '\n')
	return( tryCatch({
		seed_file = paste(paste(working_dir, "seed_",sep="/"),iter,sep="")
		save(seedList, file=seed_file)
		return(list(logL = seedList$logL, file = iter))
	}, error = function(e) {
		if(v) print(paste("error: in seed",iter))
		return(list(logL = -Inf, file = NA))
	}))
}

seedinit = function(zvalues, covariates, distances.included, A11, A22, alternativeDistribution, alternativeCompartmentNumber, hypothesis, dependency, v)
{
	NUM = length(zvalues)
	
	gammA = matrix(rep(0, NUM*2), NUM, 2, byrow=TRUE)
	omega = c()
	if(alternativeDistribution != "kernel" & alternativeDistribution != "kernel.symetric" & alternativeCompartmentNumber > 1)
	{
		omega     = t(rmultinom(NUM, alternativeCompartmentNumber, rep(1/alternativeCompartmentNumber, alternativeCompartmentNumber)))
		omega     = omega[,]/alternativeCompartmentNumber
	}
	
	if(dependency == "none")
	{
		for( i in 1:NUM)
		{
			proba = 1
			while(proba + 1e-4 >= 1 | proba - 1e-4 <= 0)
			{
				proba = rbeta(1, A11, 1-A11)
			}
			gammA[i,1] = proba
		}
		gammA[,2] = 1-gammA[,1]
		return(list(gammA = gammA, omega = omega))
	}
	
	A         = array(0,c(2,2, NUM-1))
	trans.par = array(0,c(2,3+dim(covariates)[2]))
	
	tmp = NA
	while(sum(is.na(tmp))>0 & length(tmp) == 1)
	{
		A_11 = 1
		A_22 = 1
		while(A_11 + 1e-4 >= 1 | A_11 - 1e-4 <= 0)
		{
			A_11 = rbeta(1, A11, 1-A11)
		}
		while(A_22 + 1e-4 >= 1 | A_22 - 1e-4 <= 0)
		{
			A_22 = rbeta(1, A22, 1-A22)
		}
		A[1,1,] = A_11
		A[1,2,] = 1 - A[1,1,1]
		A[2,2,] = A_22
		A[2,1,] = 1 - A[2,2,1]
		tmp = try(inverse.rle( list(values=rep(1:2,NUM) ,lengths=1+rgeom( 2*NUM, rep( c( A[1,2,1], A[2,1,1] ), NUM) )))[1:NUM] - 1)
	}
	gammA[,1] = tmp
	gammA[,2] = 1-gammA[,1]
	
	dgammA    = array(0,c(2,2, NUM-1))
	dgammA[1,1,] = as.numeric(gammA[-dim(gammA)[1],2] == 0 & gammA[-1,1] == 0)
	dgammA[2,1,] = as.numeric(gammA[-dim(gammA)[1],2] == 1 & gammA[-1,1] == 0)
	dgammA[1,2,] = as.numeric(gammA[-dim(gammA)[1],2] == 0 & gammA[-1,2] == 1)
	dgammA[2,2,] = as.numeric(gammA[-dim(gammA)[1],2] == 1 & gammA[-1,2] == 1)
	
	if(length(covariates)>0)
	{
		for(i in 1:2)
		{
			tmp = glm(dgammA[i,2,] ~ covariates[-1,], family=binomial("logit"))
			trans.par[2,1+i]     = tmp$coefficient[1]
			trans.par[2,-c(1:3)] = tmp$coefficient[-1]
		}
		if(distances.included)
		{
			trans.par[2,4] = abs(trans.par[2,4])
		}
	}
	else
	{
		gammA[(gammA[,1] == 1),1] = 0.9
		gammA[(gammA[,1] == 0),1] = 0.1
		gammA[(gammA[,2] == 1),2] = 0.9
		gammA[(gammA[,2] == 0),2] = 0.1
	}
	
#	if(hypothesis != "two.sided")
#	{
#		alternativeCompartmentNumber = alternativeCompartmentNumber*2
#	}
	
	if(length(covariates)>0)
	{
		return(list(gammA = gammA, dgammA = dgammA, omega = omega, trans.par = trans.par))
	}
	else
	{
		return(list(gammA = gammA, dgammA = dgammA, omega = omega))
	}
}

