initialisation = function(parameters)
{
	seedList = mclapply(c(1:parameters[['seed_number']]), FUN = function(x, parameters){ return(runseed(iter = x, parameters = parameters)) }, mc.cores = parameters[['thread_number']], parameters = parameters)
	logL = c()
	for(i in 1:parameters[['seed_number']])
	{
		logL[i] = ifelse(length(seedList[[i]]) > 1, (seedList[[i]])$logL, -Inf)
		if(parameters[['v']]) print(logL[i])
	}
	seedList_ordered = list()
	j = 1
	for(i in order(logL, decreasing = T))
	{
		seedList_ordered[[j]] = seedList[[i]]
		j = j + 1
	}
	
	return(seedList_ordered)
}

loadseed = function(seedList, working_dir)
{
	load(seedList$file)
	return(seedList)
}

runseed = function(iter, parameters)
{
	seedList = list(logL=-Inf, Mvar = -1)
	while(length(seedList) <= 2)
	{
		# the EM can fail, and in this case we set a logL of -Inf for the seed
		seedList = tryCatch({
				ExpectationMaximisation(parameters = parameters, Mvar = Maximisation(parameters = parameters, Evar = seedinit(parameters = parameters)), burn = TRUE)
			}, error = function(e) {
				list(logL=-Inf, Mvar = -1)
			})
	}
	logL = seedList$logL
	# as seed can eat a loot of memory we save them on the disk
	if(parameters[['v']]) print(paste("seed : ",iter,"/",parameters[['seed_number']],"   logL :", logL, sep=""))
	return( tryCatch({
		seed_file = tempfile(pattern = "seed_", tmpdir = tempdir(), fileext = ".RData")
		save(seedList, file=seed_file)
		return(list(logL = logL, file = seed_file))
	}, error = function(e) {
		if(parameters[['v']]) print(paste("error: in seed",iter))
		return(list(logL = -Inf, file = NA))
	}))
}

seedinit = function(parameters)
{
	# gammA is the transition matrix of the markov chain
	# for a NHHM instead of computing the effect of covariate for the 
	# stransition probability between each point, we compute once a transition 
	# matrix for each point.
	gammA = matrix(rep(0, parameters[['NUM']]*2), parameters[['NUM']], 2, byrow=TRUE)
	omega = c()
	
	# omega is the vector of weight of each compartment of the f1 distribution if i's not kernel and as more than one compartment
	if(parameters[['f1']] != "kernel" & parameters[['f1']] != "kernel.symetric" & parameters[['f1_compartiments']] > 1)
	{
		omega     = t(rmultinom(parameters[['NUM']], parameters[['f1_compartiments']], rep(1/parameters[['f1_compartiments']], parameters[['f1_compartiments']])))
		omega     = omega[,]/parameters[['f1_compartiments']]
	}
	
	# in the independent case the transition matrix is always the same and 
	# correspond to a binomial distribution for the two states
	if(parameters[['dependency']] == "none")
	{
		gammA[,1] = ifelse(rbinom(parameters[['NUM']], 1, parameters[['pi_0']]) == 1, 1, 0)
		if(length(gammA[gammA[,1] == 0.1,1] ) <= 1)
			gammA[1:10,1] = 0
		gammA[,2] = 1-gammA[,1]
		return(list(gammA = gammA, omega = omega))
	}
	
	if(length(parameters[['covariates']])>0)
	{
		for(i in 1:2)
		{
			tmp = list(converged = FALSE)
			# in the case of a NHMM the initialiation may not converge (glm)
			# thus we loop until it does
			while(!tmp$converged)
			{
				tmp = initATransParGammADgammA(parameters[['NUM']], parameters[['zvalues']], parameters[['covariates']], parameters[['A11']], parameters[['A22']])
				A = tmp$A
				trans.par = tmp$trans.par
				gammA = tmp$gammA
				dgammA = tmp$dgammA
				# we fit the transition parameter to the HMM with a general 
				# linear model which may not converge
				tmp = tryCatch({
					glm(dgammA[i,2,] ~ parameters[['covariates']][-1,], family=binomial("logit"))
				}, warning = function(e) {list(converged = FALSE)}, error = function(e) {list(converged = FALSE)})
			}
			
			trans.par[2,1+i]     = tmp$coefficient[1]
			trans.par[2,-c(1:3)] = tmp$coefficient[-1]
		}
		# we force the transition parameter associated with the distance to
		# be positif
		if(parameters[['distances_included']])
		{
			trans.par[2,4] = abs(trans.par[2,4])
		}
		return(list(gammA = gammA, dgammA = dgammA, omega = omega, trans.par = trans.par))
	}
	else
	{
		tmp = initATransParGammADgammA(parameters[['NUM']], parameters[['zvalues']], parameters[['covariates']], parameters[['A11']], parameters[['A22']])
		A = tmp$A
		gammA = tmp$gammA
		dgammA = tmp$dgammA
		gammA[(gammA[,1] == 1),1] = 1
		gammA[(gammA[,1] == 0),1] = 0
		gammA[(gammA[,2] == 1),2] = 1
		gammA[(gammA[,2] == 0),2] = 0
		return(list(gammA = gammA, dgammA = dgammA, omega = omega))
	}
}

# function to generate HMM state with random transition probability
initATransParGammADgammA = function(NUM, zvalues, covariates, A11, A22)
{
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
		# the fastess way to generate HMM states in R
		tmp = try(inverse.rle( list(values=rep(1:2,NUM) ,lengths=1+rgeom( 2*NUM, rep( c( A[1,2,1], A[2,1,1] ), NUM) )))[1:NUM] - 1)
	}
	gammA = matrix(rep(0, NUM*2), NUM, 2, byrow=TRUE)
	gammA[,1] = tmp
	gammA[zvalues == 0, 1] = 1
	gammA[,2] = 1-gammA[,1]
	
	dgammA    = array(0,c(2,2, NUM-1))
	dgammA[1,1,] = as.numeric(gammA[-dim(gammA)[1],2] == 0 & gammA[-1,1] == 0)
	dgammA[2,1,] = as.numeric(gammA[-dim(gammA)[1],2] == 1 & gammA[-1,1] == 0)
	dgammA[1,2,] = as.numeric(gammA[-dim(gammA)[1],2] == 0 & gammA[-1,2] == 1)
	dgammA[2,2,] = as.numeric(gammA[-dim(gammA)[1],2] == 1 & gammA[-1,2] == 1)
	
	return(list(A = A, trans.par = trans.par, gammA = gammA, dgammA = dgammA))
}
