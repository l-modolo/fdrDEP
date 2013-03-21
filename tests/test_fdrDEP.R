#!/usr/bin/Rscript
options(width = 320)

system("R CMD build /home/vmiele/Collaborations/Modolo/fdrDEP && mv fdrDEP_1.0.0.tar.gz /tmp")
try(remove.packages("fdrDEP"))
install.packages("/tmp/fdrDEP_1.0.0.tar.gz")
library(fdrDEP)

# simul data
compute.A.nhmm = function (Z, trans.par1, trans.par2, dist.included = TRUE) 
{
	p = dim(Z)[2]
	pii = rep(0, 2)
	A = array(0, c(2, 2, dim(Z)[1] - 1))
	pii[1] = sum(exp(trans.par1[1] + sum(trans.par1[-c(1:3)] * 
	Z[1, ])))/(sum(exp(trans.par1[1] + sum(trans.par1[-c(1:3)] * 
	Z[1, ]))) + sum(exp(trans.par2[1] + sum(trans.par2[-c(1:3)] * 
	Z[1, ]))))
	pii[2] = 1 - pii[1]
	tmp12 = t(trans.par2[-c(1:3)] * t(Z[-1, ]))
	tmp22 = t(trans.par2[-c(1:3)] * t(Z[-1, ]))
	if (dist.included == TRUE) {
		tmp12 = t(c(trans.par2[4], trans.par2[-c(1:4)]) * t(Z[-1, 
		    ]))
		tmp22 = t(c(-trans.par2[4], trans.par2[-c(1:4)]) * t(Z[-1, 
		    ]))
	}
	denom12 = (denom22 = rep(0, dim(Z)[1] - 1))
	for (i in 1:p) {
		denom12 = denom12 + tmp12[, i]
		denom22 = denom22 + tmp22[, i]
	}
	num11 = exp(trans.par1[2])
	num21 = exp(trans.par1[3])
	num12 = exp(trans.par2[2] + denom12)
	num22 = exp(trans.par2[3] + denom22)
	A[1, 1, ] = num11/(num11 + num12)
	A[1, 2, ] = num12/(num11 + num12)
	A[2, 1, ] = num21/(num21 + num22)
	A[2, 2, ] = num22/(num21 + num22)
	return(list(A = A, pii = pii))
}

simdata.nhmm = function (NUM, pii, A, f0, pc, f1) 
{
	theta = rep(0, NUM)
	x = rep(0, NUM)
	nc = length(pc)
	if (nc == 1) {
		data = simdata1.nhmm(NUM, pii, A, f0, f1)
	}
	else {
		theta[1] = rbinom(1, 1, pii[2])
		for (i in 2:NUM) {
			if (theta[i - 1] == 0) 
			theta[i] = rbinom(1, 1, A[1, 2, i - 1])
			else theta[i] = rbinom(1, 1, A[2, 2, i - 1])
		}
		for (i in 1:NUM) {
			if (theta[i] == 0) {
				mu0 = f0[1]
				sd0 = f0[2]
				x[i] = rnorm(1, mean = mu0, sd = sd0)
			}
			else {
				c = sample(1:nc, 1, prob = pc)
				mu1 = f1[c, 1]
				sd1 = f1[c, 2]
				x[i] = rnorm(1, mean = mu1, sd = sd1)
			}
		}
		data = list(s = theta, o = x)
	}
	return(data)
}

# display data partition

f0f1 = function (zval, LIS, k=F, title="")
{
#	plot(density(zval[zval != 0]))
	hist(zval[zval != 0], nclass=sqrt(length(zval)), main=title, freq=F, xlim=c(-10,10), ylim=c(0,1), xlab="z-value", cex.main = 3, cex.lab = 2, cex.axis = 2.5)
	
	legend("topright", legend=c("F0", "F1", "F0 + F1"), 
	col=c("blue", "red", "green"), lty=1, lwd=4, cex=3, bty="n")
	
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
}

display_data = function(x, fit, k = F)
{
	states = plis(fit$LIS, fdr=0.1) # fdr = 0.1 -> mean(emp.fdr) = 0.1
	sig.test = which(states == 1)
	TP = length(which(states == 1 & theta1 == 1))
	FP = length(which(states == 1 & theta1 == 0))
	emp.fdr = FP/length(sig.test)
	cat("fdr : ",emp.fdr, " TP : ", TP, " FP : ", FP,'\n')
	f0f1(x, fit, k=k)
}


################################ simulated data ################################

# size of the data
NUM1 = 2000 # 2000 = quick; 20000 = slow

SIMUL=FALSE
if (SIMUL){
	Z = rnorm(NUM1)
	Z = matrix(Z,ncol=1)
	Z = apply(Z,2,scale)
	Dist = runif(NUM1,min=0,max=1)
	log2Dist = log2(Dist+2)
	covariate = cbind(log2Dist,Z)
	trans.par1.true = c(0,0,0,0,0)
	trans.par2.true = c(0.6,-3,-1,2,-1)
	
# compute the transition probabilities
# given transition parameters
	A.true = compute.A.nhmm(covariate, trans.par1.true, trans.par2.true, dist.included=TRUE)$A
	pii.true = compute.A.nhmm(covariate, trans.par1.true, trans.par2.true, dist.included=TRUE)$pii
	
# null distribution
	f0 = c(0, 1)
# alternative distribution
	f1 = matrix(c(3, -2,2, -4, 1, 1,1.2,1.4),ncol=2)
	### simulate NHMM data
	simdat = simdata.nhmm(NUM1, pii.true, A.true, f0, c(0.3,0.25,0.2,0.25), f1)
	### simulated observed values
	x = simdat$o
	### simulated unobserved true states
	theta1 = simdat$s
	table(theta1)
	
	save.image(file = paste("/home/vmiele/Collaborations/Modolo/fdrDEP/tests/databenchs/bench", NUM1, ".RData", sep=""))
} else{
	load(paste("/home/vmiele/Collaborations/Modolo/fdrDEP/tests/databenchs/bench", NUM1, ".RData", sep=""))
}

# Model fitting
# complexity : fit.none.k < fit.hmm.k < fit.nhmm.k
# complexity : fit.none.mix < fit.hmm.mix < fit.nhmm.mix
# complexity CPU : kernel < mixnormal
# complexity memory : mixnormal < kernel
NBCORES=1
Rprof("/home/vmiele/Collaborations/Modolo/fdrDEP/tests/profiling/fitnonek.out")
fit.none.k = fdrDEP(pvalues = x, covariates = NULL, distances = NULL, observerdValues = Null, hypothesis = "two.sided", threshold = NULL, alternativeDistribution = 'kernel', alternativeCompartmentNumber = 2, dependency = 'none', seedNumber = 20, burn = 20, ptol = 1e-3, core = NBCORES, maxiter=1000, iter.CG = 1000, v = T, trans = T)
Rprof()
summaryRprof("/home/vmiele/Collaborations/Modolo/fdrDEP/tests/profiling/fitnonek.out")

Rprof("/home/vmiele/Collaborations/Modolo/fdrDEP/tests/profiling/fithmmk.out")
fit.hmm.k = fdrDEP(pvalues = x, covariates = NULL, distances = NULL, observerdValues = Null, hypothesis = "two.sided", threshold = NULL, alternativeDistribution = 'kernel', alternativeCompartmentNumber = 2, dependency = 'HMM', seedNumber = 20, burn = 20, ptol = 1e-3, core = NBCORES, maxiter=1000, iter.CG = 1000, v = T, trans = T)
Rprof()
summaryRprof("/home/vmiele/Collaborations/Modolo/fdrDEP/tests/profiling/fithmmk.out")

Rprof("/home/vmiele/Collaborations/Modolo/fdrDEP/tests/profiling/fitnhmmk.out")
fit.nhmm.k = fdrDEP(pvalues = x, covariates = Z, distances = Dist, observerdValues = Null, hypothesis = "two.sided", threshold = NULL, alternativeDistribution = 'kernel', alternativeCompartmentNumber = 2, dependency = 'NHMM', seedNumber = 20, burn = 20, ptol = 1e-3, core = NBCORES, maxiter=1000, iter.CG = 1000, v = T, trans = T)
Rprof()
summaryRprof("/home/vmiele/Collaborations/Modolo/fdrDEP/tests/profiling/fitnhmmk.out")

Rprof("/home/vmiele/Collaborations/Modolo/fdrDEP/tests/profiling/fitnonemix.out")
fit.none.mix = fdrDEP(pvalues = x, covariates = NULL, distances = NULL, observerdValues = Null, hypothesis = "two.sided", threshold = NULL, alternativeDistribution = 'mixnormal', alternativeCompartmentNumber = 2, dependency = 'none', seedNumber = 20, burn = 20, ptol = 1e-3, core = NBCORES, maxiter=1000, iter.CG = 1000, v = T, trans = T)
Rprof()
summaryRprof("/home/vmiele/Collaborations/Modolo/fdrDEP/tests/profiling/fitnonemix.out")

Rprof("/home/vmiele/Collaborations/Modolo/fdrDEP/tests/profiling/fithmmmix.out")
fit.hmm.mix = fdrDEP(pvalues = x, covariates = NULL, distances = NULL, observerdValues = Null, hypothesis = "two.sided", threshold = NULL, alternativeDistribution = 'mixnormal', alternativeCompartmentNumber = 2, dependency = 'HMM', seedNumber = 20, burn = 20, ptol = 1e-3, core = NBCORES, maxiter=1000, iter.CG = 1000, v = T, trans = T)
Rprof()
summaryRprof("/home/vmiele/Collaborations/Modolo/fdrDEP/tests/profiling/fithmmmix.out")

Rprof("/home/vmiele/Collaborations/Modolo/fdrDEP/tests/profiling/fitnhmmmix.out")
fit.nhmm.mix = fdrDEP(pvalues = x, covariates = Z, distances = Dist, observerdValues = Null, hypothesis = "two.sided", threshold = NULL, alternativeDistribution = 'mixnormal', alternativeCompartmentNumber = 2, dependency = 'NHMM', seedNumber = 20, burn = 20, ptol = 1e-3, core = NBCORES, maxiter=1000, iter.CG = 1000, v = T, trans = T)
Rprof()
summaryRprof("/home/vmiele/Collaborations/Modolo/fdrDEP/tests/profiling/fitnhmmmix.out")


############################# display results ##################################
par(mfrow=c(3,2))
# empirical fdr and graph depdency = none, f1 = kernel
display_data(x, fit.none.k, k=T)
# empirical fdr and graph depdency = hmm, f1 = kernel
display_data(x, fit.hmm.k, k=T)
# empirical fdr and graph depdency = nhmm, f1 = kernel
display_data(x, fit.nhmm.k, k=T)
# empirical fdr and graph depdency = none, f1 = mixnormal
display_data(x, fit.none.mix, k=F)
# empirical fdr and graph depdency = hmm, f1 = mixnormal
display_data(x, fit.hmm.mix, k=F)
# empirical fdr and graph depdency = nhmm, f1 = mixnormal
display_data(x, fit.nhmm.mix, k=F)

