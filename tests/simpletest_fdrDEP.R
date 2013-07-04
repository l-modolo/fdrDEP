#!/usr/bin/Rscript
options(width = 320)

system("R CMD build .. && mv fdrDEP_1.0.0.tar.gz /tmp")
try(remove.packages("fdrDEP"))
install.packages("/tmp/fdrDEP_1.0.0.tar.gz")
library(fdrDEP)

################################ simulated data ################################

# size of the data
#NUM1 = 2000 # 2000 = quick; 20000 = slow
#load(paste("/home/vmiele/Collaborations/Modolo/fdrDEP/tests/databenchs/bench", NUM1, ".RData", sep=""))
#NBCORES=1
#fit.nhmm.k = fdrDEP(pvalues = x, covariates = Z, distances = Dist, observerdValues = Null, hypothesis = "two.sided", threshold = NULL, alternativeDistribution = 'kernel', alternativeCompartmentNumber = 2, dependency = 'NHMM', seedNumber = 1, burn = 20, ptol = 1e-3, core = NBCORES, maxiter=1000, iter.CG = 1000, v = T, trans = T)


################################ focus on ConjugateGradient ################################

NUM1 = 2000 # 2000 = quick; 20000 = slow
load(paste("databenchs/inputComputeCG",NUM1,".RData",sep=""))

################################
distances.included = FALSE

tmp.trans.prob = pii.A(covariates, distances.included, trans.par, v)
cat(tmp.trans.prob$pii,"\n")
cat(tmp.trans.prob$A[1,1,1], tmp.trans.prob$A[1,2,1], tmp.trans.prob$A[2,1,1], tmp.trans.prob$A[2,2,1], tmp.trans.prob$A[1,1,3], tmp.trans.prob$A[1,1,4], tmp.trans.prob$A[1,2,4],"\n")

tmp.trans.prob.C = pii.A.C(covariates, distances.included, trans.par, v)
cat(tmp.trans.prob.C$pii,"\n")
cat(tmp.trans.prob.C$A[1,1,1], tmp.trans.prob.C$A[1,2,1], tmp.trans.prob.C$A[2,1,1], tmp.trans.prob.C$A[2,2,1], tmp.trans.prob.C$A[1,1,3], tmp.trans.prob.C$A[1,1,4], tmp.trans.prob.C$A[1,2,4],"\n")

################################
distances.included = FALSE

tmp.gr = ComputeGradient(covariates, distances.included, dgammA, gammA, trans.par, v)
cat(tmp.gr,"\n")

tmp.gr.C = ComputeGradient.C(covariates, distances.included, dgammA, gammA, trans.par, v)
cat(tmp.gr.C,"\n")

################################
distances.included = FALSE

gradient.old = ComputeGradient.C(covariates, distances.included, dgammA, gammA, trans.par, v)
phi = rbind(rep(0, length(gradient.old)),-gradient.old)

tmp.ls = LineSearch(covariates, distances.included, dgammA, gammA, trans.par, phi, iter.CG, ptol, v)
cat(tmp.ls$nu,"\n",tmp.ls$trans.par,"\n")

tmp.ls.C = LineSearch.C(covariates, distances.included, dgammA, gammA, trans.par, phi, iter.CG, ptol, v)
cat(tmp.ls.C$nu,"\n",tmp.ls.C$trans.par,"\n")
