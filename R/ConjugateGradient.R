ComputeCG = function(covariates, distances.included, dgammA, gammA, trans.par, iter.CG, ptol, v)
{
#if(v) print("CG step")
# cat("Saving input of ComputeCG\n")
#save(covariates, distances.included, dgammA, gammA, trans.par, iter.CG, ptol, v, file=paste("/home/vmiele/Collaborations/Modolo/fdrDEP/tests/databenchs/inputComputeCG",dim(covariates)[1],".RData",sep=""))
#stop("Done")

gradient.old = ComputeGradient.C(covariates, distances.included, dgammA, gammA, trans.par, v = v)
trans.par.old = trans.par
phi = rbind(rep(0, length(gradient.old)),-gradient.old)
difference = 1
niter = 0
while(difference > ptol & niter < iter.CG)
{
trans.par.old = trans.par
niter = niter + 1
tmp = tryCatch({ LineSearch(covariates, distances.included, dgammA, gammA, trans.par, phi, iter.CG, ptol, v = v) }, warning = function(e) {print(e)}, error = function(e) {print(e)})
if(is.na(tmp$nu))
{
break
}
trans.par = tmp$trans.par
gradient.new = ComputeGradient.C(covariates, distances.included, dgammA, gammA, trans.par, v = v)
PR = sum((gradient.new - gradient.old)*gradient.new) / sum(gradient.old^2)
if(is.nan(PR) | PR < 0)
{
PR = 0
}

phi = rbind(rep(0, length(gradient.new)), -gradient.new + PR * phi[2,])
gradient.old = gradient.new
if(distances.included & trans.par[2,4]<0 )
{
trans.par[2,4] = abs(trans.par[2,4])
}
difference = max(abs(trans.par[2,-1] - trans.par.old[2,-1]))
}
if(!is.na(tmp$nu))
{
tmp = pii.A.C(covariates, distances.included, trans.par, v = v)
return(list(pii = tmp$pii, A = tmp$A, trans.par = trans.par))
}else
{
return(-1)
}
}

LineSearch = function(covariates, distances.included, dgammA, gammA, trans.par, phi, iter.CG, ptol, v)
{
N = dim(covariates)[1] - 1
nu = 0
nu.old = 1
niter = 0
difference = 1
phi[1,] = 0

trans.par.new = trans.par

fixe_1 = phi[,1] + sum( phi[,-c(1:3)] * covariates[1,] )
if(distances.included)
{
fixe_1 = phi[,1] + sum( c(-phi[,4],phi[,-c(1:4)]) * covariates[1,] )
}
fixe_2 = array(0,c(2,2,N))
for(i in 1:2)
{
for(j in 1:2)
{
if(distances.included & i == 2)
fixe_2[i,j,] = ( phi[j,i+1] + rowSums( t(c(-phi[j,4],phi[j,-c(1:4)]) * t(covariates[-1,])) ) )
else
fixe_2[i,j,] = ( phi[j,i+1] + rowSums( t(phi[j,-c(1:3)] * t(covariates[-1,])) ) )
}
}

while(difference > ptol & niter < iter.CG)
{
niter = niter + 1
trans.par.new[1,] = trans.par[1,] + nu * phi[1,]
trans.par.new[2,] = trans.par[2,] + nu * phi[2,]
trans.prob = pii.A.C(covariates, distances.included, trans.par.new, v = v)
dQ = sum( fixe_1 * (gammA[1,] - trans.prob$pii) )
dQ2 = -sum( fixe_1^2 * trans.prob$pii * (1 - trans.prob$pii) )

for(i in 1:2)
{
for(j in 1:2)
{
dQ = dQ + sum( fixe_2[i,j,] * ( dgammA[i,j,] - ( gammA[-dim(gammA)[1],i] * trans.prob$A[i,j,] ) ) )
dQ2 = dQ2 - sum( fixe_2[i,j,]^2 * trans.prob$A[i,j,] * ( 1 - trans.prob$A[i,j,] ) * gammA[-dim(gammA)[1], i] )
}
}

nu = nu - dQ / dQ2
difference = abs(dQ/ dQ2)
if(is.na(difference) | niter > 100){
nu = NaN
break
}
}
if(is.nan(dQ2))
{
if(v) print("Error in line search")
}
trans.par.new[1,] = trans.par[1,] + nu * phi[1,]
trans.par.new[2,] = trans.par[2,] + nu * phi[2,]
return(list(nu = nu, trans.par = trans.par.new))
}

LineSearch.C = function(covariates, distances.included, dgammA, gammA, trans.par, phi, iter.CG, ptol, v)
{
nu = 0
trans.par.new = trans.par
res <- .C('C_LineSearch', as.integer(dim(covariates)[1]), as.integer(dim(covariates)[2]),
as.double(covariates),
as.integer(distances.included),
as.double(dgammA), as.double(gammA),
as.double(trans.par), as.double(phi), as.integer(iter.CG),
as.double(ptol), as.integer(v),
resnu = as.double(nu), restrans.par.new = as.double(trans.par.new))
return(list(nu = res$resnu, trans.par = res$restrans.par.new))
}

ComputeGradient = function(covariates, distances.included, dgammA, gammA, trans.par, v)
{
gradient <- rep(0,3+dim(covariates)[2])

tmp.trans.prob = pii.A.C(covariates, distances.included, trans.par, v = v)

pii <- tmp.trans.prob$pii
A <- tmp.trans.prob$A

gradient[1] <- gammA[1,2] - pii[2]	# G1
gradient[2] <- sum(dgammA[1,2,] - gammA[-dim(gammA)[1],1]*A[1,2,])	# G2
gradient[3] <- sum(dgammA[2,2,] - gammA[-dim(gammA)[1],2]*A[2,2,])	# G3

tmp <- rep(0,dim(dgammA)[3])

for(i in 1:2){
tmp <- tmp + gammA[-dim(covariates)[1],i]*A[i,2,] # G422
}
gradient[-c(1:3)] <- (gammA[1,2] - pii[2])*covariates[1,] +
apply(matrix((gammA[-1,2] - tmp)*covariates[-1,],ncol=dim(covariates)[2]),2,sum) # G41 + G42

if(distances.included==TRUE){
gradient[4] <- (gammA[1,2] - pii[2])*covariates[1,1] +
sum(dgammA[1,2,]*covariates[-1,1]) -
sum(dgammA[2,2,]*covariates[-1,1]) -
sum(gammA[-dim(covariates)[1],1]*A[1,2,]*covariates[-1,1]) +
sum(gammA[-dim(covariates)[1],2]*A[2,2,]*covariates[-1,1])
}

return(-gradient)
}

ComputeGradient.C = function(covariates, distances.included, dgammA, gammA, trans.par, v)
{
gradient <- rep(0,3+dim(covariates)[2])
pii <- rep(0,2)
A <- array(0,c(2,2,dim(covariates)[1]-1))
res <- .C('C_ComputeGradient', as.integer(dim(covariates)[1]), as.integer(dim(covariates)[2]),
as.double(covariates),
as.integer(distances.included),
as.double(dgammA), as.double(gammA),
as.double(trans.par), as.integer(v),
as.double(pii), as.double(A),
resgradient = as.double(gradient))
return(-res$resgradient)
}

pii.A = function(covariates, distances.included, trans.par, v)
{
# covariates is a matrix of size m x p
p <- dim(covariates)[2]
trans.par1 = trans.par[1,]
trans.par2 = trans.par[2,]

pii <- rep(0,2)
A <- array(0,c(2,2,dim(covariates)[1]-1))

pii[1] <- sum(exp(trans.par1[1] + sum(trans.par1[-c(1:3)]*covariates[1,])))/
( sum(exp(trans.par1[1] + sum(trans.par1[-c(1:3)]*covariates[1,])))
+ sum(exp(trans.par2[1] + sum(trans.par2[-c(1:3)]*covariates[1,]))) )
pii[2] = 1 - pii[1]

tmp12 <- t(trans.par2[-c(1:3)]*t(covariates[-1,]))
tmp22 <- t(trans.par2[-c(1:3)]*t(covariates[-1,]))

if(distances.included==TRUE)
{
tmp12 <- t(c(trans.par2[4],trans.par2[-c(1:4)])*t(covariates[-1,]))
tmp22 <- t(c(-trans.par2[4],trans.par2[-c(1:4)])*t(covariates[-1,]))
}

denom12 = denom22 <- rep(0,dim(covariates)[1]-1)
for(i in 1:p)
{
denom12 <- denom12 + tmp12[,i]
denom22 <- denom22 + tmp22[,i]
}

num11 <- exp(trans.par1[2])
num21 <- exp(trans.par1[3])
num12 <- exp(trans.par2[2] + denom12)
num22 <- exp(trans.par2[3] + denom22)


A[1,1,] <- num11/(num11 + num12)
A[1,2,] <- num12/(num11 + num12)
A[2,1,] <- num21/(num21 + num22)
A[2,2,] <- num22/(num21 + num22)

return(list(A = A, pii = pii))
}

pii.A.C = function(covariates, distances.included, trans.par, v)
{
pii <- rep(0,2)
A <- array(0,c(2,2,dim(covariates)[1]-1))
res <- .C('C_piiA', as.integer(dim(covariates)[1]), as.integer(dim(covariates)[2]), as.double(covariates),
as.integer(distances.included), as.double(trans.par), as.integer(v),
respii = as.double(pii), resA = as.double(A))
return(list(A = array(res$resA, dim(A)), pii = res$respii))
}

pii.A.Call = function(covariates, distances.included, trans.par, v)
{
pii <- rep(0,2)
A <- array(0,c(2,2,dim(covariates)[1]-1))
.Call('Call_piiA',A);
return(list(A = array(res$resA, dim(A)), pii = res$respii))
}


