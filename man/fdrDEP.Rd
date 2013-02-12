\name{fdrDEP}
\alias{fdrDEP}
\title{
FDR control based for multiple hypothesis tests under dependence
}
\description{
Main function which fits the NHMM and HMM based fdr control. It offers numerous options for FDR control as described in Kuan et al. (2011) and implementation for one-sided hypothesis case.

}
\usage{
fdrDEP(pvalues = x, covariates = NULL, distances = NULL, observerdValues = NULL, hypothesis = "two.sided", threshold = 0, alternativeDistribution = 'kernel', alternativeCompartmentNumber = 2, dependency = 'NHMM', seedNumber = 20, burn = 20, ptol = 1e-3, core = 2, maxiter=1000, iter.CG = 1000, v = F, transformed = F)
}
\arguments{
  \item{pvalues}{
A vector of the observed pvalues data. It is assumed that the pvalues has been ordered by time or position.
}
\item{covariates}{
A matrix of covariates EXCLUDING spacing/distance between probes. Each column corresponds to one covariate. If there is only one covariate
}
\item{distances}{
A vector of spacing/distance between data
}
\item{hypothesis}{
"\code{c("two.sided", "one.sided"}" apply the correct z-values transform and LIS model for one-sided hypothesis testing or two sided hypothesis testing
}
\item{alternativeDistribution}{
Type of estimating for alternative hypothesis f1(x). Available choices are "\code{mixnormal}" for gaussian mixtures or "\code{kernel}" or "\code{kernel.symetric}" for non-parametric kernel density estimates.
}
\item{alternativeCompartmentNumber}{
Number of mixture component for \code{alttype = "mixnormal"}. If  \code{alttype = "kernel"}, this is irrelevant.
}
\item{maxiter}{
Maximum iterations in the EM algorithm to speed up computation.
}
\item{nulltype}{
Type of null hypothesis assumed in estimating f0(z), for use in the fdr calculations. \code{0} is the theoretical null N(0,1), \code{1} is maximum likelihood estimation, \code{2} is central matching estimation. This method is imported from R package \code{locfdr}.\\ NOTE: Recommended to use nulltype = 0, i.e., theoretical null to avoid double correction of correlation structure.
}
\item{dependency}{
Types of dependence structure. Available choices are "\code{none}", "\code{HMM}" or "\code{NHMM}".
}
\item{transformed}{
\code{TRUE} if the p-values provided are allready transformed : skip the z-values transform
}
\item{seedNumber}{
number of seed for the initialisation step
}
\item{burn}{
number of iterations of the EM algorithm for each seed
}
\item{core}{
number of core to use for the seed generation
}
\item{maxiter}{
maximum number of iterations for the EM algorithm
}
\item{iter.CG}{
maximum number of iterations for the CG algorithm
}
\item{v}{
\code{TRUE} verbose mode
}
}

\value{

\item{LIS}{Local index of significance (analog of p-values). This will be used for FDR control in \code{\link{LIS.adjust}}.}
\item{BIC}{Bayesian information criterion (BIC).}
\item{pii}{Initial probabilities.}
\item{A}{Transition probability matrix.}
\item{f0}{Null hypothesis.}
\item{f1}{Alternative hypothesis.}
\item{logL}{Log likelihood.}
\item{trans.par}{Transition parameters for State 1.}
}
\references{
P.F. Kuan and D.Y. Chiang (2012). Integrating Prior Knowledge in Multiple Testing Under Dependence with Applications in Detecting Differential DNA Methylation. Biometrics, doi: 10.1111/j.1541-0420.2011.01730.x.

B. Efron (2004). Large-scale simultaneous hypothesis testing: the choice of a null hypothesis. Journal of the American Stat. Assoc. 99, 96-104.

W. Sun and T. Cai (2009). Large-scale multiple testing under dependence. J. R.Stat. Soc B. 71, 393-424.
}
\author{Laurent MODOLO}
\seealso{ \code{\link{plis}} }
\examples{

}

