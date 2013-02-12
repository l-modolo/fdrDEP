\name{plis}
\alias{plis}
\title{
Determine statistically significant tests at a user-specified FDR level.
}
\description{
Compute adjusted local index of significance (LIS) and determine the tests which are statistically significant at a user-specified FDR level.

}
\usage{
plis(LIS, fdr = 0.001, adjust = FALSE)
}
\arguments{
  \item{LIS}{
LIS scores from \code{\link{fdrDEP}}.
}
  \item{fdr}{
Desirable FDR level. Usually 0.05. 
}
  \item{adjust}{
Logical value. \code{TRUE} if adjusted LIS is to be computed. See \code{details}.
}
}
\details{
If \code{adjust = TRUE}, the adjusted LIS will be computed. This will be useful if the user wants to control at different FDR level. The significant tests corresponds to those adjusted LIS which are less than the FDR level.
}
\value{
If \code{adjust = TRUE}
\item{states}{A binary valued vector. Takes value 1 if the test or observation is declared to be significant.}
\item{aLIS}{Adjusted LIS.}
else
A binary valued vector. Takes value 1 if the test or observation is declared to be significant.
}
\references{
W. Sun and T. Cai (2009). Large-scale multiple testing under dependence. J. R.Stat. Soc B. 71, 393-424.
}
\author{Laurent Modolo}
\seealso{ \code{\link{fdrDEP}}}
\examples{

}

