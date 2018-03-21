\name{Synthetic_purity}
\alias{Synthetic_purity}
\title{Performance evaluation in estimating tumor-purity based on synthetic data.}
\usage{
    Synthetic_purity(p, Rep, n, cluster, sameModuleN, dataseed,  
                    topology, betaModel=TRUE, betaSize=Inf)
}
\arguments{
\item{p}{Number of genes.}
\item{Rep}{Number of replicate data sets.}
\item{n}{Number of observations for each replicate.}
\item{cluster}{Number of disjoint network components.}
\item{sameModuleN}{Number of overlapping modules between tumor and non-tumor components.}
\item{dataseed}{Random seed.}
\item{topology}{Which topology between "Star" and "Power" should be utilized to simulate tumor specific and non-tumor specific networks.}
\item{betaModel}{TRUE if the observed tumor purity (prior) is generated from a beta distribution. FALSE if a logistic regression needs to be utilized.}
\item{betaSize}{Parameter controlling the variance of the distribution of the prior (observed) tumor purity. The mean parameter is equal to the true tumor purity. Default value equal to Inf (i.e., observed tumor purity is equal to the true tumor purity).}
}
\value{List object including:
\item{TSNet.purity}{TSNet's estimate of tumor purity. Vector of dimension (Rep * n).}
\item{priorP}{Prior value of tumor purity. Vector of dimension (Rep * n).}
\item{realP}{True value of tumor purity. Vector of dimension (Rep * n).}
\item{cor.rep}{Correlation between estimated purity and prior value of tumor purity. Vector of dimension rep.}
\item{cor.true}{Correlation between estimated purity and true value of tumor purity. Vector of dimension rep.}
}
\description{Generate synthetic data sets and estimate tumor purity for each of the synthetic data sets.}
\examples{

n=200    # -- number of samples for each replicate
Rep=5   # -- number of replicates
p=500    # -- network dimension
cluster=2   # -- number of disjoint network components
betaSize=1  

#Synthetic_purity(p=p,Rep=Rep,n=n,cluster=cluster,sameModuleN=2,
#dataseed=1,topology="Power", betaSize=betaSize) 

}
\references{
Francesca Petralia, Li Wang , Jie Peng , Arthur Yan,  Jun Zhu and Pei Wang, A new method for constructing tumor specific gene co-expression networks based on samples with tumor purity heterogeneity (Submitted).

Jerome Friedman, Trevor Hastie, and Robert Tibshirani. Sparse inverse covariance estimation with the graphical lasso. \emph{Biostatistics}, 9(3):432-441, 2008.

Patrick Danaher, Pei Wang, and Daniela M Witten. The joint graphical lasso for inverse covariance estimation across multiple classes. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 76(2):373-397, 2014.}
