\name{Synthetic_data}
\alias{Synthetic_data}
\title{Generate independent datasets from a mixed expression model.}
\usage{
Synthetic_data(p,rep,n,cluster,parameterseed,dataseed,meanY=NULL,
               meanZ=NULL,sdY=NULL,sdZ=NULL,P.v=NULL,meanP=NULL,
               varP=NULL,sameModuleN=0,topology="Power",
               betaModel=TRUE,betaSize=Inf)
}
\arguments{
\item{p}{Number of genes (i.e., number of nodes in the network).}
\item{rep}{Number of independent data sets (i.e., replicates).}
\item{n}{Number of observations for each replicate.}
\item{cluster}{Number of disjoint components in the network.}
\item{parameterseed}{Random seed for parameter sampling.}
\item{dataseed}{Random seed for observations sampling.}
\item{meanY}{\code{p} dimensional vector containing the mean parameter of the Gaussian distribution of the gene expression in tumor component. When \code{NULL}, \code{meanY=0}}
\item{meanZ}{p dimensional vector containing the mean parameter of the Gaussian distribution of the gene expression in non-tumor component. When \code{NULL}, \code{meanZ=0}.}
\item{sdY}{p dimensional vector containing the variance parameter of the Gaussian distribution of the gene expression in tumor component. When \code{NULL}, \code{sdY=1}.}
\item{sdZ}{p dimensional vector containing the variance parameter of the Gaussian distribution of the gene expression in non-tumor component. When \code{NULL}, \code{sdZ=1}.}
\item{P.v}{True level of tumor-purity. When \code{NULL}, the value of tumor purity is simulated from a Beta distribution with mean meanP and variance varP.}
\item{meanP}{Mean of the \code{Beta} distribution utilized to simulate tumor purity. When \code{NULL}, meanP=0.5.}
\item{varP}{Variance of the \code{Beta} distribution utilized to simulate tumor purity. When \code{NULL}, varP=0.04.}
\item{sameModuleN}{Number of clusters shared between the networks of tumor and non-tumor components.}
\item{topology}{Which topology between Star and Power should be utilized to simulate tumor-specific and non-tumor specific networks. Default value is Power.}
\item{betaModel}{\code{TRUE} if the observed tumor purity (prior) is drawn from a beta distribution. \code{FALSE} if a logistic regression is utilized.}
\item{betaSize}{Parameter controlling the variance of the distribution of the prior (observed) tumor purity. The mean parameter is equal to the true tumor purity. Default value equal to Inf (i.e., observed tumor purity is equal to the true tumor purity)}
}

\value{List object containing the following values:
\item{X}{(p x (n * rep)) matrix of gene expression from mixed model.}
\item{Uy}{(p x 1) vector containing mean value for each gene in tumor component.}
\item{Uz}{(p x 1) vector containing mean value for each gene in non-tumor component.}
\item{Sigmay}{(p x p) matrix containing covarian matrix of tumor component.}
\item{Sigmaz}{(p x p) matrix containing covarian matrix of non-tumor component.}
\item{P}{(n x 1) vector containing observed values of tumor purity.}
\item{Y}{(p x (n * rep)) matrix containing expression profile of tumor component.}
\item{Z}{(p x (n * rep)) matrix containing expression profile of non-tumor component.}
\item{realP}{((n * rep) x 1) vector containing true values of tumor purity.}
}
\description{This function is utilized to sample mixed expression profiles.}
\examples{

n=100    # -- number of samples for each replicate
rep=5   # -- number of replicates
p=500    # -- network dimension
betaSize=1

#Data<-Synthetic_data(p=p,rep=rep,n=n,cluster=2,parameterseed=1,
#              dataseed=1,meanY=NULL, meanZ=NULL,sdY=NULL,sdZ=NULL, 
#              P.v=NULL,meanP=NULL,varP=NULL, sameModuleN=0,
#              betaModel=TRUE,topology="Power",betaSize=betaSize)

}
\references{
Francesca Petralia, Li Wang , Jie Peng , Arthur Yan,  Jun Zhu and Pei Wang, A new method for constructing tumor specific gene co-expression networks based on samples with tumor purity heterogeneity (Submitted).

Jerome Friedman, Trevor Hastie, and Robert Tibshirani. Sparse inverse covariance estimation with the graphical lasso. \emph{Biostatistics}, 9(3):432-441, 2008.

Patrick Danaher, Pei Wang, and Daniela M Witten. The joint graphical lasso for inverse covariance estimation across multiple classes. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 76(2):373-397, 2014.}
