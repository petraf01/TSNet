\name{deNet_purity}
\alias{deNet_purity}
\title{Function to estimate tumor purity from mixed expression data.}
\usage{
deNet_purity(exprM,purity)
}
\arguments{
\item{exprM}{A \code{(n x p)} matrix of gene expression with p being the number of genes and n the number of samples.}
\item{purity}{An \code{n}-dimensional vector containing tumor-purity for each sample.}}
\value{
\item{Epurity}{Estimated tumor purity via \code{TSNet}.} 
\item{paraM}{\code{(p x 2)} matrix containing mean parameter of tumor and non-tumor components for each gene. First column contains mean of tumor component, while the second column the mean of non-tumor component.}
\item{convergence}{value indicating whether or not the algorithm converged}
\item{time}{total running time.}
}
\description{TSNet function to estimate tumor purity from mixed expression profiles.}
\examples{

set.seed(1)
p<-500; n=200;
cluster=4
betaSize=1

# --- Sample Data
Data<-Synthetic_data(p=p,rep=1,n=n,cluster=cluster,parameterseed=4,
              dataseed=1,meanY=NULL, meanZ=NULL,sdY=NULL,sdZ=NULL, 
              P.v=NULL,meanP=NULL,varP=NULL, sameModuleN=2,
              betaModel=TRUE,topology="Power",betaSize=betaSize)

# --- Estimate Purity
set.seed(1)
purity<-deNet_purity(exprM=Data$X,purity=Data$P)

}
\references{
Francesca Petralia, Li Wang , Jie Peng , Arthur Yan,  Jun Zhu and Pei Wang, A new method for constructing tumor specific gene co-expression networks based on samples with tumor purity heterogeneity (Submitted).

Jerome Friedman, Trevor Hastie, and Robert Tibshirani. Sparse inverse covariance estimation with the graphical lasso. \emph{Biostatistics}, 9(3):432-441, 2008.

Patrick Danaher, Pei Wang, and Daniela M Witten. The joint graphical lasso for inverse covariance estimation across multiple classes. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 76(2):373-397, 2014.}
