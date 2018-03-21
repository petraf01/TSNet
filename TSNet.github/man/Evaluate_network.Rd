\name{Evaluate_network}
\alias{Evaluate_network}
\title{Performance evaluation in estimating tumor specific and non-tumor specific networks.}
\usage{
Evaluate_network(TSNet,Data)
}
\arguments{
\item{TSNet}{Object returned by function \code{deNet}.}
\item{Data}{Object returned by function \code{Synthetic_data}.}
}
\value{
\item{TPTumor}{True Positive Edges estimated for tumor component.}
\item{FPTumor}{False Positive Edges estimated for tumor component.}
\item{TPNormal}{True Positive Edges estimated for non-tumor component.}
\item{FPNormal}{False Positive Edges estimated for non-tumor component.}
}
\description{Calculation of true positive and false positive edges in network estimation.}
\examples{

set.seed(1)
n=200    # -- number of samples for each replicate
rep=1   # -- number of replicates
p=500    # -- network dimension
betaSize=1

# --- Sample Data
Data<-Synthetic_data(p=p,rep=rep,n=n,cluster=4,
         parameterseed=1,dataseed=1,meanY=NULL,
         meanZ=NULL,sdY=NULL,sdZ=NULL, P.v=NULL,
         meanP=NULL,varP=NULL,sameModuleN=0,betaModel=TRUE,
         topology="Power",betaSize=betaSize)

# --- Estimate Tumor purity
#purity<-deNet_purity(exprM=as.matrix(Data$X),purity=Data$P)

# --- Estimate Network
#rhoy<-rhoz<-.8 # -- penalty parameters graphical lasso
#TSNet<-deNet(exprM=t(Data$X),purity=purity$Epurity,
#           paraM=purity$paraM,rhoy=rhoy,rhoz=rhoz)

# --- Evaluate Network
#Evaluate_network(TSNet,Data) 

}
\references{
Francesca Petralia, Li Wang , Jie Peng , Arthur Yan,  Jun Zhu and Pei Wang, A new method for constructing tumor specific gene co-expression networks based on samples with tumor purity heterogeneity (Submitted).

Jerome Friedman, Trevor Hastie, and Robert Tibshirani. Sparse inverse covariance estimation with the graphical lasso. \emph{Biostatistics}, 9(3):432-441, 2008.

Patrick Danaher, Pei Wang, and Daniela M Witten. The joint graphical lasso for inverse covariance estimation across multiple classes. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 76(2):373-397, 2014.}
