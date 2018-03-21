\name{deNet}
\alias{deNet}
\title{Function to estimate co-expression networks for tumor and non-tumor components from mixed expression data.}
\usage{
deNet(exprM,purity,paraM,rhoy,rhoz)
}
\arguments{
\item{exprM}{A \code{(p x n)} matrix of expression with p being the number of genes and n the number of samples.}
\item{purity}{An n-dimensional vector containing estimated purity for each sample.}
\item{paraM}{A \code{(p x 4)} matrix. The first and second columns contain the initial value of the mean of genes in tumor and non-tumor cells, respectively, and the third and fourth columns contain standard deviations for each gene in tumor and non-tumor cells, respectively.}
\item{rhoy}{Graphical lasso penalty parameter for covariance matrix in tumor cells.}
\item{rhoz}{Graphical lasso penalty parameter for covariance matrix in normal cells.}
}
\value{
\item{r}{list containing parameters estimated by the model (i.e., \code{Sigmay}, \code{Sigmaz}, \code{Uy}, \code{Uz})}
\item{LL.temp}{log-likelihood of the data for different \code{EM} iterations.}
\item{BIC}{Bayesian Information Criteria.}
\item{LL}{Log-likelihood at convergence.}
\item{EY}{Number of edges of tumor network.}
\item{EZ}{Number of edges of non-tumor network.}
\item{convergence}{value indicating whether or not the algorithm converged}
\item{purity}{estimated tumor purity.}
\item{data}{expression data matrix.}
\item{time}{total running time.}
}
\description{Estimation of co-expression networks for tumor and non-tumor components for a certain value of penalty parameters.}
\examples{

set.seed(1)
p<-20; # -- number of genes
n=10;  # -- number of samples
rhoy<-rhoz<-.8 # -- penalty parameters
paraM<-cbind(matrix(0,p,2),matrix(1,p,2))

# --- Sample Data
exprM=matrix(rnorm(p*n),p,n)
purity=runif(n)

# --- Estimate Network
deNet(exprM,purity,paraM,rhoy,rhoz)
}
\references{
Francesca Petralia, Li Wang , Jie Peng , Arthur Yan,  Jun Zhu and Pei Wang, A new method for constructing tumor specific gene co-expression networks based on samples with tumor purity heterogeneity (Submitted).

Jerome Friedman, Trevor Hastie, and Robert Tibshirani. Sparse inverse covariance estimation with the graphical lasso. \emph{Biostatistics}, 9(3):432-441, 2008.

Patrick Danaher, Pei Wang, and Daniela M Witten. The joint graphical lasso for inverse covariance estimation across multiple classes. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}, 76(2):373-397, 2014.}

