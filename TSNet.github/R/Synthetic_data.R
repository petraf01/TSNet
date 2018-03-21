
  # -- Generate synthetic data
  # Input:
  # -- rep: number of replicates
  # -- p: number of genes
  # -- n: number of samples for each replicate
  # -- cluster: number of clusters
  # -- sameModuleN: number of overlapping modules
  # -- dataseed: seedd
  # -- topology: which topology between "Star" and "Power" should be utilized to simulate tumor-specific and non-tumor specific networks 
  # -- beta model: TRUE if prior purity is modeled as a beta distribution; FALSE if a logistic regression is utilized

# Output:
# X: Synthetic data matrix containing n*p rows corresponding to different observations and p columns corresponding to different genes
# P: prior value of tumor purity
# realP: true value of tumor purity
# Uy: Mean of the multivariate Guassian of the tumor component 
# Uz: Mean of the multivariate Guassian of the non-tumor component 
# Sigmay: Covariance matrix of the multivariate Guassian of the tumor component 
# Sigmaz: Covariance matrix of the multivariate Guassian of the non-tumor component 

##### --------- Input Paramers
 p=500;
 rep=1;
 n=200
 cluster=4;
 parameterseed=1
 dataseed=1
 meanY=NULL;meanZ=NULL;sdY=NULL;sdZ=NULL;
 P.v=NULL;meanP=NULL;varP=NULL; sameModuleN=2
 topology="Power";
 betaModel=TRUE; betaSize=1

 # --- Load Function
 source("/sc/orga/work/petraf01/TSNet.github/R/TSNet_function.R") 

 # --- Simulate Data
if (topology=="Power") type=1 else type=2

d<-simulateDataN(dn=p,cluster=cluster,sn=n*rep, parameterseed=parameterseed,dataseed=dataseed,
                 meanY=meanY,meanZ=meanZ,sdY=sdY,sdZ=sdZ,
                 P.v=P.v,meanP=meanP,varP=varP, sameModuleN=sameModuleN,betaModel=betaModel,type=type,
                 betaSize)





