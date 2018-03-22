##### --- Code to generate synthetic data and implent TSNet for tumor-purity and co-expression networks estimation

##### --------- Input Paramers for Synthetic Data Generation 
 p=200;           # -- Number of genes (i.e., number of nodes in the network).
 rep=1;           # -- Number of independent data sets (i.e., replicates).
 n=100            # -- Number of observations for each replicate.
 cluster=4;       # -- cluster: number of disjoint components in tumor and non-tumor networks
 parameterseed=1  # -- Random seed for parameter sampling.
 dataseed=1       # -- Random seed for observations sampling.
 sameModuleN=0         # -- sameModuleN: number of overlapping modules between tumor and normal components
 topology="Power";     # -- topology: which topology between "Star" and "Power" should be utilized to simulate tumor-specific and non-tumor specific networks 
 betaModel=TRUE;       # -- beta model: TRUE if prior purity is modeled as a beta distribution; FALSE if a logistic regression is utilized
 betaSize=1            # -- Variance parameter of the beta distribution
 meanY=NULL;    # -- p dimensional vector containing the mean parameter of the Gaussian distribution of the gene expression in tumor component. When NULL, meanY=0
 meanZ=NULL;    # -- p dimensional vector containing the mean parameter of the Gaussian distribution of the gene expression in non-tumor component. When NULL, meanZ=0
 sdY=NULL;      # -- p dimensional vector containing the variance parameter of the Gaussian distribution of the gene expression in tumor component. When NULL, sdY=1
 sdZ=NULL;      # -- p dimensional vector containing the variance parameter of the Gaussian distribution of the gene expression in non-tumor component. When NULL, sdZ=1      
 P.v=NULL;      # -- True level of tumor-purity. When NULL, the value of tumor purity is simulated from a Beta distribution with mean meanP and variance varP.
 meanP=NULL;    # -- Mean of the Beta distribution utilized to simulate tumor purity. When NULL, meanP=0.5.
 varP=NULL;     # -- Variance of the Beta distribution utilized to simulate tumor purity. When NULL, varP=0.04.
 
 # --- Load Functions used for Tumor purity and Network inference
 source(".../R/TSNet_function.R") 

 # --- Simulate Data using function Synthetic_data (check TSNet-manual for more information about this function)
 Data<-Synthetic_data(p=p,rep=rep,n=n,cluster=2,parameterseed=1,
                      dataseed=1,meanY=NULL, meanZ=NULL,sdY=NULL,sdZ=NULL, 
                      P.v=NULL,meanP=NULL,varP=NULL, sameModuleN=0,
                      betaModel=TRUE,topology=topology,betaSize=betaSize)
 
# --- Estimate tumor purity using function deNet_purity (check TSNet-manual for more information about this function)

 purity<-deNet_purity(exprM=d$X[seq(1,n),],purity=d$P[seq(1,n)])
 print(paste("Correlation between true and estimated tumor purity = ",round(cor(purity$Epurity,d$realP),3)))

 # --- Estimate network using function deNet (check TSNet-manual for more information about this function)
 rhoy=.25;  # -- penalty parameters for the tumor component 
 rhoz=.25   # -- penalty parameters for the non-tumor component

 out<-deNet(t(d$X[seq(1,n),]),purity=purity$Epurity,paraM=purity$paraM,rhoy=rhoy,rhoz=rhoz)

 # ---- Evaluate network (TP and FP rate in estimating networks) (check TSNet-manual for more information about this function)
 Evaluate_network(out,d)
