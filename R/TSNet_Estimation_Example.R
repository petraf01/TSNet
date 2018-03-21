##### --- R Code to implement TSNet on Synthetic Data Examples

##### --------- Input Paramers for Synthetic Data Example
 p=200;     # -- p: number of genes
 rep=1;     # -- rep: number of replicates
 n=100      # -- n: number of samples for each replicate
 cluster=4;  # -- cluster: number of clusters
 parameterseed=1
 dataseed=1
 meanY=NULL;meanZ=NULL;sdY=NULL;sdZ=NULL;
 P.v=NULL;meanP=NULL;varP=NULL; 
 sameModuleN=0         # -- sameModuleN: number of overlapping modules between tumor and normal components
 topology="Power";     # -- topology: which topology between "Star" and "Power" should be utilized to simulate tumor-specific and non-tumor specific networks 
 betaModel=TRUE;       # -- beta model: TRUE if prior purity is modeled as a beta distribution; FALSE if a logistic regression is utilized
 betaSize=1            # -- Variance parameter of the beta distribution

 # --- Load Function
 source(".../R/TSNet_function.R") 

 # --- Simulate Data
 Data<-Synthetic_data(p=p,rep=rep,n=n,cluster=2,parameterseed=1,
                      dataseed=1,meanY=NULL, meanZ=NULL,sdY=NULL,sdZ=NULL, 
                      P.v=NULL,meanP=NULL,varP=NULL, sameModuleN=0,
                      betaModel=TRUE,topology=topology,betaSize=betaSize)
 
# --- Estimate tumor purity

 purity<-deNet_purity(d$X[seq(1,n),],d$P[seq(1,n)])
 print(paste("Correlation between true and estimated tumor purity = ",round(cor(purity$Epurity,d$realP),3)))

 # --- Estimate network
 rhoy=.25; rhoz=.25 # -- penalty parameters

 out<-deNet(t(d$X[seq(1,n),]),purity=purity$Epurity,paraM=purity$paraM,rhoy=rhoy,rhoz=rhoz)

 # ---- Evaluate network (TP and FP rate in estimating networks)
 Evaluate_network(out,d)
