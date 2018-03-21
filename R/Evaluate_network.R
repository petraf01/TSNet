
Evaluate_network<-function(TSNet,Data){ 
 
  # -- Performance evaluation in estimating tumor-purity based on synthetic data
  
  # -- INPUT
  # -- TSNet: object returned by function deNet
  # -- Data:  object returned by function Synthetic_data

  # -- OUTPUT
  # -- True Positive and False Positive in detecting network edges for both tumor and non-tumor networks
 
  # --- Edges in True Networks  
  p<-dim(Data$Sigmay)[1]
  cutoff=1E-3;
  edge.y<-(abs(ginv(Data$Sigmay)[lower.tri(Data$Sigmay,diag=FALSE)])>cutoff)
  edge.z<-(abs(ginv(Data$Sigmaz)[lower.tri(Data$Sigmaz,diag=FALSE)])>cutoff)
  
  # --- Edges identified by TSNet  
  index<-rowSums(is.na(TSNet[[1]][[3]]))
  inv.nona.T<-matrix(0,p,p)  
  inv.nona.T[index!=p,index!=p]<-ginv(TSNet[[1]][[3]][index!=p,index!=p])
  edge.Ey<-(abs(inv.nona.T[lower.tri(TSNet[[1]][[3]],diag=FALSE)])>cutoff)
  
  index<-rowSums(is.na(TSNet[[1]][[4]]))
  inv.nona.T<-matrix(0,p,p)
  inv.nona.T[index!=p,index!=p]<-ginv(TSNet[[1]][[4]][index!=p,index!=p])
  edge.Ez<-(abs(inv.nona.T[lower.tri(TSNet[[1]][[4]],diag=FALSE)])>cutoff)
  
  TP.normal<-sum(edge.Ez & edge.z)
  FP.normal<-sum(edge.Ez & edge.z==FALSE)
  
  TP.tumor<-sum(edge.Ey & edge.y)
  FP.tumor<-sum(edge.Ey & edge.y==FALSE)
  
  # --- Return TP and FP rates
  return(list(TPTumor=TP.tumor,FPTumor=FP.tumor,TPNormal=TP.normal,FPNormal=FP.normal))
}
