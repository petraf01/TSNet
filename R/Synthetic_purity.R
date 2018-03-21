#####################
############### (a) generate network
####generate a random graph based on given degrees of nodes
RanGra<-function(degree,p){
  ##p--number of nodes;
  ##degree--degrees of each node 
  ##sum(degree) must be an even number 
  
  result<-matrix(0,p,p)
  stub<-matrix(0,sum(degree),2)
  
  count<-0
  for(i in 1:length(degree)){
    cur<-(count+1):(count+degree[i])
    stub[cur,1]<-i
    stub[cur,2]<-cur
    count<-count+degree[i]
  }
  
  index<-sample(1:sum(degree),sum(degree)/2)
  
  stub.1<-stub[index,]
  stub.2<-stub[-index,]    
  for (i in 1:nrow(stub.1)){
    cur.1<-stub.1[i,1]
    cur.2<-stub.2[i,1]
    result[cur.1,cur.2]<-1
    result[cur.2,cur.1]<-1
  }
  return(result)
  
}


#####
###get partial correlation matrix based on an adjacent network 
###
GenPar<-function(net.adj,umin,umax,flip=TRUE,factor=2){
  ##net.adj<-the adjacent network
  ##unim,umax,  the range of the orginal parcor
  ##filp=T means random sign of the parcor
  p<-nrow(net.adj)
  result<-matrix(0,p,p)
  for(i in 2:p){
    for(j in 1:(i-1)){
      cur<-runif(1,umin,umax)
      sign<-1 
      if(flip) sign<-sample(c(-1,1),1)   
      cur<-cur*sign    
      result[i,j]<-cur*net.adj[i,j]
      result[j,i]<-result[i,j]
    }
  }
  
  diag.re<-matrix(apply(abs(result),1,sum),p,p)*factor+1e-6
  result<-result/diag.re
  result<-(result+t(result))/2
  diag(result)<-1
  
  return(result)
  
}



####
##get covariance matrix from a partial correlation matrix 
##make it p.d.
GenCov<-function(ParCor.m){
  temp<-(-ParCor.m)
  diag(temp)<-1
  Sig<-solve(ParCor.m)
  
  
  p<-nrow(Sig)
  for (i in 2:p){
    for (j in 1:(i-1)){
      Sig[i,j]<-Sig[i,j]/sqrt(Sig[i,i]*Sig[j,j]) ##standardize to sigma=1
      Sig[j,i]<-Sig[i,j]
    }
  }
  diag(Sig)<-1
  
  #diagonose
  D<-eigen(Sig)$values
  if(!all(D>0)){
    print("invalid covariance matrix")
    return()
  }
  
  return(Sig)
}


########  Generate powerlaw network
PowerLaw.net<-function(p, lib.loc=NULL, umin=0.5, umax=1, power=1, flip=TRUE, factor=1.5, thre=0.1)
{
  
  test<-barabasi.game(p, power =power, directed=FALSE) 
  
  adj=matrix(0, p, p)
  edge.no=ecount(test)
  
  for(i in 1:(edge.no))
  {
    temp=as.numeric(get.edge(test, i))
    adj[temp[1], temp[2]]<-1
    adj[temp[2], temp[1]]<-1
  }
  
  ParCor.gen<-GenPar(adj,umin=umin,umax=umax,flip=flip,factor=factor)
  ##generate original "paratial correlation"
  ##check whether the generatee partial correlation is p.d?
  #all(eigen(ParCor.gen)$values>0)                         ##p.d.?
  
  ##truncation the partial correlations to make small ones larger
  temp=range(ParCor.gen[upper.tri(ParCor.gen)])
  print(temp)
  ParCor.trun<-ParCor.gen
  index<-ParCor.gen>(-thre) & ParCor.gen<(-1e-6)
  ParCor.trun[index]=runif(min=-thre,max=-thre,n=sum(index))
  index<-ParCor.gen<thre  & ParCor.gen>1e-6
  ParCor.trun[index]=runif(min=thre,max=thre,n=sum(index))
  #all(eigen(ParCor.trun)$values>0)                        ##still p.d.?
  
  ##get covariance matrix and save in Sig.Rdata
  diag(ParCor.trun)<-1
  Sig<-GenCov(ParCor.trun)
  
  ##
  return(Sig)
}


########## generate network based on an empirical distribution

Empi.net<-function(p, umin=0.5,umax=1, factor=1.5){
  
  temp<-c(1167,  634,  350,  198,  125, 69, 47, 45, 22, 11, 13, 14, 7, 5, 4, 3, 2,  1, 1, 1, 1, 1)
  emp.dist<-temp/sum(temp) ### from Rosetta network
  
  degree.gen<-sample(1:length(emp.dist), p, prob=emp.dist, replace=T)
  
  if(sum(degree.gen)/2 != round(sum(degree.gen)/2)) 
    degree.gen[p]<-degree.gen[p]+1
  
  net.adj<-RanGra(degree.gen,p)
  diag(net.adj)<-0
  
  ParCor.gen<-GenPar(net.adj,umin,umax,flip=TRUE,factor=factor) 
  ##generate original "paratial correlation"
  
  ##check whether the generatee partial correlation is p.d?
  all(eigen(ParCor.gen)$values>0)                         ##p.d.?
  
  ##truncation the partial correlations to make small ones larger
  thre<-0.1
  ParCor.trun<-ParCor.gen
  index<-ParCor.gen>(-thre) & ParCor.gen<(-1e-6)
  ParCor.trun[index]=runif(min=-thre,max=-thre,n=sum(index))
  index<-ParCor.gen<thre  & ParCor.gen>1e-6
  ParCor.trun[index]=runif(min=thre,max=thre,n=sum(index))
  all(eigen(ParCor.trun)$values>0)                        ##still p.d.?
  
  ##get covariance matrix and save in Sig.Rdata
  Sig<-GenCov(ParCor.trun) 
  
  ##
  return(Sig)
}



######## generate starr based network: three hubs 
Star.net<-function(p,hub.no=3,hub.size=16,umin=0.5,umax=1){
  
  degree.gen<-sample(1:3,p,replace=T)
  degree.gen[1:hub.no]<-hub.size          ##hub size
  
  if(sum(degree.gen)/2 != round(sum(degree.gen)/2)) 
    degree.gen[p]<-degree.gen[p]+1
  
  net.adj<-RanGra(degree.gen,p)
  diag(net.adj)<-0
  
  
  ParCor.gen<-GenPar(net.adj,umin,umax,flip=TRUE,factor=1.5) 
  ##generate original "paratial correlation"
  
  ##check whether the generatee partial correlation is p.d?
  all(eigen(ParCor.gen)$values>0)                         ##p.d.?
  
  ##truncation the partial correlations to make small ones larger
  thre<-0.1
  
  
  ParCor.trun<-ParCor.gen
  ParCor.trun[abs(ParCor.gen)<thre&abs(ParCor.gen)>0]<-thre
  all(eigen(ParCor.trun)$values>0)                        ##still p.d.?
  
  ##get covariance matrix and save in Sig.Rdata
  Sig<-GenCov(ParCor.trun) 
  
  ##
  return(Sig)
}



#######
######generate MB network: simulate connecting matrix
MB.net.connect<-function(p, max.edge=4)
{
  # p=10
  
  coordi.m<-matrix(runif(p*2, min=0, max=1), ncol=2)
  d.m<-as.matrix(dist(coordi.m))
  
  prob.m<-pnorm(d.m[upper.tri(d.m)]/p^0.5)
  connect.v<-apply(as.matrix(prob.m), 1, rbinom, n=1, size=1)
  connect.m<-matrix(0, p,p)
  connect.m[upper.tri(connect.m)]<-connect.v
  connect.m<-connect.m+t(connect.m)
  
  edge.count<-apply(connect.m>0, 1, sum)
  for(i in 1:(p*p))
  {
    cur.r<-which.max(edge.count)[1]
    cur.c<-sample((1:p)[connect.m[cur.r, ]==1], max(edge.count)-max.edge)
    connect.m[cur.r, cur.c]=0
    connect.m[cur.c, cur.r]=0
    edge.count<-apply(connect.m>0, 1, sum)
    
    if(max(edge.count)<=max.edge)
      break
  }
  
  return(connect.m)
}

###simulate partial correlation matrix
MB.GenPar<-function(net.adj,u)
{
  p<-nrow(net.adj)
  result<-matrix(0,p,p)
  result[net.adj>0]<-u
  diag(result)<-1
  
  return(result)
}



simulateDataN<-function(dn,cluster,sn,parameterseed,dataseed,meanY=NULL,
                        meanZ=NULL,sdY=NULL,sdZ=NULL,
                        P.v=NULL,meanP=NULL,varP=NULL, sameModuleN=0,
                        betaModel=TRUE,type=1,betaSize=Inf){
  ## number of clusters
  ## dn: dimension (number of genes)
  ## sn: sample size

  sizeL<-rep(floor(dn/cluster),cluster)
  sigmaFactor=1; sigmaRatio=1;
  ####simulate parameters
  if (is.null(meanY)) meanY=0
  if (is.null(meanZ)) meanZ=0
  if (is.null(sdY)) sdY=1
  if (is.null(sdZ)) sdZ=1
  set.seed(parameterseed)
  Uy<-rnorm(dn,meanY,sdY)  # -- sample mean parameter from standard normal
  Uz<-rnorm(dn,meanZ,sdZ)  # -- sample mean parameter from standard normal
  
  Sigmay<-matrix(0,dn,dn);
  Sigmaz<-Sigmay;
  
  pos<-1;
  for(i in 1:length(sizeL)) {
    while(1)
    {
      if (type==1)  x<-try(PowerLaw.net(sizeL[i], umin=0.5, umax=1, power=0.5, flip=TRUE, factor=1.5, thre=0.1))
      if (type==2) x<-try(Star.net(sizeL[i], hub.no=2,hub.size=16,umin=0.5,umax=1))
      if(length(x)>1) break()
    }       
    Sigmay[pos:(pos+sizeL[i]-1),pos:(pos+sizeL[i]-1)]<-x;
    pos<-pos+sizeL[i];
  }
  
  pos<-1;
  for(i in 1:length(sizeL))
  {
    while(1)
    {
      if (type==1) x<-try(PowerLaw.net(sizeL[i], umin=0.5, umax=1, power=0.5, flip=TRUE, factor=1.5, thre=0.1))
      if (type==2) x<-try(Star.net(sizeL[i], hub.no=2,hub.size=16,umin=0.5,umax=1))
      
      if(length(x)>1) break()
    }       
    Sigmaz[pos:(pos+sizeL[i]-1),pos:(pos+sizeL[i]-1)]<-x;
    if(i<=sameModuleN) { # --- normal = tumor
      Sigmaz[pos:(pos+sizeL[i]-1),pos:(pos+sizeL[i]-1)]<-Sigmay[pos:(pos+sizeL[i]-1),pos:(pos+sizeL[i]-1)];
      Uz[pos:(pos+sizeL[i]-1)]<-Uy[pos:(pos+sizeL[i]-1)];
    }
    pos<-pos+sizeL[i];
  }
  
  Sigmay<-Sigmay*sigmaFactor*sqrt(sigmaRatio)
  Sigmaz<-Sigmaz*sigmaFactor/sqrt(sigmaRatio)
  
  if(is.null(P.v))
  {
    if (is.null(meanP)) meanP<-0.5;
    if (is.null(varP)) varP<-0.04;
    sizeP<-meanP*(1-meanP)/varP-1
    realP<-rbeta(sn,meanP*sizeP,(1-meanP)*sizeP)
    
  } else { realP=P.v}
  
  
  ## -- Simulate observed data
  set.seed(dataseed)
  Y<-mvrnorm(sn,Uy,Sigmay)
  Z<-mvrnorm(sn,Uz,Sigmaz)
  X<-sweep(Y,1,realP,"*")+sweep(Z,1,1-realP,"*");
  
  if(betaSize==Inf) P<-realP
  if(betaSize!=Inf & betaModel==FALSE){d<-rnorm(length(realP),log(realP/(1-realP)),sqrt(betaSize)); P<-exp(d)/(1+exp(d))} # -- noisy P
  if(betaSize!=Inf & betaModel==TRUE){P<-rbeta(length(realP),betaSize*realP,betaSize*(1-realP))} # -- noisy P
  
  P[P>.99]<-.99; P[P<.01]<-.01
  print(cor(P,realP))
  
  c<-5;
  A<-P/c;
  
  # -- P is the noisy tumor purity (prior estimate)
  # -- realP the true tumor purity
  return(list(X=X,Uy=Uy,Uz=Uz,Sigmay=Sigmay,Sigmaz=Sigmaz,P=P,Y=Y,Z=Z,realP=realP));
}


Synthetic_purity<-function(p,Rep,n,cluster,sameModuleN,dataseed,topology,betaModel=TRUE,betaSize=Inf){ 
 
  # -- Performance evaluation in estimating tumor-purity based on synthetic data
  
  # -- Input
  # -- p: number of genes
  # -- Rep: number of replicates
  # -- n: number of observations for each replicate
  # -- cluster: number of disjoint components
  # -- sameModuleN: number of overlapping modules between tumor and non-tumor components
  # -- dataseed: random seed
  # -- topology: which topology between "Star" and "Power" should be utilized to simulate tumor-specific and non-tumor specific networks 

  set.seed(dataseed)
  Sim=dataseed  # -- random seed
  time<-cor.rep<-cor.true<-rep(0,Rep)
  if (topology=="Power") type=1 else type=2
  
  d<-simulateDataN(dn=p,cluster=cluster,sn=n*Rep, parameterseed=dataseed,dataseed=1,
                   sameModuleN=sameModuleN,type=type,betaSize=betaSize)
  
  betaModel=TRUE; realP.final<-TSNet.purity<-estimateP.final<-NULL;
  
  # --- Estimate purity
  for (rep in 1:Rep){ 
    realP<-d$realP[seq((rep-1)*n+1,rep*n)]
    realP.final<-c(realP.final,realP)
    estimateP<-d$P[seq((rep-1)*n+1,rep*n)]
    estimateP.final<-c(estimateP.final,estimateP)
    purity<-deNet_purity(exprM=d$X[seq((rep-1)*n+1,rep*n),],purity=estimateP);
    cor.rep[rep]<-cor(purity[[1]],realP)
    cor.true[rep]<-cor(estimateP,realP)
    time[rep]<-purity[[4]][3]
    TSNet.purity<-c(TSNet.purity,purity[[1]])
}

realP<-realP.final
priorP<-estimateP.final
return(list(TSNet.purity,priorP,realP,cor.rep,cor.true))
}
