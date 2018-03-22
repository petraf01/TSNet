#################################### Function used to simulate data
simulateDataN<-function(dn,cluster,sn,parameterseed,dataseed,meanY=NULL,meanZ=NULL,sdY=NULL,sdZ=NULL,
                        P.v=NULL,meanP=NULL,varP=NULL, sameModuleN=0,betaModel=TRUE,type=1,betaSize=Inf){
  
  
  library(igraph)
  library(MASS)
  ## number of clusters
  ## dn: dimension (number of genes)
  ## sn: sample size
  ## c: a scale; unknown parameter to scale the percentage estimation
  sizeL<-rep(floor(dn/cluster),cluster)
  sigmaFactor=1; sigmaRatio=1
  ####simulate parameters
  set.seed(parameterseed)
  Uy<-rnorm(dn,0,1)  # -- sample mean parameter from standard normal
  Uz<-rnorm(dn,0,1)  # -- sample mean parameter from standard normal
  
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
    meanP<-0.5;
    varP<-0.04;
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
  
  # -- P is the noisy purity
  # -- realP the real purity
  return(list(X=X,A=A,Uy=Uy,Uz=Uz,Sigmay=Sigmay,Sigmaz=Sigmaz,c=c,P=P,Y=Y,Z=Z,realP=realP));
}

####################################

#################################
######09-17-07
######R functions for the JSRM project:  Part III: generate data and networks for simulation studies

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
  
  library(igraph, lib.loc=lib.loc)
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
  temp=range(ParCor.gen[upper.tri(ParCor.gen)])
  print(temp)
  ParCor.trun<-ParCor.gen
  index<-ParCor.gen>(-thre) & ParCor.gen<(-1e-6)
  ParCor.trun[index]=runif(min=-thre,max=-thre,n=sum(index))
  index<-ParCor.gen<thre  & ParCor.gen>1e-6
  ParCor.trun[index]=runif(min=thre,max=thre,n=sum(index))
  diag(ParCor.trun)<-1
  Sig<-GenCov(ParCor.trun)
  
  ##
  return(Sig)
}

######## Plot network

plot.net.adj<-function(adj.m,lib.loc=NULL, v.size=8)
{
  library(igraph, lib.loc=lib.loc)
  temp=graph.adjacency(adjmatrix=adj.m, mode="undirected")
  V(temp)$color=(degree(temp)>9)+3
  plot(temp, vertex.size=v.size, vertex.frame.color="white",layout=layout.fruchterman.reingold, vertex.label=NA, edge.color=grey(0.5))
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
  all(eigen(ParCor.gen)$values>0)                         ##p.d.?
  
  thre<-0.1
  ParCor.trun<-ParCor.gen
  index<-ParCor.gen>(-thre) & ParCor.gen<(-1e-6)
  ParCor.trun[index]=runif(min=-thre,max=-thre,n=sum(index))
  index<-ParCor.gen<thre  & ParCor.gen>1e-6
  ParCor.trun[index]=runif(min=thre,max=thre,n=sum(index))
  all(eigen(ParCor.trun)$values>0)                        ##still p.d.?
  
  Sig<-GenCov(ParCor.trun) 
  
  return(Sig)
}



######## generate star based network: three hubs 
Star.net<-function(p,hub.no=3,hub.size=16,umin=0.5,umax=1){
  
  degree.gen<-sample(1:3,p,replace=T)
  degree.gen[1:hub.no]<-hub.size          ##hub size
  
  if(sum(degree.gen)/2 != round(sum(degree.gen)/2)) 
    degree.gen[p]<-degree.gen[p]+1
  
  net.adj<-RanGra(degree.gen,p)
  diag(net.adj)<-0
  
  
  ParCor.gen<-GenPar(net.adj,umin,umax,flip=TRUE,factor=1.5) 
  all(eigen(ParCor.gen)$values>0)                         ##p.d.?
  
  thre<-0.1
  
  
  ParCor.trun<-ParCor.gen
  ParCor.trun[abs(ParCor.gen)<thre&abs(ParCor.gen)>0]<-thre
  all(eigen(ParCor.trun)$values>0)                        ##still p.d.?
  
  Sig<-GenCov(ParCor.trun) 
  
  ##
  return(Sig)
}



#######
######generate MB network: simulate connecting matrix
MB.net.connect<-function(p, max.edge=4)
{
  
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

####
###simulate partial correlation matrix
MB.GenPar<-function(net.adj,u)
{
  p<-nrow(net.adj)
  result<-matrix(0,p,p)
  result[net.adj>0]<-u
  diag(result)<-1
  
  return(result)
}



############################
################### (b): generate data
###generate data based on 1.given network (atm);2. partial correaltion (pr); 
###3. number of irrelavent nodes (m); 4. sample size (n); and 5. Gaussian assumption  
Gen.Data<-function(n,m,pr,net.adj=atm.adj){
  ###covariance matrix
  #
  SigInv1<-net.adj*pr
  diag(SigInv1)<-1
  
  #add noise
  SigInv2<-diag(rep(1,m))
  cross<-matrix(0,nrow(SigInv1),ncol(SigInv2))
  SigInv1<-cbind(SigInv1,cross)
  SigInv2<-cbind(t(cross),SigInv2)
  SigInv<-rbind(SigInv1,SigInv2)
  
  #
  Sig<-solve(SigInv)
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
  
  
  ### generate data 
  R<-t(chol(Sig))
  X<-matrix(rnorm(n*p),nrow=p)
  X<-R%*%X
  rownames(X)<-paste("g",1:p)
  return(list(X,Sig))
}




#################################### Functions called by deNet_purity and deNet

Estep.noloop<-function(Uy,Uz,Sigmay,Sigmaz,c.scale,X,A)
{
  ### Uy: vector of length p
  ### Uz: vector of length p
  ### Sigmay: covariance matrix of dimension pxp
  ### Sigmaz: covariance matrix of dimension pxp
  ### X: data matrix of dimnsion nxp
  ### A: abundance based percentage estimate; a vector of length n
  
  #Uy=r$Uy;Uz=r$Uz; Sigmay=r$Sigmay; Sigmaz=r$Sigmaz; c.scale=1; X=exprM; A=Epurity
  
  P<-c.scale*A  ### percentage estimate
  P.m=diag(P)
  p=ncol(X)
  n=nrow(X)
  
  
  #begin=proc.time()
  
  Ey<-matrix(NA,n, p)
  Ez<-Ey
  Esigmay<-matrix(0,p,p)
  Esigmaz<-Esigmay
  
  temp=eigen(Sigmay, symmetric=TRUE)
  M=temp$values
  U=t(temp$vectors)
  ### temp1=t(U)%*%diag(M)%*%U
  ### max(abs(temp1-Sigmay))
  Sy.half=t(U)%*%diag(M^0.5)
  Sy.half.inv=diag(1/M^0.5)%*%U
  ### max(abs(Sy.half%*%t(Sy.half)-Sigmay))
  
  W=Sy.half.inv%*%Sigmaz%*%t(Sy.half.inv)
  temp1=eigen(W, symmetric=TRUE)
  V=t(temp1$vectors)
  D=temp1$values
  ### temp2=t(V)%*%diag(D)%*%V
  ### max(abs(temp2-W))
  
  A=Sy.half%*%t(V)
  B=V%*%Sy.half.inv ## p by p 
  
  BX=B%*%t(X)%*%P.m ##p by n
  BUy=(B%*%Uy)%*%matrix(P^2, 1, n) ## p by n 
  BUz=(B%*%Uz)%*%matrix(P*(1-P),1,n) ## p by n
  tempB=BX-BUy-BUz ## p by n
  
  tempP2=matrix(P^2, p,n, byrow=TRUE) ##p by n
  tempC=tempP2+matrix(D, p,1)%*%matrix((1-P)^2, 1,n) ## p by n 
  tempC=1/tempC
  
  tempBC=tempB*tempC ##p by n 
  tempABC=A%*%tempBC ## p by n
  
  Ey=t(matrix(Uy, p,n)+tempABC) ##  n by p
  
  Cbar= diag(apply(tempP2*tempC, 1, mean)) ##p by p diagonal 
  Esigmay=Sigmay-A%*%Cbar%*%t(A)
  
  P.2m=diag(1/(1-P))
  Ez=P.2m%*%(X-P.m%*%Ey)
  tempP4=matrix(P^4/(1-P)^2, p,n, byrow=TRUE) ##p by n
  CbarZ= diag(apply(tempP4*tempC, 1, mean)) ##p by p 
  Esigmaz=mean(P^2/(1-P)^2)*Sigmay-A%*%CbarZ%*%t(A) ## p by p 
  
  
  Esigmay<-(Esigmay+t(Esigmay))/2  ### garantee its symmetry
  Esigmaz<-(Esigmaz+t(Esigmaz))/2  ### garantee its symmetry
  
  #end=proc.time()  
  list(Ey=Ey,Ez=Ez,Esigmay=Esigmay,Esigmaz=Esigmaz)
}



###
Estep<-function(Uy,Uz,Sigmay,Sigmaz,c.scale,X,A)
{
  ### Uy: vector of length p
  ### Uz: vector of length p
  ### Sigmay: covariance matrix of dimension pxp
  ### Sigmaz: covariance matrix of dimension pxp
  ### X: data matrix of dimnsion nxp
  ### A: abundance based percentage estimate; a vector of length n
  
  P<-c.scale*A  ### percentage estimate
  
  p=ncol(X)
  n=nrow(X)
  
  Ey<-matrix(NA,n, p)
  Ez<-Ey
  Esigmay<-matrix(0,p,p)
  Esigmaz<-Esigmay
  
  for(i in 1:n)
  {
    Ux<-P[i]*Uy+(1-P[i])*Uz;
    Sigmax<-P[i]^2*Sigmay+(1-P[i])^2*Sigmaz
    Sigmaxinv<-ginv(Sigmax)
    Sigmaxy<-P[i]*Sigmay
    Ey[i,]<-Uy+Sigmaxy%*%Sigmaxinv%*%(X[i,]-Ux)
    Esigmayi<-Sigmay-Sigmaxy%*%Sigmaxinv%*%t(Sigmaxy)
    Esigmay<-Esigmay+Esigmayi
    Ez[i,]<-(X[i,]-P[i]*Ey[i,])/(1-P[i]);
    Esigmaz<-Esigmaz+P[i]^2/(1-P[i])^2*Esigmayi;
  }
  Esigmay<-Esigmay/nrow(X);
  Esigmaz<-Esigmaz/nrow(X);
  
  Esigmay<-(Esigmay+t(Esigmay))/2  ### garantee its symmetry
  Esigmaz<-(Esigmaz+t(Esigmaz))/2  ### garantee its symmetry
  
  list(Ey=Ey,Ez=Ez,Esigmay=Esigmay,Esigmaz=Esigmaz)
}


#############################################################################################################################
Mstep.v2<-function(Ey,Ez,Esigmay,Esigmaz,rhoy=rhoy, rhoz=rhoz, Uy,Uz)
{
  ### Ey: matrix of dimension nxp
  ### Ez: matrix of dimension nxp
  ### Esigmay: Conditional Expectation of the covariance matrix of dimension pxp
  ### Esigmaz: Conditional Expectation of the covariance matrix of dimension pxp
  
  p=ncol(Ey);
  n=nrow(Ey);
  
  if(missing(Uy))Uy=apply(Ey,2,mean);
  if(missing(Uz))Uz=apply(Ez,2,mean);
  
  temp=matrix(Uy, nrow=n, ncol=p, byrow=T)
  Sy=t(Ey-temp)%*%(Ey-temp)/n+Esigmay
  Sigmay<-glasso(Sy,rhoy, penalize.diagonal=FALSE)$w;
  
  temp=matrix(Uz, nrow=n, ncol=p, byrow=T)
  Sz=t(Ez-temp)%*%(Ez-temp)/n+Esigmaz
  Sigmaz<-glasso(Sz,rhoz, penalize.diagonal=FALSE)$w;
  
  list(Uy=Uy,Uz=Uz,Sigmay=Sigmay,Sigmaz=Sigmaz)
}


######
Mstep<-function(Ey,Ez,Esigmay,Esigmaz,rho=rho,Uy,Uz)
{
  ### Ey: matrix of dimension nxp
  ### Ez: matrix of dimension nxp
  ### Esigmay: Conditional Expectation of the covariance matrix of dimension pxp
  ### Esigmaz: Conditional Expectation of the covariance matrix of dimension pxp
  
  p=ncol(Ey);
  n=nrow(Ey);
  
  if(missing(Uy))Uy=apply(Ey,2,mean);
  if(missing(Uz))Uz=apply(Ez,2,mean);
  
  temp<-matrix(0,p, p)
  for(i in 1:n)
  {
    temp<-temp+(Ey[i,]-Uy)%*%t(Ey[i,]-Uy)
  }
  Sy<-temp/n+Esigmay;
  Sigmay<-glasso(Sy,rho)$w;
  
  temp<-matrix(0,p, p)
  for(i in 1:n)
  {
    temp<-temp+(Ez[i,]-Uz)%*%t(Ez[i,]-Uz)
  }
  Sz<-temp/n+Esigmaz;
  Sigmaz<-glasso(Sz,rho)$w;
  
  list(Uy=Uy,Uz=Uz,Sigmay=Sigmay,Sigmaz=Sigmaz)
}


calculateLL<-function(Uy,Uz,Sigmay,Sigmaz,c.scale,X,A)
{
  P<-c.scale*A
  sum(sapply(1:nrow(X),function(i)
  {
    Ux<-P[i]*Uy+(1-P[i])*Uz;
    Sigmax<-P[i]^2*Sigmay+(1-P[i])^2*Sigmaz
    dmvnorm(X[i,], Ux, Sigmax, log = TRUE)
  }))
}

#######################
Estep.1d<-function(uy,uz,vy,vz, c,X,A)
{
  P<-c*A
  Ux<-P*uy+(1-P)*uz;
  Vx<-P^2*vy+(1-P)^2*vz
  Sxy<-P*vy
  Ey<-uy+(X-Ux)*Sxy/Vx
  Vy<-vy-Sxy^2/Vx
  
  Ey2<-Vy+Ey^2
  
  Sxz<-(1-P)*vz
  Ez<-uz+(X-Ux)*Sxz/Vx
  Vz<-vz-Sxz^2/Vx
  Ez2<-Vz+Ez^2
  #Ez<-(X-Ey*P)/(1-P);
  #Ez2<-(X^2-2*P*X*Ey+P^2*Ey2)/(1-P)^2;
  
  list(Ey=Ey,Ey2=Ey2, Ez=Ez, Ez2=Ez2)
}

Mstep.1d<-function(Ey,Ey2,Ez,Ez2)
{
  uy<-mean(Ey)
  uz<-mean(Ez)
  vy<-mean(Ey2)-uy^2
  vz<-mean(Ez2)-uz^2
  list(uy=uy,uz=uz,vy=vy,vz=vz)
}

####
calculateLL.1d<-function(uy,uz,vy,vz,c,X,A)
{
  P<-c*A
  Ux<-P*uy+(1-P)*uz;
  Vx<-P^2*vy+(1-P)^2*vz
  sum(log( sapply(1:length(X),function(i)dnorm(X[i],Ux[i],sqrt(Vx[i])))))
}

#####
calMSE<-function(est, tru)
{
  result=rep(0, 4)
  result[1]=mean((est$Uy-tru$Uy)^2)
  result[2]=mean((est$Uz-tru$Uz)^2)
  result[3]=mean((est$Sigmay-tru$Sigmay)^2)
  result[4]=mean((est$Sigmaz-tru$Sigmaz)^2)
  return(result)
}


deNet_purity<-function(exprM,purity)
{
  
  library("MASS")
  library("mvtnorm")
  library("clusterGeneration")
  library("glasso")
  library("doMC")
  library("foreach")
  
  estimateP=TRUE; bootstrap=0; betaModel=TRUE
  t<-proc.time()
  registerDoMC(10)
  options(warn=-1)
  convergecutoff=0.001 ### convergence criterion for EM
  if(bootstrap>0){
    set.seed(bootstrap);
    rindex<-sample(1:ncol(exprM),ncol(exprM),replace=TRUE)
    exprM<-exprM[,rindex];
    purity<-purity[rindex];
    geneNames<-rownames(exprM);
    exprM<-exprM[apply(exprM,1,sd)!=0,];
    curgeneNames<-rownames(exprM);
    ###standardization
  }
  exprM<-t(scale((exprM)))
  geneNames<-rownames(exprM);
  curgeneNames<-rownames(exprM);
  
  if(estimateP)
  {
    for(j in 1:1000)  # --- over EM iterations 
    {
      
      if(j==1){
        cur.P<-purity;
        cur.b<-mean(cur.P)*(1-mean(cur.P))/var(cur.P)-1
      }else{
        print(paste("Epurity",j,":",cor(purity,cur.P)))
      }
      
      ## -- Estimate u and v (mean and variance of latent variables, i.e., uy uz vy vz)
      
      paraM<-foreach(pi= 1:nrow(exprM))%dopar%  # -- for each gene
{
  x.v=exprM[pi,] 
  for(i in 1:10000) # -- (A Loop) over EM iterations
  {
    if(i==1)
    {
      rnew<-list(uy=mean(x.v),uz=mean(x.v),vy=var(x.v),vz=var(x.v),c=1)
    }
    r<-rnew;
    temp<-Estep.1d(r$uy,r$uz,r$vy,r$vz,1,x.v,cur.P)
    rnew<-Mstep.1d(temp$Ey,temp$Ey2,temp$Ez, temp$Ez2)
    if(all(abs(c(rnew$uy-r$uy,rnew$uz-r$uz,rnew$vy-r$vy,rnew$vz-r$vz))<0.01))
    {#print(paste("converge in step",i));
      break()}
  }
  c(r$uy,r$uz,r$vy,r$vz)
}

paraM<-do.call("rbind",paraM)
colnames(paraM)<-c("uy","uz","vy","vz")
print("EM1d")
before.P<-cur.P;
## -- Estimate purity (P) and delta parameter in beta distribution
for(i in 1:10000) {                          
  
  if (betaModel==FALSE){ 
    
    new.P<-sapply(1:ncol(exprM),function(ni){
      f<-function(x){
        V<-paraM[,"vy"]*x^2+paraM[,"vz"]*(1-x)^2;
        U<-paraM[,"uy"]*x+paraM[,"uz"]*(1-x);
        0.5*sum(log(V))+0.5*sum((exprM[,ni]-U)^2/V)- log(dnorm(log(purity[ni]/(1-purity[ni])),log(x/(1-x)),sqrt(cur.b))) # purity
      }
      optim(cur.P[ni],f,lower=0,upper=1,method="Brent")$par})
    new.b<-sum((log(purity/(1-purity))-log(new.P/(1-new.P)))^2)/length(purity)
    
  } else {
    
    new.P<-sapply(1:ncol(exprM),function(ni){
      f<-function(x) {
        V<-paraM[,"vy"]*x^2+paraM[,"vz"]*(1-x)^2;
        U<-paraM[,"uy"]*x+paraM[,"uz"]*(1-x);
        0.5*sum(log(V))+0.5*sum((exprM[,ni]-U)^2/V)- log(dbeta(purity[ni],x*cur.b,(1-x)*cur.b)) # purity
      }
      optim(cur.P[ni],f,lower=0,upper=1,method="Brent")$par})
    f<-function(x) -sum(log(dbeta(purity,new.P*x,(1-new.P)*x)))
    new.b<-optim(cur.b,f,lower=0)$par        
  }
  
  if(all(abs(c(new.P-cur.P))<0.001))
  {#print(paste("converge in step",i));
    break()}
  cur.b<-new.b;
  cur.P<-new.P
  print(cur.b)
  print(mean(abs(cur.P-purity)))
}
print(paste("cur.b",cur.b))


if(all((before.P-cur.P)<0.001))break();
print(mean(abs(cur.P-before.P)))
    }

if(j==10000){ stop("P not converge"); convergence="P not converge"}else{
  print("Epurity converged"); convergence="P converged"}

Epurity<-cur.P
  }else{Epurity<-purity};

time<-(proc.time()-t)
return(list(Epurity=Epurity,paraM=paraM,convergence=convergence,time=time))

}



deNet<-function(exprM,purity,paraM,rhoy,rhoz)
{
  
  # -- exprM: (p x n) matrix of expression with p being the number of genes and n the number of samples
  # --- purity: n-dimensional vector containing estimated purity for each sample
  # --- paraM: (p x 4) matrix. The first and second columns contain the initial value of the mean of genes in tumor and normal cells respectively, the third and fourth columns contain standard deviations for each gene in tumor and normal cells respectively
  # --- rhoy: graphical lasso penalty parameter for covariance matrix in tumor cells
  # --- rhoz: graphical lasso penalty parameter for covariance matrix in normal cells
  
  library("MASS")
  library("mvtnorm")
  library("clusterGeneration")
  library("glasso")
  library("doMC")
  library("foreach")
  
  estimateP=NULL; bootstrap=NULL
  label="test";
  CV=10; CVIndex=0; insilico="NO"; BIC=TRUE;
  if (is.null(estimateP)) estimateP=FALSE
  if (is.null(bootstrap)) bootstrap=FALSE
  
  
  t<-proc.time()
  
  registerDoMC(10)
  options(warn=-1)
  convergecutoff=0.1 ### convergence criterion for EM
  if(bootstrap>0){
    set.seed(bootstrap);
    rindex<-sample(1:ncol(exprM),ncol(exprM),replace=TRUE)
    exprM<-exprM[,rindex];
    purity<-purity[rindex];
    geneNames<-rownames(exprM);
    exprM<-exprM[apply(exprM,1,sd)!=0,];
    curgeneNames<-rownames(exprM);
    ###standardization
    exprM<-t(scale(t(exprM)))
  }
  geneNames<-rownames(exprM);
  curgeneNames<-rownames(exprM);
  
  if(CVIndex>0)
  {
    set.seed(0);
    rIndex<-sample(1:ncol(exprM));
    size<-floor(ncol(exprM)/CV);
    left<-ncol(exprM)%%CV;
    if(left>0)index<-c(rep(1:CV,size),1:left);
    if(left==0)index<-rep(1:CV,size)
    testIndex<-rIndex[which(index==CVIndex)];
    trainIndex<-rIndex[which(index!=CVIndex)];
    
    exprM<-exprM[apply(exprM[,trainIndex],1,sd)!=0,];
    curgeneNames<-rownames(exprM);
    ###standardization
    exprM<-t(scale(t(exprM)))
  }
  
  if (CVIndex==0) exprM<-t(scale(t(exprM)))
  
  if(estimateP){
    for(j in 1:1000)  # --- over EM iterations
    {
      
      if(j==1){
        cur.P<-purity;
        cur.b<-mean(cur.P)*(1-mean(cur.P))/var(cur.P)-1
      }else{
        
        #   print(paste("Epurity",j,":",cor(purity,cur.P)))
      }
      
      ## -- Estimate u and v (mean and variance of latent variables, i.e., uy uz vy vz)
      
      paraM<-foreach(pi= 1:nrow(exprM))%dopar%  # -- for each gene
{
  x.v=exprM[pi,]
  for(i in 1:10000) # -- (A Loop) over EM iterations
  {
    if(i==1)
    {
      rnew<-list(uy=mean(x.v),uz=mean(x.v),vy=var(x.v),vz=var(x.v),c=1)
    }
    r<-rnew;
    temp<-Estep.1d(r$uy,r$uz,r$vy,r$vz,1,x.v,cur.P)
    rnew<-Mstep.1d(temp$Ey,temp$Ey2,temp$Ez, temp$Ez2)
    if(all(abs(c(rnew$uy-r$uy,rnew$uz-r$uz,rnew$vy-r$vy,rnew$vz-r$vz))<0.01))
    {#print(paste("converge in step",i));
      break()}
  }
  c(r$uy,r$uz,r$vy,r$vz)
}

paraM<-do.call("rbind",paraM)
colnames(paraM)<-c("uy","uz","vy","vz")
#      print("EM1d")
before.P<-cur.P;

## -- Estimate purity (P) and delta parameter in beta distribution
for(i in 1:10000)  # -- (B Loop) CHECK this loop should be the same as (A Loop)
{
  new.P<-sapply(1:ncol(exprM),function(ni){
    f<-function(x)
    {
      V<-paraM[,"vy"]*x^2+paraM[,"vz"]*(1-x)^2;
      U<-paraM[,"uy"]*x+paraM[,"uz"]*(1-x);
      
      0.5*sum(log(V))+0.5*sum((exprM[,ni]-U)^2/V)-
        log(dbeta(purity[ni],x*cur.b,(1-x)*cur.b)) # purity
    }
    optim(cur.P[ni],f,lower=0.1,upper=0.9,method="Brent")$par})
  
  
  f<-function(x)sum(-log(dbeta(purity,new.P*x,(1-new.P)*x))) # CHECK: purity and new.P should be inverted??
  new.b<-optim(cur.b,f,lower=0)$par
  if(all(abs(c(new.P-cur.P))<0.001))
  {#print(paste("converge in step",i));
    break()}
  cur.b<-new.b;
  cur.P<-new.P
  #   print(cur.b)
  #    print(mean(abs(cur.P-purity)))
}
# print(paste("cur.b",cur.b))


if(all((before.P-cur.P)<0.001))break();
#  print(mean(abs(cur.P-before.P)))
    }

if(j==10000){ stop("P not converge")}else{print("Epurity converged")}

Epurity<-cur.P


  }else{
    # ---- ************************************************************************ ----#
    ## -- Estimate u and v using univariate model
    paraM<-foreach(pi= 1:nrow(exprM))%dopar%  # -- for each gene
{
  x.v=exprM[pi,]
  for(i in 1:10000) # -- (A Loop) over EM iterations
  {
    if(i==1)
    {
      rnew<-list(uy=mean(x.v),uz=mean(x.v),vy=var(x.v),vz=var(x.v),c=1)
    }
    r<-rnew;
    temp<-Estep.1d(r$uy,r$uz,r$vy,r$vz,1,x.v,purity)
    rnew<-Mstep.1d(temp$Ey,temp$Ey2,temp$Ez, temp$Ez2)
    if(all(abs(c(rnew$uy-r$uy,rnew$uz-r$uz,rnew$vy-r$vy,rnew$vz-r$vz))<0.01))
    {break()}
  }
  c(r$uy,r$uz,r$vy,r$vz)
}

paraM<-do.call("rbind",paraM)
colnames(paraM)<-c("uy","uz","vy","vz")
#   print("EM1d")
Epurity<-purity

# ---- ************************************************************************ ----#
  };

exprM<-t(exprM)
## ---  Analyze through Multiv-EM
## --- Estimate mean and covariance matrix for y and z
LL.temp=c()
i<-1
while(i<5000)
{
  #   print(paste("Multiv-EM",i))
  if(i==1)
  {
    # -- initialize values
    rnew<-list(Uy=paraM[,"uy"],Uz=paraM[,"uz"],Sigmay=diag(paraM[,"vy"]),Sigmaz=diag(paraM[,"vz"]))
    LL.temp<-calculateLL(rnew$Uy,rnew$Uz,rnew$Sigmay,rnew$Sigmaz,c.scale=1,exprM,Epurity)
  }
  r<-rnew;
  cur.para<-Estep.noloop(r$Uy,r$Uz,r$Sigmay,r$Sigmaz,c.scale=1,exprM,Epurity)
  
  rnew<-Mstep.v2(cur.para$Ey,cur.para$Ez,cur.para$Esigmay,cur.para$Esigmaz,rhoy=rhoy, rhoz=rhoz)
  
  # --- remove bad nodes
  v<-cbind(diag(rnew$Sigmay),diag(rnew$Sigmaz))
  badpoints<-which(v[,1]>20|v[,2]>20|is.na(v[,1])|is.na(v[,2]))
  
  if(length(badpoints)>0)
  {
    
     exprM<-exprM[,-c(badpoints)];
    curgeneNames<-curgeneNames[-c(badpoints)];
    rnew$Uy<-rnew$Uy[-c(badpoints)];
    rnew$Uz<-rnew$Uz[-c(badpoints)];
    rnew$Sigmay<-rnew$Sigmay[-c(badpoints),-c(badpoints)];
    rnew$Sigmaz<-rnew$Sigmaz[-c(badpoints),-c(badpoints)];
    
    r$Uy<-r$Uy[-c(badpoints)];
    r$Uz<-r$Uz[-c(badpoints)];
    r$Sigmay<-r$Sigmay[-c(badpoints),-c(badpoints)];
    r$Sigmaz<-r$Sigmaz[-c(badpoints),-c(badpoints)];
    #        next;
  }
  
  LL.temp<-c(LL.temp, calculateLL(rnew$Uy,rnew$Uz,rnew$Sigmay,rnew$Sigmaz,c.scale=1,exprM,Epurity))
  
  if(abs(LL.temp[i+1]-LL.temp[i])<convergecutoff) # -- break when log likelihood converges
  { break()}else{
    
  }
  i<-i+1;
  
}
if(i==5000) {converge<-FALSE}else{converge<-TRUE}

if(length(curgeneNames)!=length(geneNames)){
  t<-rep(NA,length(geneNames));
  names(t)<-geneNames;
  t[curgeneNames]<-r$Uy;
  r$Uy<-t;
  t<-rep(NA,length(geneNames));
  names(t)<-geneNames;
  t[curgeneNames]<-r$Uz;
  r$Uz<-t;
  t<-matrix(NA,length(geneNames),length(geneNames),dimnames=list(geneNames,geneNames));
  t[curgeneNames,curgeneNames]<-r$Sigmay;
  r$Sigmay<-t;
  t<-matrix(NA,length(geneNames),length(geneNames),dimnames=list(geneNames,geneNames));
  t[curgeneNames,curgeneNames]<-r$Sigmaz;
  r$Sigmaz<-t;
}


time=(proc.time()-t)

edge.y<-sum(abs(ginv(rnew$Sigmay)[upper.tri(rnew$Sigmay)])>0.001)
edge.z<-sum(abs(ginv(rnew$Sigmaz)[upper.tri(rnew$Sigmaz)])>0.001)
BIC<-(-2*LL.temp[i]+(edge.y+edge.z)*log(dim(exprM)[1]))

return(list(r=r,LL.temp,BIC=BIC,LL=LL.temp[i],EY=edge.y,EZ=edge.z,
            converge=converge,purity=Epurity,data=exprM,time=time))
}


Synthetic_data=function(p,rep,n,cluster,parameterseed,dataseed,meanY=NULL,meanZ=NULL,sdY=NULL,sdZ=NULL,
                        P.v=NULL,meanP=NULL,varP=NULL, sameModuleN=0,topology="Power",betaModel=TRUE,
                        betaSize=Inf){
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
  
  # --- Simulate Data
  if (topology=="Power") type=1 else type=2
  
  d<-simulateDataN(dn=p,cluster=cluster,sn=n*rep, parameterseed=parameterseed,dataseed=dataseed,
                   meanY=meanY,meanZ=meanZ,sdY=sdY,sdZ=sdZ,
                   P.v=P.v,meanP=meanP,varP=varP, sameModuleN=sameModuleN,betaModel=betaModel,type=type,
                   betaSize)
  return(d)
  
  
}

Evaluate_network<-function(TSNet,Data){ 
  
  # -- Performance evaluation in estimating tumor-purity based on synthetic data
  # -- TSNet: object returned by function deNet
  # -- Data:  object returned by function Synthetic_data
  
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
  
  TP.normal<-sum(edge.Ez & edge.z)/sum(edge.z)
  FP.normal<-sum(edge.Ez & edge.z==FALSE)/sum(edge.z==FALSE)
  
  TP.tumor<-sum(edge.Ey & edge.y)/sum(edge.y)
  FP.tumor<-sum(edge.Ey & edge.y==FALSE)/sum(edge.y==FALSE)
  
  # --- Return TP and FP rates
  return(list(TPTumor=TP.tumor,FPTumor=FP.tumor,TPNormal=TP.normal,FPNormal=FP.normal))
}
