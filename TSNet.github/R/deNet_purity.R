## -- function to be used
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


calculateLL.1d<-function(uy,uz,vy,vz,c,X,A)
{
  P<-c*A
  Ux<-P*uy+(1-P)*uz;
  Vx<-P^2*vy+(1-P)^2*vz
  sum(log( sapply(1:length(X),function(i)dnorm(X[i],Ux[i],sqrt(Vx[i])))))
}


#' @export
deNet_purity<-function(exprM,purity)
{  
  #####
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
return(list(Epurity,paraM,convergence,time))

}



# 
# deNet_purity<-function(exprM,purity)
# {
#   t<-proc.time()
#   registerDoMC(10)
#   options(warn=-1)
#   convergecutoff=0.001 ### convergence criterion for EM
#   estimateP=TRUE; bootstrap=0; betaModel=TRUE
#   
#   if(bootstrap>0){
#     set.seed(bootstrap);
#     rindex<-sample(1:ncol(exprM),ncol(exprM),replace=TRUE)
#     exprM<-exprM[,rindex];
#     purity<-purity[rindex];
#     geneNames<-rownames(exprM);
#     exprM<-exprM[apply(exprM,1,sd)!=0,];
#     curgeneNames<-rownames(exprM);
#     ###standardization
#   }
#   exprM<-t(scale((exprM)))
#   geneNames<-rownames(exprM);
#   
#   curgeneNames<-rownames(exprM);
#   
#   if(estimateP)
#   {
#     for(j in 1:1000)  # --- over EM iterations 
#     {
#       
#       if(j==1){
#         cur.P<-purity;
#         cur.b<-mean(cur.P)*(1-mean(cur.P))/var(cur.P)-1
#       }else{
#         print(paste("Epurity",j,":",cor(purity,cur.P)))
#       }
#       
#       ## -- Estimate u and v (mean and variance of latent variables, i.e., uy uz vy vz)
#       
#       paraM<-foreach(pi= 1:nrow(exprM))%dopar%  # -- for each gene
# {
#   x.v=exprM[pi,] 
#   for(i in 1:10000) # -- (A Loop) over EM iterations
#   {
#     if(i==1)
#     {
#       rnew<-list(uy=mean(x.v),uz=mean(x.v),vy=var(x.v),vz=var(x.v),c=1)
#     }
#     r<-rnew;
#     temp<-Estep.1d(r$uy,r$uz,r$vy,r$vz,1,x.v,cur.P)
#     rnew<-Mstep.1d(temp$Ey,temp$Ey2,temp$Ez, temp$Ez2)
#     if(all(abs(c(rnew$uy-r$uy,rnew$uz-r$uz,rnew$vy-r$vy,rnew$vz-r$vz))<0.01))
#     {#print(paste("converge in step",i));
#       break()}
#   }
#   c(r$uy,r$uz,r$vy,r$vz)
# }
# 
# paraM<-do.call("rbind",paraM)
# colnames(paraM)<-c("uy","uz","vy","vz")
# print("EM1d")
# before.P<-cur.P;
# 
# ## -- Estimate purity (P) and delta parameter in beta distribution
# for(i in 1:10000) {                          
#   
#   if (betaModel==FALSE){ 
#     
#     new.P<-sapply(1:ncol(exprM),function(ni){
#       f<-function(x){
#         V<-paraM[,"vy"]*x^2+paraM[,"vz"]*(1-x)^2;
#         U<-paraM[,"uy"]*x+paraM[,"uz"]*(1-x);
#         0.5*sum(log(V))+0.5*sum((exprM[,ni]-U)^2/V)- log(dnorm(log(purity[ni]/(1-purity[ni])),log(x/(1-x)),sqrt(cur.b))) # purity
#       }
#       optim(cur.P[ni],f,lower=0,upper=1,method="Brent")$par})
#     new.b<-sum((log(purity/(1-purity))-log(new.P/(1-new.P)))^2)/length(purity)
#     
#   } else {
#     
#     new.P<-sapply(1:ncol(exprM),function(ni){
#       f<-function(x) {
#         V<-paraM[,"vy"]*x^2+paraM[,"vz"]*(1-x)^2;
#         U<-paraM[,"uy"]*x+paraM[,"uz"]*(1-x);
#         0.5*sum(log(V))+0.5*sum((exprM[,ni]-U)^2/V)- log(dbeta(purity[ni],x*cur.b,(1-x)*cur.b)) # purity
#       }
#       optim(cur.P[ni],f,lower=0,upper=1,method="Brent")$par})
#     f<-function(x) -sum(log(dbeta(purity,new.P*x,(1-new.P)*x)))
#     new.b<-optim(cur.b,f,lower=0)$par        
#   }
#   
#   if(all(abs(c(new.P-cur.P))<0.001))
#   {#print(paste("converge in step",i));
#     break()}
#   cur.b<-new.b;
#   cur.P<-new.P
#   print(cur.b)
#   print(mean(abs(cur.P-purity)))
# }
# print(paste("cur.b",cur.b))
# 
# 
# if(all((before.P-cur.P)<0.001))break();
# print(mean(abs(cur.P-before.P)))
#     }
# 
# if(j==10000){ stop("P not converge"); convergence="P not converge"}else{
#   print("Epurity converged"); convergence="P converged"}
# 
# Epurity<-cur.P
#   }else{Epurity<-purity};
# 
# time<-(proc.time()-t)
# return(list(Epuirty=Epurity,paraM=paraM,convergence=convergence,time=time))
# 
# }
