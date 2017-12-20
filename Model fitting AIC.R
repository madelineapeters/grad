fit.fun<-function(parameters=parameters, constraints=constraints,timedata=timedata,type){
  timevec<-seq(0, 1, 0.01); #vector of times for which integration is evaluated
  paravec<-parameters #Para.corr(parameters,constraints)
  
  if (type == "Donor") {
    
    Ringfun<-function(t){
      integrand<-Vectorize(function(x){
        mat.piecewise(x,t=t,l=paravec[1],B=paravec[5],c=paravec[3],x.c=paravec[4],m=constraints[9],s=constraints[10],P0=constraints[11])
      })
      integrated<-integrate(integrand,lower=0,upper=0.5)$value
      return(integrated)
    }
    
    Latefun<-function(t){
      integrand<-Vectorize(function(x){
        mat.piecewise(x,t=t,l=paravec[1],B=paravec[5],c=paravec[3],x.c=paravec[4],m=constraints[9],s=constraints[10],P0=constraints[11])
      })
      integrated<-integrate(integrand,lower=0.5,upper=1)$value
      return(integrated)
    }
    
    #For ring curve
    Ring.out<-function(t) {sapply(t,Ringfun)}
    
    #For total parasite curve
    Late.out<-function(t) {sapply(t,Latefun)}
    ###############
    Rings<-Ring.out(seq(0,1,0.01))
    Late<-Late.out(seq(0,1,0.01))
    
    parafit<-c(Rings[c(1,17,30,42,58,75,88,100)]*100,Late[c(1,17,30,42,58,75,88,100)]*100)
    
    return(parafit)
    
  } else if (type == "Total") {
    
    Totfun<-function(t){
      integrand<-Vectorize(function(x){
        mat.piecewise(x,t=t,l=paravec[1],B=paravec[2],c=paravec[3],x.c=paravec[4],m=constraints[9],s=constraints[10],P0=constraints[12])
      })
      integrated<-integrate(integrand,lower=0,upper=1)$value
      return(integrated)
    }
    
    Tot.out<-function(t) {sapply(t,Totfun)}
    ###############
    Total<-Tot.out(seq(0,1,0.01))
    
    parafit<-c(Total[c(1,17,30,42,58,75,88,100)]*100)

    return(parafit)
    
  } else if (type == "Both") {
    
    Ringfun<-function(t){
      integrand<-Vectorize(function(x){
        mat.piecewise(x,t=t,l=paravec[1],B=paravec[5],c=paravec[3],x.c=paravec[4],m=constraints[9],s=constraints[10],P0=constraints[11])
      })
      integrated<-integrate(integrand,lower=0,upper=0.5)$value
      return(integrated)
    }
    
    Latefun<-function(t){
      integrand<-Vectorize(function(x){
        mat.piecewise(x,t=t,l=paravec[1],B=paravec[5],c=paravec[3],x.c=paravec[4],m=constraints[9],s=constraints[10],P0=constraints[11])
      })
      integrated<-integrate(integrand,lower=0.5,upper=1)$value
      return(integrated)
    }
    
    #For ring curve
    Ring.out<-function(t) {sapply(t,Ringfun)}
    
    #For total parasite curve
    Late.out<-function(t) {sapply(t,Latefun)}
    ###############
    Rings<-Ring.out(seq(0,1,0.01))
    Late<-Late.out(seq(0,1,0.01))
    
    Totfun<-function(t){
      integrand<-Vectorize(function(x){
        mat.piecewise(x,t=t,l=paravec[1],B=paravec[2],c=paravec[3],x.c=paravec[4],m=constraints[9],s=constraints[10],P0=constraints[12])
      })
      integrated<-integrate(integrand,lower=0,upper=1)$value
      return(integrated)
    }
    
    Tot.out<-function(t) {sapply(t,Totfun)}
    ###############
    Total<-Tot.out(seq(0,1,0.01))
    
    parafit<-c(Rings[c(1,17,30,42,58,75,88,100)]*100,Late[c(1,17,30,42,58,75,88,100)]*100)
    parafit2<-c(Total[c(1,17,30,42,58,75,88,100)]*100)
    parafit3<-c(parafit,parafit2)
    
    return(parafit3)
    
  } #end ifelse
} #end function

paradata<-c(N_ring.m$Percent,N_late.m$Percent)
paradata2<-N_para.m$Percent
timedata<-c(N_ring.m$Time)

fit1<-fit.fun(c(0.65,9,0.89,0,0.44),constraints,timedata,"Both")
fit2<-fit.fun(c(0.9044457,	7.750137,	1.307188,	0.5398042,	0.4772489),constraints,timedata,"Both")

#AICC (less than 40 observations) using sum squared error
n = length(c(paradata,paradata2))
k = length(c(0.65,9,0.89,0,0.44))

SSE.1 = -ode.fun(c(1,9,0.89,0,0.05),constraints,timedata,paradata,paradata2,"Both")
MSE.1 = (1/n)*SSE.1
AIC.1 = n*log(MSE.1)+2*k
AICc.1 = AIC.1 + 2*k*(k+1)/(n-k-1)

SSE.2 = -ode.fun(c(0.9816544, 11.76833, 3.060536, 0.635045, 0.1867809),constraintsA,timedata,paradata,paradata2,"Both")
MSE.2 = (1/n)*SSE.2
AIC.2 = n*log(MSE.2)+2*k
AICc.2 = AIC.2 + 2*k*(k+1)/(n-k-1)

SSE.3 = -ode.fun(c(1.038197, 8.992885, 3.778515, 0.7446529, 0.05329522),constraints,timedata,paradata,paradata2,"Both")
MSE.3 = (1/n)*SSE.3
AIC.3 = n*log(MSE.3)+2*k
AICc.3 = AIC.3 + 2*k*(k+1)/(n-k-1)

#AIC dataframe
AIC.df<-as.data.frame(matrix(nrow=3,ncol=4))
names(AIC.df)<-c("ID","MSE","k","AICc")
AIC.df[1,]<-c(1,MSE.1,5,as.double(AICc.1))
AIC.df[2,]<-c(2,MSE.2,5,as.double(AICc.2))
AIC.df[3,]<-c(3,MSE.3,5,as.double(AICc.3))
AIC.df <- AIC.df %>% mutate(delt.AICc=AICc-min(AICc)) %>% bind_cols(., as.data.frame(Weights(c(AICc.1,AICc.2, AICc.3))))
names(AIC.df)[6]<-"w"


