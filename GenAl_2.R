ode.fun<-function(parameters=parameters, constraints=constraints,timedata=timedata,paradata=paradata,paradata2=paradata2,type){
  timevec<-seq(0, 1, 0.01); #vector of times for which integration is evaluated
  paravec<-parameters #Para.corr(parameters,constraints)
  
  if (type == "Donor") {
    P0<-sum(paradata[1],paradata[9])/100
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
    
    #plot data and current fit to watch fitting in real time
    plot(timedata[1:8]/24,paradata[1:8],type="p",xlim=c(0,1),ylim=c(0,5),col="blue");
    lines(timedata[1:8]/24,paradata[9:16],type="p",col="red");
    lines(timedata[1:8]/24,parafit[1:8],type="l",col="blue");
    lines(timedata[1:8]/24,parafit[9:16],type="l",col="red");
    
    out<-(-sum((parafit-paradata)^2))
    return(out)
    
  } else if (type == "Total") {
    P1<-paradata2[1]/100
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
    
    #plot data and current fit to watch fitting in real time
    plot(timedata[1:8]/24,paradata2[1:8],type="p",xlim=c(0,1),col="green");
    lines(timedata[1:8]/24,parafit[1:8],type="l",col="green");
    
    out<-(-sum((parafit-paradata2)^2)*10)
    return(out)
    
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
    #plot data and current fit to watch fitting in real time
    plot(timedata[1:8]/24,paradata[1:8],type="p",xlim=c(0,1),ylim=c(0,5),col="blue");
    lines(timedata[1:8]/24,paradata[9:16],type="p",col="red");
    lines(timedata[1:8]/24,parafit[1:8],type="l",col="blue");
    lines(timedata[1:8]/24,parafit[9:16],type="l",col="red");
    lines(timedata[1:8]/24,paradata2[1:8],type="p",xlim=c(0,1),col="green");
    lines(timedata[1:8]/24,parafit2[1:8],type="l",col="green");
    
    out<-(-sum(c((parafit-paradata)^2),((parafit2-paradata2)^2)*10))
    return(out)
    
  } #end ifelse
} #end function

ode.fun(parameters=c(0.9,	5.39,	1,	0.6,	0.1), constraintsA,timedata,paradata,paradata2,"Both")
constraintsA<-c(0,1,0,10,0.75,0.95,0,0.042,0.42,0.28,0.0477,0.0018) #l, B, c, x.c, m, s, P0, P1 (P1 total)
constraintsN<-c(0,1,0,10,0.75,0.95,0,0.042,0.32,0.29,0.0477,0.0018) #l, B, c, x.c, m, s, P0, P1 (P1 total)

#Putting parasite percent data and times into vectors to be used later
paradata<-c(N_ring.m$Percent,N_late.m$Percent)
timedata<-c(N_ring.m$Time)

paradata2<-c(N_para.m$Percent)
timedata2<-c(A_para.m$Time)

###################################################################
#the main part, which calls the fit function 
###################################################################
theta_min <- c(l = 0.5, B = 5, c = 3, x.c = 0.6, BD = 0)
theta_max <- c(l = 1.1, B =20, c = 5, x.c = 0.9, BD = 2)

# Run the genetic algorithm
results <- ga(type = "real-valued", fitness = ode.fun,constraintsN,timedata,paradata,paradata2,"Both", 
              names = names(theta_min), 
              min = theta_min, max = theta_max,
              popSize = 100, maxiter = 50)

