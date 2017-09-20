library(dplyr)

gen<-120

E.pref<-0.01 #preference for erythrocytes
R.pref<-1.5 #preferences for reticulocytes

mat<-24 #average number of hours for parasite to mature in RBC

p<-5 #parasites produced from a burst donor RBC

R.0<-9150000000 #starting concentration of RBCs
R.norm<-9150000000 #normal number of RBCs in healthy host

E.R.0<-0.97 #total proportion of RBCs in circulation that are erythrocytes
R.D.0<-0.03 #total proportion that are reticulocytes

gamma<-0 #adjusts susceptibility based on age

a<-0.25 #adjustment parameter for erythropoesis in response to anaemia
age.sub<-R.norm*0.01/24
response<-"no"

kill.rate<-0
clear.rate<-0.0095

#Dataframe holding proportion of uninfected recipient RBCs of age x hours (row number)
U.df<-as.data.frame(matrix(nrow=4, ncol=1))
I.df<-as.data.frame(matrix(nrow=4, ncol=mat))

RBC.sum<-as.data.frame(matrix(nrow=gen, ncol=2))
names(RBC.sum)<-c("U","I")

P.sum<-as.data.frame(matrix(nrow=gen, ncol=1))
names(P.sum)<-c("New.inf")

infect.step<-function(df1, R.pref, E.pref){
  
    temp<-as.data.frame(matrix(nrow=4, ncol=4))
    names(temp)<-c("index", "count", "pref", "prob")
    temp[1:4,1]<-1:4
    temp[1:4,2]<-df1[1:4,1]
    temp[1:3,3]<-R.pref
    temp[4,3]<-E.pref

    RBC.total<-sum(temp[,2])
    temp[,2]<-temp[,2]/RBC.total
    pref.total<-sum(temp[,3])
    temp[,3]<-temp[,3]/pref.total
    temp[,4]<-temp[,2]*temp[,3]
    temp[,4]<-temp[,4]/sum(temp[,4])
  
  prop.inf<-as.data.frame(temp)
  
  return(prop.inf)
  
}

infect.step2<-function(df1, change, P){
    para.df1<-as.data.frame(matrix(nrow=4, ncol=1))
    para.df1[1:4,]<-P
    names(para.df1)<-"P"
    names(df1)<-"V1"
    df1<-df1 %>% 
        bind_cols(change, para.df1) %>% 
        mutate(para=P*prob) %>% 
        mutate(rem=V1-para) %>% 
        mutate(para.rem=para-V1)
    df1[df1<0]<-0
    df1<-mutate(df1,next.para=para-para.rem)
    
    remainder<-sum(df1$para.rem)
    if (remainder > 0) {
        df1.temp<-filter(df1, rem>0)
        df1.temp$P<-remainder
    
        for (i in 1:dim(df1.temp)[1]){
            df1.temp$prob[i]<-df1.temp$rem[i]*df1.temp$pref[i]
        }
        df1.temp$prob<-df1.temp$prob/sum(df1.temp$prob)
      
        for (i in 1:dim(df1.temp)[1]){
            df1.temp$para[i]<-df1.temp$P[i]*df1.temp$prob[i]
        }
        for (i in 1:dim(df1.temp)[1]){
            df1.temp$next.para[i]<-df1.temp$next.para[i]+df1.temp$para[i]
        }
    
        df1.temp$rem<-df1.temp$rem-df1.temp$P
        df1.temp$para.rem<-df1.temp$P-df1.temp$rem
        df1.2<-anti_join(df1, df1.temp, by="index")
        df1.final<-union(df1.2,df1.temp) %>% arrange(index)
        df1.final[df1.final<0]<-0
    } else {df1.final<-df1}
  return(df1.final[,c(2,8,10)])
}


retic.step<-function(res, df1, df2){
  
  retic.df<-as.data.frame(matrix(nrow=1, ncol=1))
  if (res == "yes"){
    retic.df[1,1]<-(R.norm-sum(df1))*a
  } else if (res == "no"){
    retic.df[1,1]<-((0.01*R.norm)/(24))
  }
  
}

for (t in 1:gen){
  
  if (t == 1){
    #Total concentration of recipient RBCs
    R.con<-R.0
    inf.con<-1000000
    
    #Reticulocytes
    U.df[1:3,1]<-(0.03*R.con)/3
    I.df[1:3,]<-(inf.con*0.01)/mat #assume equal distribution of parasite development stages across RBCs
    
    #Normocytes
    U.df[4,1]<-(0.97*R.con) #fills in concentration for 49+ hour age classes
    I.df[4,]<-(inf.con*0.97)/mat

  }
  
      #RBC burst of infected donor cells (for t=1 no infected cells to burst)
      Burst.tot<-sum(I.df[,mat])
      P<-Burst.tot*p

      #Infect uninfected recipient and donor RBCs
      inf.change<-as.data.frame(infect.step(U.df, R.pref, E.pref))
      inf.change<-as.data.frame(infect.step2(U.df, inf.change, P))

      #Removes newly infected RBCs from the uninfected to infected class
      U.adj<-as.data.frame(inf.change[,2])

      #Create dataframes with newly infected cells
      next.I<-as.data.frame(inf.change[,3])

      #Shifts columns so that parasites age
      I.adj<-bind_cols(next.I, I.df[,1:(mat-1)])

        if (t < 25){
            adj1<-U.adj[1,1]*(1/24)
            adj2<-U.adj[2,1]*(1/24)
            adj3<-U.adj[3,1]*(1/24)
            adj4<-U.adj[4,1]*(1/24)*(1/60)

            next1<-as.data.frame(retic.step(response, U.df, I.df))[1,1]
            next2<-adj1
            next3<-adj2
            next4<-adj3

            U.adj[1,1]<-U.adj[1,1]-adj1+next1
            U.adj[2,1]<-U.adj[2,1]-adj2+next2
            U.adj[3,1]<-U.adj[3,1]-adj3+next3
            U.adj[4,1]<-(U.adj[4,1]*(1-kill.rate))-adj4+next4
          } else if (t > 24){
            adj1<-U.age[t-24,1]*(1/24)
            adj2<-U.age[t-24,2]*(1/24)
            adj3<-U.age[t-24,3]*(1/24)
            adj4<-U.age[t-24,4]*(1/24)*(1/60)

            next1<-as.data.frame(retic.step(response, U.df, I.df))[1,1]
            next2<-adj1
            next3<-adj2
            next4<-adj3

            U.adj[1,1]<-U.adj[1,1]-adj1+next1
            U.adj[2,1]<-U.adj[2,1]-adj2+next2
            U.adj[3,1]<-U.adj[3,1]-adj3+next3
            U.adj[4,1]<-(U.adj[4,1]*(1-kill.rate))-adj4+next4
          }

        for (i in 1:mat){  
          if (t < 25){
            adj1<-I.adj[1,i]*(1/24)
            adj2<-I.adj[2,i]*(1/24)
            adj3<-I.adj[3,i]*(1/24)
            adj4<-I.adj[4,i]*(1/24)*(1/60)

            next1<-0
            next2<-adj1
            next3<-adj2
            next4<-adj3

            I.adj[1,i]<-(I.adj[1,i]*(1-clear.rate))-adj1+next1
            I.adj[2,i]<-(I.adj[2,i]*(1-clear.rate))-adj2+next2
            I.adj[3,i]<-(I.adj[3,i]*(1-clear.rate))-adj3+next3
            I.adj[4,i]<-(I.adj[4,i]*(1-clear.rate))-adj4+next4
          } else if (t > 24){
            adj1<-I.age[((t-25)*4+1),i]*(1/24)
            adj2<-I.age[((t-25)*4+2),i]*(1/24)
            adj3<-I.age[((t-25)*4+3),i]*(1/24)
            adj4<-I.age[((t-25)*4+4),i]*(1/24)*(1/60)

            next1<-0
            next2<-adj1
            next3<-adj2
            next4<-adj3

            I.adj[1,i]<-(I.adj[1,i]*(1-clear.rate))-adj1+next1
            I.adj[2,i]<-(I.adj[2,i]*(1-clear.rate))-adj2+next2
            I.adj[3,i]<-(I.adj[3,i]*(1-clear.rate))-adj3+next3
            I.adj[4,i]<-(I.adj[4,i]*(1-clear.rate))-adj4+next4
          }
        } #end i loop

      if (t == 1){
        #columns "U","I"
        RBC.sum[t,1]<-sum(U.adj)
        RBC.sum[t,2]<-sum(I.adj)
        
        I.age<-I.adj
        names(I.age)<-c(1:mat)
        U.age<-t(U.adj)
        names(U.age)<-c(1:4)
      } else {
        #Set adjusted dataframes as next-generation dataframes
        U.df<-U.adj
        temp.U<-t(U.df)
        names(temp.U)<-c(1:4)
        I.df<-I.adj
        names(I.df)<-c(1:mat)

        #columns "U","I"
        RBC.sum[t,1]<-sum(U.adj)
        RBC.sum[t,2]<-sum(I.adj)

        I.age<-bind_rows(I.age, I.df)
        names(I.age)<-c(1:mat)
        U.age<-bind_rows(U.age, temp.U)
        names(U.age)<-c(1:4)
      }
        P.sum[t,1]<-P

        print(t)


} #end loop over t
