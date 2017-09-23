
library(dplyr)

gen<-15

E.pref<-0.01 #preference for erythrocytes
R.pref<-1.5 #preferences for reticulocytes

BR<-0.5
BD<-0.5

mat.R<-24
mat.D<-24 #average number of hours for parasite to mature in RBC

p<-5 #parasites produced from a burst donor RBC

R.0<-9150000000 #starting concentration of RBCs
R.norm<-9150000000 #normal number of RBCs in healthy host
tran.count<-1250000000

E.R.0<-0.97 #total proportion of RBCs in circulation that are erythrocytes
R.D.0<-0.03 #total proportion that are reticulocytes

gamma<-0 #adjusts susceptibility based on age

a<-0.25 #adjustment parameter for erythropoesis in response to anaemia
age.sub<-R.norm*0.01/24
response<-"no"

kill.rate<-0
clear.rate<-0.0095

#Dataframe holding proportion of uninfected recipient RBCs of age x hours (row number)
RU.df<-as.data.frame(matrix(nrow=4, ncol=1))
RIR.df<-as.data.frame(matrix(nrow=4, ncol=mat.R))
RID.df<-as.data.frame(matrix(nrow=4, ncol=mat.D))
DU.df<-as.data.frame(matrix(nrow=4, ncol=1))
DIR.df<-as.data.frame(matrix(nrow=4, ncol=mat.R))
DID.df<-as.data.frame(matrix(nrow=4, ncol=mat.D))

RBC.sum<-as.data.frame(matrix(nrow=gen, ncol=6))
names(RBC.sum)<-c("RU", "RIR", "RID", "DU", "DIR", "DID")

P.sum<-as.data.frame(matrix(nrow=gen, ncol=3))
names(P.sum)<-c("Total", "PR", "PD")

RU.age<-as.data.frame(matrix(nrow=gen, ncol=4))
RIR.age<-as.data.frame(matrix(nrow=gen, ncol=mat.R))
RID.age<-as.data.frame(matrix(nrow=gen, ncol=mat.D))
DU.age<-as.data.frame(matrix(nrow=gen, ncol=4))
DIR.age<-as.data.frame(matrix(nrow=gen, ncol=mat.R))
DID.age<-as.data.frame(matrix(nrow=gen, ncol=mat.D))

infect.step<-function(df1, df2, R.pref, E.pref, BR, BD){
  
    temp<-as.data.frame(matrix(nrow=4, ncol=5))
    names(temp)<-c("index", "count", "pref", "B", "prob")
    temp[1:8,1]<-1:8
    temp[1:4,2]<-df1[1:4,1]
    temp[5:8,2]<-df2[1:4,1]
    temp[c(1:3,5:7),3]<-R.pref
    temp[c(4,8),3]<-E.pref
    temp[c(1:4),4]<-BR
    temp[c(5:8),4]<-BD

    RBC.total<-sum(temp[,2])
    temp[,2]<-temp[,2]/RBC.total
    pref.total<-sum(temp[,3])
    temp[,3]<-temp[,3]/pref.total
    B.total<-sum(temp[,4])
    temp[,4]<-temp[,4]/B.total
    temp[,5]<-temp[,2]*temp[,3]*temp[,4]
    temp[,5]<-temp[,5]/sum(temp[,5])
  
  prop.inf<-as.data.frame(temp)
  
  return(prop.inf)
  
}

infect.step2<-function(df1, df2, change, B.R, B.D, PR, PD){
  if (PD == PR) {PD<-PD+1}
  P.list<-c(PR, PD)
  P<-PR+PD
  for (p in P.list){
    para.df<-as.data.frame(matrix(nrow=8, ncol=1))
    para.df[1:8,]<-p
    names(para.df)<-"P"
    names(df1)<-"V1"
    names(df2)<-"V1"
    df<-bind_rows(df1,df2)
    df[,1]<-df*(p/P)
    df<-df %>% 
      bind_cols(change, para.df) %>% 
      mutate(para=P*prob) %>% 
      mutate(rem=V1-para) %>% 
      mutate(para.rem=para-V1)
    df[df<0]<-0
    df<-mutate(df,next.para=para-para.rem)
    
    remainder<-sum(df1$para.rem)
    if (remainder > 0) {
      df1.temp<-filter(df1, rem>0)
      df1.temp$P<-remainder
      
      for (i in 1:dim(df1.temp)[1]){
        df1.temp$prob[i]<-df1.temp$rem[i]*df1.temp$pref[i]*df1.temp$B
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
    } else {df1.final<-df}
    
    if (p == PR) {
      df.final<-df1.final[,c(2,9,11)]
    } else if (p == PD) {
        df.final<-bind_cols(df.final, df1.final[,c(9,11)])
        names(df.final)<-c("index", "rem.R","next.para.R","rem.D","next.para.D")
        df.final<-mutate(df.final, rem=rem.R+rem.D)
        }
  }
  return(df.final)
}

retic.step<-function(res, df1){
  
  retic.df<-as.data.frame(matrix(nrow=1, ncol=1))
  if (res == "yes"){
    retic.df[1,1]<-(R.norm-sum(df1))*a
  } else if (res == "no"){
    retic.df[1,1]<-((0.01*R.norm)/(24))
  }
  
}

for (t in 1:gen){
  
  if (t == 1){
      
    RU.df[1:4,]<-as.data.frame(t(U.age[120,1:4]))
    DU.df[1:4,]<-t(U.age[96,1:4])
    
    RIR.df[1:4,]<-I.age[477:480,1:mat.R]
    RID.df[1:4,]<-0
    DIR.df[1:4,]<-0
    DID.df[1:4,]<-I.age[381:384,1:mat.D]
      
    D.sum<-sum(DU.df, DID.df)
    DU.df<-DU.df/D.sum
    DID.df<-DID.df/D.sum
    DU.df<-DU.df*tran.count
    DID.df<-DID.df*tran.count

  }
  
    #RBC burst of infected donor cells (for t=1 no infected cells to burst)
    Burst.tot.R<-sum(RIR.df[,mat.R], DIR.df[,mat.R])
    PR<-Burst.tot.R*p
    Burst.tot.D<-sum(RID.df[,mat.D], DID.df[,mat.D])
    PD<-Burst.tot.D*p

    #Infect uninfected recipient and donor RBCs
    inf.change<-as.data.frame(infect.step(RU.df, DU.df, R.pref, E.pref, BR, BD))
    inf.change<-as.data.frame(infect.step2(RU.df, DU.df, inf.change, BR, BD, PR, PD))

    #Removes newly infected RBCs from the uninfected to infected class
    RU.adj<-as.data.frame(inf.change[1:4,6])
    DU.adj<-as.data.frame(inf.change[5:8,6])
    
    #Create dataframes with newly infected cells
    next.IR<-as.data.frame(inf.change[,3])
    names(next.IR)<-"new.R"
    next.ID<-as.data.frame(inf.change[,5])
    names(next.ID)<-"new.D"
    #Shifts columns so that parasites age
    RIR.adj<-bind_cols(as.data.frame(next.IR[1:4,]), RIR.df[,1:(mat.R-1)])
    RID.adj<-bind_cols(as.data.frame(next.ID[1:4,]), RID.df[,1:(mat.D-1)])
    DIR.adj<-bind_cols(as.data.frame(next.IR[5:8,]), DIR.df[,1:(mat.R-1)])
    DID.adj<-bind_cols(as.data.frame(next.ID[5:8,]), DID.df[,1:(mat.D-1)])

    #RBC aging (hell yeah)
    
    for (d in 1:2){
        if (d == 1){
            U.adj<-RU.adj
            U.age<-RU.age
            next1<-as.data.frame(retic.step(response, U.df))[1,1]
        } else if (d == 2) {
            U.adj<-DU.adj
            U.age<-DU.age
            next1<-0}
         if (t < 25){
            adj1<-U.adj[1,1]*(1/24)
            adj2<-U.adj[2,1]*(1/24)
            adj3<-U.adj[3,1]*(1/24)
            adj4<-U.adj[4,1]*(1/24)*(1/60)

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

            next2<-adj1
            next3<-adj2
            next4<-adj3

            U.adj[1,1]<-U.adj[1,1]-adj1+next1
            U.adj[2,1]<-U.adj[2,1]-adj2+next2
            U.adj[3,1]<-U.adj[3,1]-adj3+next3
            U.adj[4,1]<-(U.adj[4,1]*(1-kill.rate))-adj4+next4
          }
      if (d == 1) {RU.adj<-U.adj} else if (d == 2) {DU.adj<-U.adj}
    } #end loop over d
    
    for (d in 1:4){
        if (d == 1){
            U.adj<-RIR.adj
            U.age<-RIR.age
            next1<-0
            mat<-mat.R
        } else if (d == 2){
            U.adj<-RID.adj
            U.age<-RID.age
            next1<-0
            mat<-mat.D
        } else if (d == 3){
            U.adj<-DIR.adj
            U.age<-DIR.age
            next1<-0
            mat<-mat.R
        } else if (d == 4){
            U.adj<-DID.adj
            U.age<-DID.age
            next1<-0
            mat<-mat.D
        }

        for (i in 1:mat){  
          if (t < mat){
            adj1<-U.adj[1,i]*(1/24)
            adj2<-U.adj[2,i]*(1/24)
            adj3<-U.adj[3,i]*(1/24)
            adj4<-U.adj[4,i]*(1/24)*(1/60)

            next2<-adj1
            next3<-adj2
            next4<-adj3

            U.adj[1,i]<-(U.adj[1,i]*(1-clear.rate))-adj1+next1
            U.adj[2,i]<-(U.adj[2,i]*(1-clear.rate))-adj2+next2
            U.adj[3,i]<-(U.adj[3,i]*(1-clear.rate))-adj3+next3
            U.adj[4,i]<-(U.adj[4,i]*(1-clear.rate))-adj4+next4
          } else if (t > mat){
            adj1<-U.age[((t-mat+1)*4+1),i]*(1/24)
            adj2<-U.age[((t-mat+1)*4+2),i]*(1/24)
            adj3<-U.age[((t-mat+1)*4+3),i]*(1/24)
            adj4<-U.age[((t-mat+1)*4+4),i]*(1/24)*(1/60)

            next2<-adj1
            next3<-adj2
            next4<-adj3

            U.adj[1,i]<-(U.adj[1,i]*(1-clear.rate))-adj1+next1
            U.adj[2,i]<-(U.adj[2,i]*(1-clear.rate))-adj2+next2
            U.adj[3,i]<-(U.adj[3,i]*(1-clear.rate))-adj3+next3
            U.adj[4,i]<-(U.adj[4,i]*(1-clear.rate))-adj4+next4
          }
        } #end i loop
        
        if (d == 1) {
            RIR.adj<-U.adj
        } else if (d == 2) {
            RID.adj<-U.adj
        } else if (d == 3) {
            DIR.adj<-U.adj
        } else if (d == 4) {
            DID.adj<-U.adj
        } #end if/else for d
        
    } #end d loop
    
      if (t == 1){
        #columns "RU","RIR"", "RID", "DU", "DIR", "DID"
        RBC.sum[t,1]<-sum(RU.adj)
        RBC.sum[t,2]<-sum(RIR.adj)
        RBC.sum[t,3]<-sum(RID.adj)
        RBC.sum[t,4]<-sum(DU.adj)
        RBC.sum[t,5]<-sum(DIR.adj)
        RBC.sum[t,6]<-sum(DID.adj)
        
        RIR.age<-RIR.adj
        RID.age<-RID.adj
        DIR.age<-DIR.adj
        DID.age<-DID.adj
        names(RIR.age)<-c(1:mat.R)
        names(RID.age)<-c(1:mat.D)
        names(DIR.age)<-c(1:mat.R)
        names(DID.age)<-c(1:mat.D)
        RU.age<-t(RU.adj)
        names(RU.age)<-c(1:4)
        DU.age<-t(DU.adj)
        names(DU.age)<-c(1:4)
      } else {
        #Set adjusted dataframes as next-generation dataframes
        RU.df<-RU.adj
        temp.RU<-t(RU.df)
        names(temp.RU)<-c(1:4)
        DU.df<-DU.adj
        temp.DU<-t(DU.df)
        names(temp.DU)<-c(1:4)
          
        RIR.df<-RIR.adj
        names(RIR.df)<-c(1:mat.R)
        RID.df<-RID.adj
        names(RID.df)<-c(1:mat.D)
        DIR.df<-DIR.adj
        names(DIR.df)<-c(1:mat.R)
        DID.df<-DID.adj
        names(DID.df)<-c(1:mat.D)

        #columns "U","I"
        RBC.sum[t,1]<-sum(RU.adj)
        RBC.sum[t,2]<-sum(RIR.adj)
        RBC.sum[t,3]<-sum(RID.adj)
        RBC.sum[t,4]<-sum(DU.adj)
        RBC.sum[t,5]<-sum(DIR.adj)
        RBC.sum[t,6]<-sum(DID.adj)

        RIR.age<-bind_rows(RIR.age, RIR.df)
        names(RIR.age)<-c(1:mat.R)
        RID.age<-bind_rows(RID.age, RID.df)
        names(RID.age)<-c(1:mat.D)
        DIR.age<-bind_rows(DIR.age, DIR.df)
        names(DIR.age)<-c(1:mat.R)
        DID.age<-bind_rows(DID.age, DID.df)
        names(DID.age)<-c(1:mat.D)
        RU.age<-bind_rows(RU.age, temp.RU)
        names(RU.age)<-c(1:4)
        DU.age<-bind_rows(DU.age, temp.DU)
        names(DU.age)<-c(1:4)
      }
        P.sum[t,1]<-PR+PD
        P.sum[t,2]<-PR
        P.sum[t,3]<-PD

        print(t)

} #end loop over t

