library(dplyr)

E.pref<-0.01 #preference for erythrocytes
R.pref<-1.5 #preferences for reticulocytes

B.R<-0.5 #susceptibility factor for recipient RBCs
B.D<-0.5 #susceptibility factor for donor RBCs

mat.D<-24 #average number of hours for parasite to mature in donor RBC
mat.R<-24 #average number of hours for parasite to mature in recipient RBC

P.D<-8 #parasites produced from a burst donor RBC
P.R<-8 #parasites produced from a burst recipient RBC

R.0<-8000000 #starting concentration of RBCs

T.R.0<-0.97 #total proportion of RBCs in circulation that are recipient
T.D.0<-0.03 #total proportion that are donor

gamma<-0 #adjusts susceptibility based on age

a<-0 #adjustment parameter for erythropoesis in response to anaemia

gen<-120

#Dataframe holding proportion of uninfected recipient RBCs of age x hours (row number)
RU.df<-as.data.frame(matrix(nrow=4, ncol=1))
#Dataframe holding proportion of infected donor RBCs of age x hours (row number)
DI.df<-as.data.frame(matrix(nrow=4, ncol=mat.D))
#Dataframe holding proportion of infected donor RBCs of age x hours (row number)
DU.df<-as.data.frame(matrix(nrow=4, ncol=1))
#Dataframe holding proportion of infected donor RBCs of age x hours (row number)
RI.df<-as.data.frame(matrix(nrow=4, ncol=mat.R))

RBC.sum<-as.data.frame(matrix(nrow=gen, ncol=4))
names(RBC.sum)<-c("RU","RI","DI","DU")

P.sum<-as.data.frame(matrix(nrow=gen, ncol=3))
names(P.sum)<-c("Total","R","D")

RI.age<-as.data.frame(matrix(nrow=gen, ncol=4))
DI.age<-as.data.frame(matrix(nrow=gen, ncol=4))
RU.age<-as.data.frame(matrix(nrow=gen, ncol=4))
DU.age<-as.data.frame(matrix(nrow=gen, ncol=4))

infect.step<-function(unite, df1, df2, B.R, B.D, R.pref, E.pref){
  
  if (unite == "yes"){
    temp<-as.data.frame(matrix(nrow=8, ncol=5))
    names(temp)<-c("index", "count", "pref", "B", "prob")
    temp[1:8,1]<-1:8
    temp[1:4,2]<-df1[1:4,1]
    temp[5:8,2]<-df2[1:4,1]
    temp[c(1:3,5:7),3]<-R.pref
    temp[c(4,8),3]<-E.pref
    temp[1:4,4]<-B.R
    temp[5:8,4]<-B.D
    
    RBC.total<-sum(temp[,2])
    temp[,2]<-temp[,2]/RBC.total
    pref.total<-sum(temp[,3])
    temp[,3]<-temp[,3]/pref.total
    B.total<-sum(temp[,4])
    temp[,4]<-temp[,4]/B.total
    temp[,5]<-temp[,2]*temp[,3]*temp[,4]
    temp[,5]<-temp[,5]/sum(temp[,5])
  } else if (unite == "no") {
    temp<-as.data.frame(matrix(nrow=8, ncol=5))
    names(temp)<-c("index", "count", "pref", "B", "prob")
    temp[1:8,1]<-1:8
    temp[1:4,2]<-df1[1:4,1]
    temp[5:8,2]<-df2[1:4,1]
    temp[c(1:3,5:7),3]<-R.pref
    temp[c(4,8),3]<-E.pref
    temp[1:4,4]<-1
    temp[5:8,4]<-1
    
    RBC.total1<-sum(temp[1:4,2])
    RBC.total2<-sum(temp[5:8,2])
    temp[1:4,2]<-temp[1:4,2]/RBC.total1
    temp[5:8,2]<-temp[5:8,2]/RBC.total2
    pref.total<-sum(temp[,3])
    temp[,3]<-(2*temp[,3])/pref.total
    temp[1:4,5]<-temp[1:4,2]*temp[1:4,3]*temp[1:4,4]
    temp[5:8,5]<-temp[5:8,2]*temp[5:8,3]*temp[5:8,4]
    temp[1:4,5]<-temp[1:4,5]/sum(temp[1:4,5])
    temp[5:8,5]<-temp[5:8,5]/sum(temp[5:8,5])
  }
  prop.inf<-as.data.frame(temp[,c(1,5)])
  
  return(prop.inf)
  
}

infect.step2<-function(unite, df1, df2, change, PR, PD){
  df<-bind_rows(df1,df2)
  if (unite == "yes"){
    for (i in 1:8){
      temp.df<-as.data.frame(matrix(nrow=df[i,1], ncol=2))
      temp.df[,1]<-i
      temp.df[,2]<-change[i,2]
      if (i == 1) {final.df<-temp.df} else {final.df<-bind_rows(final.df, temp.df)}
    }
    names(final.df)<-c("V1","V2")
    
    inf.temp<-sample(final.df$V1, size=(PR+PD), replace = FALSE, prob = final.df$V2)
    inf.table<-as.data.frame(table(inf.temp))
    names(inf.table)<-c("V1", "V2")
    inf.table[,1]<-as.integer(inf.table[,1])
    
    inf.df<-as.data.frame(c(1:dim(df)[1]))
    names(inf.df)<-"V1"
    inf.df[,1]<-as.integer(inf.df[,1])
    inf.ab<-full_join(inf.table, inf.df, by="V1")
    
    inf.ab[is.na(inf.ab)] <- 0
    
  } else if (unite == "no"){
    for (i in 1:4){
      temp.df<-as.data.frame(matrix(nrow=df[i,1], ncol=2))
      temp.df[,1]<-i
      temp.df[,2]<-change[i,2]
      if (i == 1) {final.df1<-temp.df} else {final.df1<-bind_rows(final.df1, temp.df)}
    }
    names(final.df1)<-c("V1","V2")
    
    for (i in 5:8){
      temp.df<-as.data.frame(matrix(nrow=df[i,1], ncol=2))
      temp.df[,1]<-i
      temp.df[,2]<-change[i,2]
      if (i == 5) {final.df2<-temp.df} else {final.df2<-bind_rows(final.df2, temp.df)}
    }
    names(final.df2)<-c("V1","V2") 
    
    inf.df1<-as.data.frame(c(1:dim(df1)[1]))
    inf.df2<-as.data.frame(dim(df1)[1]+c(1:dim(df2)[1]))
    names(inf.df1)<-"V1"
    names(inf.df2)<-"V1"
    inf.df1[,1]<-as.integer(inf.df1[,1])
    inf.df2[,1]<-as.integer(inf.df2[,1])
    
    if (PR > 0){
      inf.temp1<-sample(final.df1$V1, size=(PR), replace = FALSE, prob = final.df1$V2)
      inf.table1<-as.data.frame(table(inf.temp1))
      names(inf.table1)<-c("V1", "V2")
      inf.table1[,1]<-as.integer(inf.table1[,1])
      inf.ab1<-full_join(inf.table1, inf.df1, by="V1")
    } else {
      inf.table1<-as.data.frame(matrix(nrow=dim(df1)[1], ncol=1))
      names(inf.table1)<-"V1"
      inf.table1[,]<-0
      inf.ab1<-bind_cols(inf.df1, inf.table1)
    }
    if (PD > 0){
      inf.temp2<-sample(final.df2$V1, size=(PD), replace = FALSE, prob = final.df2$V2)
      inf.table2<-as.data.frame(table(inf.temp2))
      names(inf.table2)<-c("V1", "V2")
      inf.table2[,1]<-as.integer(inf.table2[,1])
      inf.ab2<-full_join(inf.table2, inf.df1, by="V1")
    } else {
      inf.table2<-as.data.frame(matrix(nrow=dim(df2)[1], ncol=1))
      names(inf.table2)<-"V1"
      inf.table2[,]<-0
      inf.ab2<-bind_cols(inf.df2, inf.table2)
    }
    
    inf.ab<-bind_rows(inf.ab1, inf.ab2)
    
    inf.ab[is.na(inf.ab)] <- 0
  }
  
  return(inf.ab)
}

retic.step<-function(res, df1, df2){
  
  retic.df<-as.data.frame(matrix(nrow=1, ncol=1))
  if (res == "up"){
    retic.df[1,1]<-((R.norm)/(24*60))*(1+(sum(df1, df2)/R.norm)*a)
  } else if (res == "down"){
    retic.df[1,1]<-((R.norm)/(24*60))*(1-(sum(df1, df2)/R.norm)*a)
  } else if (res == "none"){
    retic.df[1,1]<-((0.03*R.norm)/3)
  }
  
}

for (t in 1:gen){
  
  if (t == 1){
    #Total concentration of recipient RBCs
    R.con<-R.0
    D.con<-R.0
    inf.con<-1000000
    
    #Reticulocytes
    RU.df[1:3,1]<-(0.03*R.con)/3
    RI.df[1:3,1]<-(inf.con*0.03)/3
    RI.df[1:3,2:mat.R]<-0 #no infected RBCs with parasites in older age classes
    DU.df[1:3,1]<-(0.03*D.con)/3
    DI.df[1:3,1]<-(inf.con*0.03)/3
    DI.df[1:3,2:mat.D]<-0
    
    #Normocytes
    RU.df[4,1]<-(0.97*R.con) #fills in concentration for 49+ hour age classes
    RI.df[4,1]<-(inf.con*0.97)
    RI.df[4,2:mat.R]<-0
    DU.df[4,1]<-(0.97*D.con)
    DI.df[4,1]<-(inf.con*0.97)
    DI.df[4,2:mat.D]<-0
  }
  
  
  #RBC burst of infected donor cells (for t=1 no infected cells to burst)
  Burst.tot.D<-sum(DI.df[,mat.D])
  Burst.tot.R<-sum(RI.df[,mat.R])
  PD<-Burst.tot.D*P.D #should be zero here
  PR<-Burst.tot.R*P.R
  
  #Infect uninfected recipient and donor RBCs
  inf.change<-as.data.frame(infect.step("no", RU.df, DU.df, B.R, B.D, R.pref, E.pref))
  inf.change<-as.data.frame(infect.step2("no", RU.df, DU.df, inf.change, PR, PD))
  
  #Removes newly infected RBCs from the uninfected to infected class
  RU.adj<-RU.df
  if (RU.adj[1,1] < inf.change[1,1]) {RU.adj[1,1]<-0} else {RU.adj[1,1]<-RU.adj[1,1]-inf.change[1,1]}
  if (RU.adj[2,1] < inf.change[2,1]) {RU.adj[2,1]<-0} else {RU.adj[2,1]<-RU.adj[2,1]-inf.change[2,1]}
  if (RU.adj[3,1] < inf.change[3,1]) {RU.adj[3,1]<-0} else {RU.adj[3,1]<-RU.adj[3,1]-inf.change[3,1]}
  if (RU.adj[4,1] < inf.change[4,1]) {RU.adj[4,1]<-0} else {RU.adj[4,1]<-RU.adj[4,1]-inf.change[4,1]}
  
  DU.adj<-DU.df
  if (DU.adj[1,1] < inf.change[5,1]) {DU.adj[1,1]<-0} else {DU.adj[1,1]<-DU.adj[1,1]-inf.change[5,1]}
  if (DU.adj[2,1] < inf.change[6,1]) {DU.adj[2,1]<-0} else {DU.adj[2,1]<-DU.adj[2,1]-inf.change[6,1]}
  if (DU.adj[3,1] < inf.change[7,1]) {DU.adj[3,1]<-0} else {DU.adj[3,1]<-DU.adj[3,1]-inf.change[7,1]}
  if (DU.adj[4,1] < inf.change[8,1]) {DU.adj[4,1]<-0} else {DU.adj[4,1]<-DU.adj[4,1]-inf.change[8,1]}
  
  #Create dataframes with newly infected cells
  next.DI<-as.data.frame(matrix(nrow=4,ncol=1))
  next.DI[1:3,1]<-(inf.change[3,2]/3)
  next.DI[4,1]<-(inf.change[4,2])
  
  next.RI<-as.data.frame(matrix(nrow=4,ncol=1))
  next.RI[1:3,1]<-(inf.change[1,2]/3)
  next.RI[4,1]<-(inf.change[2,2])
  
  #Shifts columns so that parasites age
  DI.adj<-bind_cols(next.DI, DI.df[,1:(mat.D-1)])
  RI.adj<-bind_cols(next.RI, RI.df[,1:(mat.R-1)])
  
  #Remove proportion of oldest recipient RBCs, move to next age class
  if (RU.adj[1,1] < 5555.556) {adj1<-RU.adj[1,1]} else {adj1<-5555.556}
  if (RU.adj[2,1] < 5555.556) {adj2<-RU.adj[2,1]} else {adj2<-5555.556}
  if (RU.adj[3,1] < 5555.556) {adj3<-RU.adj[3,1]} else {adj3<-5555.556}
  if (RU.adj[4,1] < 5555.556) {adj4<-RU.adj[4,1]} else {adj4<-5555.556}
  
  next1<-as.data.frame(retic.step(response, RU.df, RI.df))[1,1]
  next2<-adj1
  next3<-adj2
  next4<-adj3
  
  RU.adj[1,1]<-RU.adj[1,1]-adj1+next1
  RU.adj[2,1]<-RU.adj[2,1]-adj2+next2
  RU.adj[3,1]<-RU.adj[3,1]-adj3+next3
  RU.adj[4,1]<-RU.adj[4,1]-adj4+next4
  
  #Remove proportion of oldest donor RBCs, move to next age class
  if (DU.adj[1,1] < 5555.556) {adj1<-DU.adj[1,1]} else {adj1<-5555.556}
  if (DU.adj[2,1] < 5555.556) {adj2<-DU.adj[2,1]} else {adj2<-5555.556}
  if (DU.adj[3,1] < 5555.556) {adj3<-DU.adj[3,1]} else {adj3<-5555.556}
  if (DU.adj[4,1] < 5555.556) {adj4<-DU.adj[4,1]} else {adj4<-5555.556}
  next1<-as.data.frame(retic.step(response, DU.df, DI.df))[1,1]
  next2<-adj1
  next3<-adj2
  next4<-adj3
  
  DU.adj[1,1]<-DU.adj[1,1]-adj1+next1
  DU.adj[2,1]<-DU.adj[2,1]-adj2+next2
  DU.adj[3,1]<-DU.adj[3,1]-adj3+next3
  DU.adj[4,1]<-DU.adj[4,1]-adj4+next4
  
  #Remove proportion of oldest recipient RBCs, move to next age class
  for (i in 1:mat.R){
    
    if (RI.adj[1,i] < 5555.556) {adj1<-RI.adj[1,i]} else {adj1<-5555.556}
    if (RI.adj[2,i] < 5555.556) {adj2<-RI.adj[2,i]} else {adj2<-5555.556}
    if (RI.adj[3,i] < 5555.556) {adj3<-RI.adj[3,i]} else {adj3<-5555.556}
    if (RI.adj[4,i] < 5555.556) {adj4<-RI.adj[4,i]} else {adj4<-5555.556}
    
    next1<-0
    next2<-adj1
    next3<-adj2
    next4<-adj3
    
    RI.adj[1,i]<-RI.adj[1,i]-adj1+next1
    RI.adj[2,i]<-RI.adj[2,i]-adj2+next2
    RI.adj[3,i]<-RI.adj[3,i]-adj3+next3
    RI.adj[4,i]<-RI.adj[4,i]-adj4+next4
    
  }
  
  #Remove proportion of oldest recipient RBCs, move to next age class
  for (i in 1:mat.D){
    
    if (DI.adj[1,i] < 5555.556) {adj1<-DI.adj[1,i]} else {adj1<-5555.556}
    if (DI.adj[2,i] < 5555.556) {adj2<-DI.adj[2,i]} else {adj2<-5555.556}
    if (DI.adj[3,i] < 5555.556) {adj3<-DI.adj[3,i]} else {adj3<-5555.556}
    if (DI.adj[4,i] < 5555.556) {adj4<-DI.adj[4,i]} else {adj4<-5555.556}
    
    next1<-0
    next2<-adj1
    next3<-adj2
    next4<-adj3
    
    DI.adj[1,i]<-DI.adj[1,i]-adj1+next1
    DI.adj[2,i]<-DI.adj[2,i]-adj2+next2
    DI.adj[3,i]<-DI.adj[3,i]-adj3+next3
    DI.adj[4,i]<-DI.adj[4,i]-adj4+next4
    
  }
  
  #Set adjusted dataframes as next-generation dataframes
  RU.df<-RU.adj
  DU.df<-DU.adj
  DI.df<-DI.adj
  RI.df<-RI.adj
  
  #columns "RU","RI","DI","DU"  
  RBC.sum[t,1]<-sum(RU.df)
  RBC.sum[t,2]<-sum(RI.df)
  RBC.sum[t,3]<-sum(DI.df)
  RBC.sum[t,4]<-sum(DU.df)
  
  RI.age[t,]<-rowSums(RI.df)
  DU.age[t,]<-rowSums(DU.df)
  DI.age[t,]<-rowSums(DI.df)
  RU.age[t,]<-rowSums(RU.df)
  
  P.sum[t,1]<-(Burst.tot.D*P.D)+(Burst.tot.R*P.R)
  P.sum[t,2]<-(Burst.tot.R*P.R)
  P.sum[t,3]<-(Burst.tot.D*P.D)
  
} #end loop over t


