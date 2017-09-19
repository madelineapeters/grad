library(dplyr)

E.pref<-0.1 #preference for erythrocytes
R.pref<-1 #preferences for reticulocytes

B.R<-0.5 #susceptibility factor for recipient RBCs
B.D<-0.5 #susceptibility factor for donor RBCs

mat.D<-24 #average number of hours for parasite to mature in donor RBC
mat.R<-24 #average number of hours for parasite to mature in recipient RBC

P.D<-5 #parasites produced from a burst donor RBC
P.R<-5 #parasites produced from a burst recipient RBC

R.0<-11700000000 #starting concentration of RBCs

T.R.0<-0.97 #total proportion of RBCs in circulation that are recipient
T.D.0<-0.03 #total proportion that are donor

gamma<-0 #adjusts susceptibility based on age

a<-0.05 #adjustment parameter for erythropoesis in response to anaemia

gen<-120

response<-"up"

R.norm<-8000000

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
  prop.inf<-as.data.frame(temp)
  
  return(prop.inf)
  
}

infect.step2<-function(unite, df1, df2, change, PR, PD){
  para.df<-as.data.frame(matrix(nrow=8, ncol=1))
  if (unite == "yes"){para.df[,]<-PR+PD} else {
    para.df[1:4,]<-PR
    para.df[5:8,]<-PD
  }
  names(para.df)<-"P"
  names(df1)<-"V1"
  names(df2)<-"V1"
  df<-bind_rows(df1,df2)
  df<-df %>% 
    bind_cols(change, para.df) %>% 
    mutate(para=P*prob) %>% 
    mutate(rem=V1-para) %>% 
    mutate(para.rem=para-V1)
  df[df<0]<-0
  df<-mutate(df,next.para=para-para.rem)
  
  if (unite == "yes"){
    remainder1<-sum(df$para.rem)
    remainder2<-remainder1
  } else {
    remainder1<-sum(df$para.rem[1:4])
    remainder2<-sum(df$para.rem[5:8])
  }
  
  df.temp<-filter(df, rem>0)
  for (i in 1:dim(df.temp)[1]){
    if (df.temp$index[i] < 5){df.temp$P[i]<-remainder1} else if (df.temp$index[i] > 4) {df.temp$P[i]<-remainder2}
  }
  if (unite == "yes"){
    for (i in 1:dim(df.temp)[1]){
      df.temp$prob[i]<-df.temp$rem[i]*df.temp$pref[i]*df.temp$B[i]
    }
    df.temp$prob<-df.temp$prob/sum(df.temp$prob)
  } else {
    for (i in 1:dim(df.temp)[1]){
      df.temp$prob[i]<-df.temp$rem[i]*df.temp$pref[i]*df.temp$B[i]
    }
    temp1<-filter(df.temp, index<5)
    temp2<-filter(df.temp, index>4)
    for (i in 1:dim(df.temp)[1]){
      if (df.temp$index[i] < 5){df.temp$prob[i]<-df.temp$prob[i]/sum(temp1$prob)} else if (df.temp$index[i] > 4) {df.temp$prob[i]<-df.temp$prob[i]/sum(temp2$prob)}
    }
  }
  
  for (i in 1:dim(df.temp)[1]){
    df.temp$para[i]<-df.temp$P[i]*df.temp$prob[i]
  }
  for (i in 1:dim(df.temp)[1]){
    df.temp$next.para[i]<-df.temp$next.para[i]+df.temp$para[i]
  }
  df.temp$rem<-df.temp$rem-df.temp$P
  df.temp$para.rem<-df.temp$P-df.temp$rem
  df.2<-anti_join(df, df.temp, by="index")
  df.final<-union(df.2,df.temp) %>% arrange(index)
  df.final[df.final<0]<-0
  
  return(df.final[,c(2,9,11)])
}

retic.step<-function(res, df1, df2){
  
  retic.df<-as.data.frame(matrix(nrow=1, ncol=1))
  if (res == "up"){
    retic.df[1,1]<-((R.norm)/(24*60))*(1+(sum(df1, df2)/R.norm)*a)
  } else if (res == "down"){
    retic.df[1,1]<-((R.norm)/(24*60))*(1-(sum(df1, df2)/R.norm)*a)
  } else if (res == "none"){
    retic.df[1,1]<-((0.03*R.norm)/(24*60))
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
  RU.adj<-as.data.frame(inf.change[1:4,2])
  DU.adj<-as.data.frame(inf.change[5:8,2])
  
  #Create dataframes with newly infected cells
  next.DI<-as.data.frame(inf.change[5:8,3])
  next.RI<-as.data.frame(inf.change[1:4,3])
  
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
  
  print(t)
  
} #end loop over t


