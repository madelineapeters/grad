df1<-RU.df
df2<-DU.df
change<-inf.change
unite<-"no"
infect.step2<-function(unite, df1, df2, change, PR, PD){
para.df<-as.data.frame(matrix(nrow=8, ncol=1))
if (unite == "yes"){para.df[,]<-PR+PD} else {
  para.df[1:4,]<-PR
  para.df[5:8,]<-PD
}
names(para.df)<-"P"
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

if ((unite == "yes")&(PR+PD == 0)){
  df<-bind_rows(df1,df2)
  names(df)<-"V3"
  df.final<-as.data.frame(df.final[,2])
  names(df.final)<-"V2"
  inf.table1<-as.data.frame(matrix(nrow=dim(df)[1], ncol=1))
  names(inf.table1)<-"V1"
  inf.table1[,]<-0
  df.final<-bind_cols(df.final,df,inf.table1)
} else if ((unite == "no")&(PR == 0)){
  temp.final<-filter(df.final, index>4)
  names(df1)<-"V3"
  df.final.ca<-as.data.frame(df.final[1:4,2])
  names(df.final)<-"V2"
  inf.table1<-as.data.frame(matrix(nrow=dim(df1)[1], ncol=1))
  names(inf.table1)<-"V1"
  inf.table1[,]<-0
  df.final.ca-bind_cols(df.finalca,df1,inf.table1)
  df.final<-bind_rows(df.finalca,temp.final)
} else if ((unite == "no")&(PD == 0)){
  temp.final<-filter(df.final, index<5)
  names(df2)<-"V3"
  df.finalcb<-as.data.frame(df.final[5:8,2])
  names(df.final)<-"V2"
  inf.table1<-as.data.frame(matrix(nrow=dim(df2)[1], ncol=1))
  names(inf.table1)<-"V1"
  inf.table1[,]<-0
  df.finalcb<-bind_cols(df.finalcb,df1,inf.table1)
  df.final<-bind_rows(temp.final,df.finalcb)
}


  if (unite == "yes") {
    PR<-PR+PD
    PD<-PR+PD
  }
  temp<-as.data.frame(matrix(nrow=8,ncol=1))
  for (i in 1:4){
    if ((PR*df[i,3]) > df[i,1]) {
      temp[i,1]<-0
      remain<-(PR*df[i,3])-df[i,1]
    }
  }
  
  
  
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