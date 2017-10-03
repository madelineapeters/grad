library(dplyr)
df<-as.data.frame(matrix(nrow=4,ncol=3))
names(df)<-c("V1","V2","V3")
df[1:4,1]<-c(87029/289023,7938,9292433/289023,9337137)
df[1:2,2]<-0.000145
df[3:4,2]<-0.000000967
df[1:2,3]<-1
df[3:4,3]<-1
df <- df %>% mutate(probadj=V1*V2) %>% mutate(probraw=V1*V3)
df[,4]<-df[,4]/sum(df[,4])
df[,5]<-df[,5]/sum(df[,5])

df %>% group_by(V1,V2) %>% summarise(sum=sum(V3))

prob.donor.adj<-sum(df[c(1,3),4])
prob.recipient.adj<-1-prob.donor.adj

prob.donor.raw<-sum(df[c(1,3),5])
prob.recipient.raw<-1-prob.donor.raw

1000*0.00927872*0.000150
(1000-(1000*0.00927872))*0.000001

1000*0.000854202*0.000150
(1000-(1000*0.000854202))*0.000001


0.00927872*0.00150+0.990721*0.00001

0.000854202*0.00150+0.999146*0.00001
