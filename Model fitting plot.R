library(ggplot2)

tot.df1<- A_ring.m %>% 
  mutate(., type = "Ring-stage donors") %>% 
  bind_rows(.,mutate(A_late.m, type = "Late-stage donors")) %>% 
  bind_rows(., mutate(A_para.m, type = "Total parasitaemia")) %>% 
  mutate(., host = "Acute")
tot.df2<- N_ring.m %>% 
  mutate(., type = "Ring-stage donors") %>% 
  bind_rows(.,mutate(N_late.m, type = "Late-stage donors")) %>% 
  bind_rows(., mutate(N_para.m, type = "Total parasitaemia")) %>% 
  mutate(., host = "Naive")
tot.df<-bind_rows(tot.df1,tot.df2)

A.ring<-as.data.frame(fit.fun(c(0.65,9,0.89,0,0.44),constraints,timedata,"Both")[1:8]) %>% 
  bind_cols(as.data.frame(timedata),.) %>% 
  mutate(., type = "Ring-stage donors") %>% 
  mutate(., Fit = "Khoury et al. 2017")
names(A.ring)<-c("Time","Percent","type","Fit")
A.late<-as.data.frame(fit.fun(c(0.65,9,0.89,0,0.44),constraints,timedata,"Both")[9:16]) %>% 
  bind_cols(as.data.frame(timedata),.) %>% 
  mutate(., type = "Late-stage donors") %>% 
  mutate(., Fit = "Khoury et al. 2017")
names(A.late)<-c("Time","Percent","type","Fit")
A.para<-as.data.frame(fit.fun(c(0.65,9,0.89,0,0.44),constraints,timedata,"Both")[17:24]) %>% 
  bind_cols(as.data.frame(timedata),.) %>% 
  mutate(., type = "Total parasitaemia") %>% 
  mutate(., Fit = "Khoury et al. 2017")
names(A.para)<-c("Time","Percent","type","Fit")
A.df1<-bind_rows(A.ring,A.late) %>% bind_rows(.,A.para)

N.ring<-as.data.frame(fit.fun(c(1,9,0.89,0,0.05),constraints,timedata,"Both")[1:8]) %>% 
  bind_cols(as.data.frame(timedata),.) %>% 
  mutate(., type = "Ring-stage donors") %>% 
  mutate(., Fit = "Khoury et al. 2017")
names(N.ring)<-c("Time","Percent","type","Fit")
N.late<-as.data.frame(fit.fun(c(1,9,0.89,0,0.05),constraints,timedata,"Both")[9:16]) %>% 
  bind_cols(as.data.frame(timedata),.) %>% 
  mutate(., type = "Late-stage donors") %>% 
  mutate(., Fit = "Khoury et al. 2017")
names(N.late)<-c("Time","Percent","type","Fit")
N.para<-as.data.frame(fit.fun(c(1,9,0.89,0,0.05),constraints,timedata,"Both")[17:24]) %>% 
  bind_cols(as.data.frame(timedata),.) %>% 
  mutate(., type = "Total parasitaemia") %>% 
  mutate(., Fit = "Khoury et al. 2017")
names(N.para)<-c("Time","Percent","type","Fit")
N.df1<-bind_rows(N.ring,N.late) %>% bind_rows(.,N.para)

A.ring<-as.data.frame(fit.fun(c(0.61, 10.91, 2.94, 0.69, 0.68),constraintsA,timedata,"Both")[1:8]) %>% 
  bind_cols(as.data.frame(timedata),.) %>% 
  mutate(., type = "Ring-stage donors") %>% 
  mutate(., Fit = "2")
names(A.ring)<-c("Time","Percent","type","Fit")
A.late<-as.data.frame(fit.fun(c(0.61, 10.91, 2.94, 0.69, 0.68),constraintsA,timedata,"Both")[9:16]) %>% 
  bind_cols(as.data.frame(timedata),.) %>% 
  mutate(., type = "Late-stage donors") %>% 
  mutate(., Fit = "2")
names(A.late)<-c("Time","Percent","type","Fit")
A.para<-as.data.frame(fit.fun(c(0.61, 10.91, 2.94, 0.69, 0.68),constraintsA,timedata,"Both")[17:24]) %>% 
  bind_cols(as.data.frame(timedata),.) %>% 
  mutate(., type = "Total parasitaemia") %>% 
  mutate(., Fit = "2")
  names(A.para)<-c("Time","Percent","type","Fit")
A.df2<-bind_rows(A.ring,A.late) %>% bind_rows(.,A.para)

N.ring<-as.data.frame(fit.fun(c(0.98, 11.77, 3.06, 0.64, 0.19),constraintsN,timedata,"Both")[1:8]) %>% 
  bind_cols(as.data.frame(timedata),.) %>% 
  mutate(., type = "Ring-stage donors") %>% 
  mutate(., Fit = "2")
names(N.ring)<-c("Time","Percent","type","Fit")
N.late<-as.data.frame(fit.fun(c(0.98, 11.77, 3.06, 0.64, 0.19),constraintsN,timedata,"Both")[9:16]) %>% 
  bind_cols(as.data.frame(timedata),.) %>% 
  mutate(., type = "Late-stage donors") %>% 
  mutate(., Fit = "2")
names(N.late)<-c("Time","Percent","type","Fit")
N.para<-as.data.frame(fit.fun(c(0.98, 11.77, 3.06, 0.64, 0.19),constraintsN,timedata,"Both")[17:24]) %>% 
  bind_cols(as.data.frame(timedata),.) %>% 
  mutate(., type = "Total parasitaemia") %>% 
  mutate(., Fit = "2")
names(N.para)<-c("Time","Percent","type","Fit")
N.df2<-bind_rows(N.ring,N.late) %>% bind_rows(.,N.para)

A.ring<-as.data.frame(fit.fun(c(0.98, 11.77, 3.06, 0.64, 0.19),constraintsA,timedata,"Both")[1:8]) %>% 
  bind_cols(as.data.frame(timedata),.) %>% 
  mutate(., type = "Ring-stage donors") %>% 
  mutate(., Fit = "3")
names(A.ring)<-c("Time","Percent","type","Fit")
A.late<-as.data.frame(fit.fun(c(0.55,	6.239398,	8.38,	0.83, 0.048, 0.54),constraintsA,timedata,"Both")[9:16]) %>% 
  bind_cols(as.data.frame(timedata),.) %>% 
  mutate(., type = "Late-stage donors") %>% 
  mutate(., Fit = "3")
names(A.late)<-c("Time","Percent","type","Fit")
A.para<-as.data.frame(fit.fun(c(0.55,	6.239398,	8.38,	0.83, 0.048, 0.54),constraintsA,timedata,"Both")[17:24]) %>% 
  bind_cols(as.data.frame(timedata),.) %>% 
  mutate(., type = "Total parasitaemia") %>% 
  mutate(., Fit = "3")
names(A.para)<-c("Time","Percent","type","Fit")
A.df3<-bind_rows(A.ring,A.late) %>% bind_rows(.,A.para)

N.ring<-as.data.frame(fit.fun(c(1.038197, 8.992885, 3.778515, 0.7446529, 0.05329522),constraintsN,timedata,"Both")[1:8]) %>% 
  bind_cols(as.data.frame(timedata),.) %>% 
  mutate(., type = "Ring-stage donors") %>% 
  mutate(., Fit = "3")
names(N.ring)<-c("Time","Percent","type","Fit")
N.late<-as.data.frame(fit.fun(c(1.038197, 8.992885, 3.778515, 0.7446529, 0.05329522),constraintsN,timedata,"Both")[9:16]) %>% 
  bind_cols(as.data.frame(timedata),.) %>% 
  mutate(., type = "Late-stage donors") %>% 
  mutate(., Fit = "3")
names(N.late)<-c("Time","Percent","type","Fit")
N.para<-as.data.frame(fit.fun(c(1.038197, 8.992885, 3.778515, 0.7446529, 0.05329522),constraintsN,timedata,"Both")[17:24]) %>% 
  bind_cols(as.data.frame(timedata),.) %>% 
  mutate(., type = "Total parasitaemia") %>% 
  mutate(., Fit = "3")
names(N.para)<-c("Time","Percent","type","Fit")
N.df3<-bind_rows(N.ring,N.late) %>% bind_rows(.,N.para)

fit.df1<-A.df1 %>% bind_rows(.,A.df2) %>% bind_rows(., A.df3) %>% mutate(., host = "Acute")
fit.df2<-N.df1 %>% bind_rows(., N.df2) %>% bind_rows(., N.df3) %>% mutate(., host = "Naive")
#fit.df1<-A.df1 %>% bind_rows(.,A.df2) %>%  mutate(., host = "Acute")
#fit.df2<-N.df1 %>% bind_rows(., N.df2) %>% mutate(., host = "Naive")
fit.df<-bind_rows(fit.df1, fit.df2)

acute.plot<-ggplot()+geom_point(data=tot.df,aes(x=Time, y=Percent))+
  geom_line(data=fit.df,aes(x=Time, y=Percent, colour=Fit))+
  facet_grid(host ~ type)+
  xlab("Hours")+ylab("Percent pRBC")+
  theme_bw()
