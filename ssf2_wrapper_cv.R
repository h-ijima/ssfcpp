setwd("/Users/ijima/Desktop/Future_projection/projection_cons_F1214")

library(r4ss)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Rcpp)
library(data.table)
library(gridExtra)
library(reshape2)
source("make_input_data.R") #Read r functions to make input data for the future projection
Sys.setenv("PKG_CXXFLAGS"="-std=c++0x") 
sourceCpp("ssf2_V2.cpp") #Read future projection program
load("res.Rdata") #Read the result of SS3

#Make input data sets******************************************************************************************************************************************
d1 = get_naa(myreplist)    #Get n at age
d2 = get_caa(myreplist)    #Get catch at age 
d3 = get_faa(d1$naa,d2,myreplist)    #Get f at age
#Average f 2010-2012
tmp = d3  %>% filter(Yr>=2012,Yr<=2014) %>% group_by(Gender,age) %>% summarise(f=mean(faa))
f1214 = tmp$f
rm(tmp)
d4 = get_maa(myreplist)    #Get maturty at age
d5 = get_rec(myreplist)    #Get spawner recruitment relationship
d6 = get_ac(myreplist)     #Get autocorrelation data
#Average catch weight 2013-2015
#d7 = get_cat_b(myreplist)
#tmp = d7 %>% filter(Year>=2012,Year<=2014)  %>% summarise(catb=mean(Retain_Catch_B))
#catb1315 = tmp$catb
#rm(tmp)
catb1014 = 82432 #This is actual catch.
d8 = get_cv(myreplist)


int_n = list()
tmp = rep(0,(d4$age+1)*2)
maxage = (d4$age+1)*d4$gender
for (i in 1:500) {
  for (j in 1:maxage) {
    tmp[j] = rnorm(1, d1$int_n[j], d1$int_n[j]*d8)
  }
  int_n[[i]] = tmp  
}

rm(tmp)

#Stochastic simulation******************************************************************************************************************************************************************
setwd("/Users/ijima/Desktop/Future_projection/projection_cons_F1214/res_cv")

for (i in 1:500) {
  
  d = list(int_n=int_n[[i]], #Set initial population number
           f=f1214, 　　　　 #Set constant Fishing mortality (define selectivity)
           catb=catb1014,　　#Set target catch weight
           maa=d4$maa,       #Set natural mortality at age
           waa=d4$waa,       #Set weight at age
           mat=d4$mat,       #Set maturity at age
           Rzero=d5$Rzero,   #Set R0
           SBzero=d5$SBzero, #Set SB0
           h=d5$h,           #Set Steepness
           sigmar=d5$sigmar, #Set sigma R
           rho=d6,           #Set autocorrelation parameter
           age=d4$age,       #Set maximum Age
           gender=d4$gender, #Set gender type 1:one gender, 2:two gender
           imax=1000,        #Set iteration number
           tmax=10,          #Set projectoin year
           recfun=2,         #Set type of recruitment uncertainty 1:normal 2:autocorrelation
           manage=1          #Set management scenario 1:Constant F, 2:Constant catch
  )

  set.seed(66)
  #assign(paste("res", i, sep=""), ssfcpp(d))
  res_list = ssfcpp(d)
  save(res_list, file=paste("res", i,"." ,"RData", sep=""))
}

rm(res_list)

#Plot figures**************************************************************************************************************************************************************************
#Calculation biological reference point
Dynamic_Bzero = myreplist$Dynamic_Bzero
SBF0 = Dynamic_Bzero %>% filter(Yr>=2012&Yr<=2014) %>%summarise(mean(SPB_nofishing)) 
SBF20 = 0.2*data.frame(SBF0)$mean.SPB_nofishing.
SBF30 = 0.3*data.frame(SBF0)$mean.SPB_nofishing.
SBF40 = 0.4*data.frame(SBF0)$mean.SPB_nofishing.
SBF50 = 0.5*data.frame(SBF0)$mean.SPB_nofishing.

#Female spawning biomass (Second quarter)
datfile = dir()
cnt = grep("", datfile)

tmp1 = 0
for (i in 1:length(cnt)){
  load(datfile[i])
  tmp1 = cbind(tmp1, data.frame(t(res_list$SB)))
}

tmp2 = tmp1[,-1]

prob16 = t(data.table(tmp2[2,])) %>% data.frame() %>% filter(.<=SBF20) %>% group_by(.) %>% summarize(count = n())
prob17 = t(data.table(tmp2[3,])) %>% data.frame() %>% filter(.<=SBF20) %>% group_by(.) %>% summarize(count = n())
prob18 = t(data.table(tmp2[4,])) %>% data.frame() %>% filter(.<=SBF20) %>% group_by(.) %>% summarize(count = n())
prob19 = t(data.table(tmp2[5,])) %>% data.frame() %>% filter(.<=SBF20) %>% group_by(.) %>% summarize(count = n())
prob20 = t(data.table(tmp2[6,])) %>% data.frame() %>% filter(.<=SBF20) %>% group_by(.) %>% summarize(count = n())
prob21 = t(data.table(tmp2[7,])) %>% data.frame() %>% filter(.<=SBF20) %>% group_by(.) %>% summarize(count = n())
prob22 = t(data.table(tmp2[8,])) %>% data.frame() %>% filter(.<=SBF20) %>% group_by(.) %>% summarize(count = n())
prob23 = t(data.table(tmp2[9,])) %>% data.frame() %>% filter(.<=SBF20) %>% group_by(.) %>% summarize(count = n())
prob24 = t(data.table(tmp2[10,])) %>% data.frame() %>% filter(.<=SBF20) %>% group_by(.) %>% summarize(count = n())
prob25 = t(data.table(tmp2[11,])) %>% data.frame() %>% filter(.<=SBF20) %>% group_by(.) %>% summarize(count = n())

prob=data.frame(prob16,prob17,prob18,prob19,prob20,prob21,prob22,prob23,prob24,prob25)

#Each run of future projection
#each_run = tmp2 %>% mutate(Yr=rep(2012:2042)) 
#each_run = melt(each_run, id=c("Yr"))
#each_run = each_run %>% rename(Run_no=variable,SpawnBio=value) %>% filter(Yr>=2012)

#Result of SS3
timeseries = myreplist$timeseries
ss_res = timeseries %>% filter(Seas==2,Yr>=1993) %>% select(Yr,SpawnBio)

tmp1 = myreplist$derived_quants
res_tile = tmp1[c(3:25),] %>% select(Value,StdDev) %>% mutate(Yr=rep(1993:2015),high=Value+1.96*StdDev,low=Value-1.96*StdDev)

rm(tmp1)

#Mean value of future projection
mean_run = data.frame(t(tmp2)) %>% summarise_each(funs(mean))
mean_run = data.frame(t(mean_run)) %>% mutate(Yr=rep(2015:2025)) %>% rename(SpawnBio=t.mean_run.) %>% filter(Yr>=2015)
mean_run[1,1] = res_tile[23,1]
#95% tile of future projection
#tmp3 = data.frame(t(tmp2)) %>% summarise_each(funs(quantile(., 0.05)))
#tmp3 = data.frame(t(tmp3))
#tmp4 = data.frame(t(tmp2)) %>% summarise_each(funs(quantile(., 0.95)))
#tmp4 = data.frame(t(tmp4))
#P_tile = tmp3 %>% mutate(Yr=rep(2012:2042),high=tmp4$t.tmp4.) %>% rename(low=t.tmp3.) %>% filter(Yr>=2012)
#rm(tmp1,tmp2,tmp3,tmp4)

tmp5 = data.frame(t(tmp2)) %>% summarise_each(funs(sd))
tmp5 = data.frame(t(tmp5))
tmp6 = data.frame(t(tmp2)) %>% summarise_each(funs(mean))
tmp6 = data.frame(t(tmp6))
tmp7 = cbind(tmp5,tmp6)
P_tile = tmp7 %>% mutate(Yr=rep(2015:2025), high=t.tmp6.+1.96*t.tmp5., low=t.tmp6.-1.96*t.tmp5.) %>% filter(Yr>=2015)
rm(tmp2,tmp5,tmp6,tmp7)

#Plot
g1 = ggplot()
#g1 = g1 + geom_line(data=each_run, aes(x=Yr, y=SpawnBio, group=Run_no),colour="#7f878f", size=0.2)
g1 = g1 + geom_ribbon(data=P_tile, aes(x=Yr, ymin=low, ymax=high),fill="#ff2800", alpha=0.3)
g1 = g1 + geom_ribbon(data=res_tile, aes(x=Yr, ymin=low, ymax=high),fill="#00a1bf", alpha=0.3)
g1 = g1 + geom_line(data=ss_res, aes(x=Yr, y=SpawnBio))
g1 = g1 + geom_line(data=mean_run,aes(x=Yr, y=SpawnBio), colour="#ff2800")
g1 = g1 + annotate("text", label="a)",x=1994,y=200000,size=4)
g1 = g1 + geom_hline(yintercept=SBF20,linetype="dashed")
g1 = g1 + annotate("text", label="0.2*SB[currentF==0]", x=1998,y=45000,size=3, parse=TRUE)
g1 = g1 + xlab("Year") + ylab("Spawning biomass (mt)")
g1

#Total catch
tmp1 = 0
for (i in 1:length(cnt)){
  load(datfile[i])
  tmp1 = cbind(tmp1, data.frame(t(res_list$C)))
}

tmp2 = tmp1[,-1]

Yr=rep(2015:2025)
tc = cbind(Yr,tmp2)
tc = melt(tc, id=c("Yr"))
tc = tc %>% rename(Run_no=variable,Catch=value) %>% filter(Yr>2015,Yr<=2025)
#tc2 = tc %>% filter(Yr==2015)
#mean(tc2$Catch)

g2 = ggplot(data=tc, aes(x=factor(Yr),y=Catch)) + geom_boxplot()
g2 = g2 + geom_hline(yintercept=catb1014,colour="#ff2800")
g2 = g2 + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))
g2 = g2 + annotate("text", label="b)",x=factor(2016),y=160000,size=4)
g2 = g2 + scale_y_continuous(breaks=seq(0,180000,by=40000),limits=c(30000,180000))
g2 = g2 + xlab("Year") + ylab("Total catch (mt)")
g2

grid.arrange(g1,g2,nrow=2)

#F multiplier
#fm = data.frame(t(res$fmulti)) %>% mutate(Yr=rep(2012:2042))
#fm = melt(fm, id=c("Yr"))
#fm = fm %>% rename(Run_no=variable,F_multiplier=value) %>% filter(Yr<=2041)

#g2 = ggplot(data=fm, aes(x=factor(Yr),y=F_multiplier)) + geom_boxplot()
#g2 = g2 + geom_hline(yintercept=1.0,colour="#ff2800")
#g2 = g2 + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))
#g2 = g2 + xlab("Year") + ylab("F multiplier")
#g2

#F at age 
F_at_age = data.frame(faa=f1214,age=rep(0:15,2),Gender=c(rep("Female",16),rep("Male",16)))

g3 = ggplot(data=F_at_age, aes(x=age,y=faa,colour=Gender)) 
g3 = g3 + geom_line()
g3 = g3 + ylim(0,0.3)
g3 = g3 + theme(legend.justification = c(1, 0), legend.position = c(1, 0))
g3 = g3 + theme(legend.background = element_rect(fill="transparent",colour=NA))
g3 = g3 + xlab("Year") + ylab(expression(paste("Instantaneous Fishing Mortality (yr"^"-1",")")))
g3
