#Calculation of N at age***************************************************************************************************************************************
get_naa = function(myreplist) {
  
  natage = myreplist$natage
  startyr = myreplist$startyr
  endyr = myreplist$endyr
  
  tmp1 = natage %>% filter(`Beg/Mid`=="B",Yr>=startyr) %>% select(-Area,-Bio_Pattern,-BirthSeas,-SubMorph,-Morph,-Time,-Era,-`Beg/Mid`)
  tmp2 = melt(tmp1, c("Yr","Gender","Seas"))
  tmp3 = tmp2 %>% rename(age=variable,pop_n=value) %>% filter(Seas==2,Yr==endyr) %>% arrange(Gender,Seas) %>% select(pop_n)
  naa = tmp2 %>% rename(age=variable,pop_n=value) %>% filter(Yr<=endyr) %>% arrange(Gender,Seas)
  int_n = tmp3$pop_n
  
  rm(tmp1,tmp2,tmp3)
  
  list(naa=naa,int_n=int_n)
  
}


#Calculation of Catch at age**********************************************************************************************************************************
get_caa = function(myreplist) {
 
  startyr = myreplist$startyr
   
  catage = myreplist$catage
  tmp = catage[,-c(1,2,4,5,6,9,10)]
  tmp = tmp %>% filter(Yr>=startyr) %>% group_by(Yr,Gender,Seas) %>% summarise_each(funs(sum),everything())
  tmp = melt(tmp, id=c("Yr","Gender","Seas"))
  caa = tmp %>% rename(age=variable,cat.n=value) %>% arrange(Gender,Seas)
  
  rm(tmp)
  
  return(caa=caa)
  
}


#Calculation of F at age**************************************************************************************************************************************
get_faa = function(naa,caa,myreplist) {
  
  startyr = myreplist$startyr 
  
  tmp0 = naa %>% filter(Seas==2) %>% mutate(age=as.numeric(age))
  tmp1 = caa %>% filter(Seas==1) %>% select(Yr,Gender,age,cat.n) %>% mutate(Yr=Yr-1)
  tmp2 = caa %>% filter(Seas>=2) %>% group_by(Yr,Gender,age) %>% summarise(cat.n=sum(cat.n))
  tmp3 = bind_rows(tmp1,tmp2)
  tmp3 = tmp3 %>% group_by(Yr,Gender,age) %>% summarise(cat.n=sum(cat.n)) %>% mutate(age=as.numeric(age))
  tmp4 = full_join(tmp0,tmp3,by=c("Yr","Gender","age")) %>% filter(Yr>=startyr )
  
  #faa was calculated by Newton method
  x = list()
  for (i in 1:length(tmp4$Yr)) {
    fx = function (x) (x/(x+0.3))*(1-exp(-x-0.3))*tmp4[i,5]-tmp4[i,6]
    res = uniroot(fx,c(0,1))
    x[i] = res$root
  }
  
  tmp4$faa = unlist(x)
  
  #rm(tmp0,tmp1,tmp2,tmp3)
  
  return(tmp4)
  
}



#M at age and Maturity at age*********************************************************************************************************************************
get_maa = function(myreplist) {
  
  endgrowt = myreplist$endgrowth #Wt_Beg=waa in season 2
  endgrowt = endgrowt %>% filter(Seas==2)
  bio = endgrowt %>% select(Wt_Beg,Age_Mat) %>% mutate(Age_Mat=ifelse(Age_Mat<=0, 0, Age_Mat))
  waa = bio$Wt_Beg
  mat = bio$Age_Mat
  age = myreplist$accuage
  gender = myreplist$nsexes
  #maa = rep(0.3,(age+1)*gender)
  maa = c(1.36, 0.56, 0.45, 0.48, 0.48, 0.48, 0.48, 0.48, 0.48, 0.48, 0.48, 0.48, 0.48, 0.48, 0.48, 0.48,
          1.36, 0.56, 0.45, 0.39, 0.39, 0.39, 0.39, 0.39, 0.39, 0.39, 0.39, 0.39, 0.39, 0.39, 0.39, 0.39)
  
  list(waa=waa,mat=mat,age=age,gender=gender,maa=maa)
  
}


#Choose R0****************************************************************************************************************************************************
get_rec = function(myreplist) {
  
  parameters = myreplist$parameters
  tmp1 = parameters %>% filter(Label=="SR_LN(R0)") %>% select(Value) %>% mutate(Value=exp(Value))
  tmp2 = parameters %>% filter(Label=="SR_BH_steep") %>% select(Value)
  tmp3 = parameters %>% filter(Label=="SR_sigmaR") %>% select(Value)
  
  Rzero = tmp1$Value
  h = tmp2$Value
  sigmar = tmp3$Value
  tmp4 = myreplist$Dynamic_Bzero %>% filter(Era=="VIRG") %>% select(SPB) 
  SBzero = tmp4$SPB
  
  rm(tmp1,tmp2,tmp3,tmp4)
  
  list(Rzero=Rzero,h=h,sigmar=sigmar,SBzero=SBzero)
  
}


#Recruitment autocorrelation**********************************************************************************************************************************
get_ac = function(myreplist) {
  
  startyr = myreplist$startyr
  
  tmp1 = myreplist$recruitpars %>% filter(Yr>=startyr)
  tmp2 = acf(tmp1$Value)
  ac = tmp2$acf[2]
  return(ac)
  
}


#Total catch weight*******************************************************************************************************************************************
get_cat_b = function(myreplist) {

  tmp1 = myreplist$sprseries
  tmp2 = tmp1[,-c(2,3)]
  cat_b = tmp2 %>% select(Year,Retain_Catch_B)
  return(cat_b)
  
  rm(tmp1,tmp2)
  
}

#Get CV of SSB************************************************************************************************************************************************
get_cv = function(myreplist) {
  tmp = myreplist$derived_quants %>% filter(LABEL=="SPB_2012") %>% mutate(cv=StdDev/Value)
  cv = tmp$cv
  return(cv)
  
  rm(tmp)
}
