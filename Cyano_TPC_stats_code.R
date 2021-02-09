# load packages
library(broom)
library(lubridate)
library(dplyr)
library(segmented)
library(bbmle)
library(ggplot2)
library(akima)
library(mleTools)
library(reshape2)
library(mgcv)
library(tidyr)
library(tidyverse)
library(broom)
library(gridExtra)
library(reshape2)

# load data 
##output from Cyano code/growth rate analysis for all species
dat<-read.csv("~/Data/Cyano_growth_data.csv") #output from Cyano_TPC_code
dat<-dat[with(dat,order(Acclim.temp)),]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#### Effect size differences ####

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

##overperformance function
opfunc <- function(x,y) {
  ((((x-y))/abs(y))*100)
}

# mean over replicates
mu.mean <- plyr::ddply(dat, c("Species","Nutrient.tx","Acclim.temp","Treatment","Temp.tx","Acute.temp","Acute.temp.id"), plyr::summarise, mu_avg = mean(mu), mu_se = sd(mu, na.rm=TRUE)/  
                          sqrt(length(mu[!is.na(mu)])))
mu.mean <- mu.mean[with(mu.mean,order(Species,Nutrient.tx,-mu_avg)),]

# anovas
caov<-aov(mu~factor(Acclim.temp), data=dat)
tab1<-summary(caov)

# fully acclimated only
acc <- subset(mu.mean,Acclim.temp==Acute.temp)

# cold acclimated only 
cold <- subset(mu.mean,Acclim.temp==15.66)

# % of total over performed
opdat <- merge(cold, acc, by=c('Species', 'Nutrient.tx', 'Acute.temp', 'Acute.temp.id'))
opdat <- mutate(opdat, op=opfunc(mu_avg.x,mu_avg.y))
table(sign(opdat$op))
47/(47+8+17) # over performance in 65% of cases

# max performance cold > acclimated
opdat2 <- subset(opdat, op>0)
op.max <- opdat2[with(opdat2,order(Species, Nutrient.tx, Acute.temp)),] # generates Table S3

# overperformance at moderate temp with adequate nutrients
op.mod <- subset(opdat, Acute.temp==22.09 & Nutrient.tx=="NORM")
opm.avg <- summarize(op.mod,mean=mean(op.mod$op),se=sd(op, na.rm=TRUE)/sqrt(length(op[!is.na(op)])))

# under performance hot < acclimated
hot <- subset(mu.mean,Temp.tx=="Hot")
updat <- merge(hot, acc, by=c('Species', 'Nutrient.tx', 'Acute.temp', 'Acute.temp.id'))
updat <- mutate(updat, up=opfunc(mu_avg.x,mu_avg.y))
table(sign(updat$up)) 
54/(54+8+10) # over performance in 75% of cases

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#### GAM model & CI comparisons ####

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# REPEAT THE FOLLOWING FOR EACH SPECIES AND EACH NUTRIENT CONDITION (listed below)

head(dat)
dat1_S<-dat%>% filter(Species=='A. flos-aquae' & Nutrient.tx=='NORM')
#dat1_S<-dat%>% filter(Species=='A. flos-aquae' & Nutrient.tx=='LOW')
#dat1_S<-dat%>% filter(Species=='M. aeruginosa (non-toxic)' & Nutrient.tx=='NORM')
#dat1_S<-dat%>% filter(Species=='M. aeruginosa (non-toxic)' & Nutrient.tx=='LOW')
#dat1_S<-dat%>% filter(Species=='M. aeruginosa (toxic)' & Nutrient.tx=='NORM')
#dat1_S<-dat%>% filter(Species=='M. aeruginosa (toxic)' & Nutrient.tx=='LOW')
#dat1_S<-dat%>% filter(Species=='P. foveolarum' & Nutrient.tx=='NORM')
#dat1_S<-dat%>% filter(Species=='P. foveolarum' & Nutrient.tx=='LOW')
head(dat1_S)

ggplot(dat1_S,aes(x=Acute.temp,y=mu))+
  geom_point(aes(colour=Temp.tx))

# reconstitute data set accounting for double-use of cold + acc and hot + acc
dat1_S_a <- subset(dat1_S, Temp.tx=="Cold" | Temp.tx=="Hot")
dat1_S_b <- subset(dat1_S, Acclim.temp==Acute.temp)
dat1_S_b$Temp.tx<-'Acc'
dat1_S<-rbind(dat1_S_a,dat1_S_b)

ggplot(dat1_S,aes(x=Acute.temp,y=mu))+
  geom_point(aes(colour=Temp.tx))+
  theme_bw()

###################################

### Question: does acclimation status affect the relationship between growth rate and temperature?

# Fit GAMs:

# one intercept and curve for all data
gm1<-gam(mu~s(Acute.temp,k=4),data=dat1_S)
summary(gm1)

# separate intercepts but same curve shapes by Temp.tx:
gm2<-gam(mu~Temp.tx+s(Acute.temp,k=4),data=dat1_S)
summary(gm2)

# separate intercepts and different curve shapes by Temp.tx:
gm3<-gam(mu~Temp.tx+s(Acute.temp,by=Temp.tx,k=4),data=dat1_S)
summary(gm3)

# Compare GAMs:

# AIC based model comparison
AICtab(gm1,gm2,gm3) #used to generate Table S2

# Visualize GAMs for comparison:

newdat<-data.frame(Acute.temp=seq(15,45,0.1))
pds_gm1<-predict(gm1,newdata=newdat,se.fit = T)
newdat$fit<-pds_gm1$fit
newdat$se<-pds_gm1$se.fit

newdat2<-expand.grid(Acute.temp=seq(15,45,0.1),Temp.tx=factor(c('Acc','Cold','Hot')))
pds_gm2<-predict(gm2,newdata=newdat2,se.fit = T)
newdat2$fit<-pds_gm2$fit
newdat2$se<-pds_gm2$se.fit

ggplot(dat1_S,aes(x=Acute.temp,y=mu))+
  geom_point(aes(colour=Temp.tx))+
  geom_ribbon(data=newdat,aes(y=fit,ymin=fit-1.96*se,ymax=fit+1.96*se),fill='gray',alpha=0.2)+
  geom_line(data=newdat,aes(y=fit))+
  geom_ribbon(data=newdat2,aes(y=fit,ymin=fit-1.96*se,ymax=fit+1.96*se,fill=Temp.tx),alpha=0.2)+
  geom_line(data=newdat2,aes(y=fit,colour=Temp.tx))+
  theme_bw()


###################################

### Question: How often does growth rate of cold/hot acclimated populations differ from perfectly acclimated?

newdat3<-unique(dat1_S[,c('Acute.temp','Temp.tx')])

pds<-predict(gm3,newdata=newdat3,se.fit=T)
newdat3$lw_ci95<-pds$fit-1.96*pds$se.fit
newdat3$up_ci95<-pds$fit+1.96*pds$se.fit

m1<-melt(newdat3,id.vars=c('Acute.temp','Temp.tx'))
head(m1)

d1<-dcast(m1,Acute.temp~Temp.tx+variable,value.var='value')

# function to check for overlaps in CIs:
# strategy: determine if intervals do NOT overlap
testCIs<-function(x_lw,x_up,y_lw,y_up){
  if(y_up < x_lw | y_lw > x_up){
    res<-1
  }else{
    res<-0
  }
  return(res)
}

# run this function for all Acute temperature treatments.
i<-1
cold_test<-rep(NA,nrow(d1))
hot_test<-rep(NA,nrow(d1))
for(i in 1:nrow(d1)){
  #print(i)
  cold_test[i]<-testCIs(d1$Acc_lw_ci95[i],d1$Acc_up_ci95[i],d1$Cold_lw_ci95[i],d1$Cold_up_ci95[i])
  hot_test[i]<-testCIs(d1$Acc_lw_ci95[i],d1$Acc_up_ci95[i],d1$Hot_lw_ci95[i],d1$Hot_up_ci95[i])
}

d1$Cold_test<-cold_test
d1$Hot_test<-hot_test
d1

# X% of experimental temperature treatments at which the cold acclimated population had growth rates with 95% CI that did not overlap with the 95% CI of the fully acclimated growth rates.
sum(d1$Cold_test)/(nrow(d1))

# X% of experimental temperature treatments at which the hot acclimated population had growth rates with 95% CI that did not overlap with the 95% CI of the fully acclimated growth rates.
sum(d1$Hot_test)/(nrow(d1))



