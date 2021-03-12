#load packages
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
library(multcompView)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#### Load and prep data ####

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# load toxin and RFU data
dat0<-read.csv("~/Dropbox/Phyto Acclimation/Manuscripts/Cyanotoxins/Data & code/Repository/Data/Cyano_toxin_data.csv")

# calculate standard curve
dat.std <- dat0 %>% filter(Nutrients=="Std") #standards only
m1<-nls(Conc~SSasymp(Abs_450, yf, y0, log_alpha), data=dat.std)
summary(m1)

abs<-seq(-0.5,2.5,.01)
m1_pred<-predict(m1,list(Abs_450 = abs))
plot(dat.std$Conc~dat.std$Abs)
lines(m1_pred,abs, col= "red", lwd=3,lty=3)

# remove stds and predict toxin concentration for each replicate (avg over duplicates)
dat1<-dat0[!is.na(dat0$RFU_raw_Final),] 
dat<-dat1 %>% group_by(Starting_Cond, Acute_temp, Nutrients, Replicate,RFU_raw_Final,RFU_raw_initial,Dilution_fact) %>%
  summarise(Abs_450=mean(Abs_450))%>%
  mutate(Mc_raw=predict(m1,list(Abs_450=Abs_450)))%>%
  mutate(Mc_adj=Dilution_fact*predict(m1,list(Abs_450=Abs_450)))%>%
ungroup()

# attach temps and temp id's
tgb.col1<-data.frame(Starting_Cond=seq(1,9),temp=c(15.7, 18.5, 22.1, 25.4, 28.6, 32.0, 35.7, 39.1, 42.6))
dat<-merge(dat,tgb.col1,by=c('Starting_Cond'))
names(dat)[ncol(dat)]<-'Acclim_temp'
dat$Acute_temp.id<-5
head(dat)

# estimate growth rate (1/day) over 5 days
dat <- dat %>% group_by(Starting_Cond,Acclim_temp,Acute_temp.id,Acute_temp,Nutrients,Replicate) %>% 
  mutate(mu = (log(RFU_raw_Final/RFU_raw_initial))/5)

dat$id <- (paste(dat$Acclim_temp,dat$Nutrients,sep="_"))

# microcystin produced per cell (ppb/cell) 

#slopes
slope1 <- data.frame(Estimate=c(3930,1000,7310,1000,1970,1000), Int=rep(0,times=6),
                     id=c("15.7_Low","15.7_Normal", "28.6_Low","28.6_Normal","39.1_Low","39.1_Normal"))

#merge dataframes, and estimate new initial and final abundances
dat<-merge(dat,slope1,by=c("id"))%>% 
  mutate(counts_initial_adj = RFU_raw_initial * Estimate+Int) %>%
  mutate(counts_final_adj = RFU_raw_Final * Estimate+Int)

#re-estimate x0 and x1 from new estimates
dat <- dat %>% group_by(Starting_Cond,Acclim_temp,Acute_temp.id,Acute_temp,Nutrients, Replicate,RFU_raw_Final) %>% 
  mutate(time=5) %>% #add in the time elapsed
  mutate(M0=0) %>% #add in starting toxin concentration
  mutate(loss=0.09) %>% # add in loss rate that describes loss rate of toxins based on a half life of ~3weeks 
  mutate(x0 =(exp(loss*time)*Mc_adj-M0)* (loss+mu) /  (-1 + exp((loss*time)+(mu*time))*counts_initial_adj)       )%>% #add in per cell toxin production
  mutate(x1= 1.56831*Mc_adj*(0.09+mu)/ (-1+exp(0.45+(5*mu))*counts_initial_adj     )   ) #should be same as above (edited)

# re-label 
# nutrients
dat$Medium<-ifelse(dat$Nutrients=="Low","Limited","COMBO")
# acclimation 
dat$Tx.label<-ifelse(dat$Acclim_temp==15.7, "Cold", ifelse(dat$Acclim_temp==39.1, "Hot", "Fully-acclimated"))
# treatment 
dat$id<-paste(dat$Acclim_temp,dat$Nutrients,sep = "_")
dat$treatment<-rep(c("A", "B", "C", "D", "E","F"), each=6) 

dat<-as.data.frame(dat)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#### Results & stats ####

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
## Growth rates (1/day) ##
#stats test & set-up for manual Tukey plot
mu<-lm(mu~Nutrients*factor(Acclim_temp),data=dat)
summary(mu)
summary(aov(mu~Nutrients*factor(Acclim_temp),data=dat))

m3<-aov(mu~as.factor(Acclim_temp),data=dat[dat$ Nutrients=='Low',])
m4<-aov(mu~as.factor(Acclim_temp),data=dat[dat$ Nutrients=='Normal',])
TukeyHSD(m3)
TukeyHSD(m4)

lab2<-data.frame("Letters"=c("α","a","α","b","β","c"),"id"=c("15.66_Low","15.66_Normal","28.63_Low","28.63_Normal","39.09_Low","39.09_Normal"))

names(lab2) <- c('Letters','id')
yvalue2 <- aggregate(dat$mu, list(dat$id), data=dat, quantile, probs=.75)  
final2 <- data.frame(lab2, yvalue2[,2])
final2$Acclim_temp<-rep(c(15.7, 28.6, 39.1), each=2) 
final2$Tx.label<-rep(c("Cold", "Fully-acclimated", "Hot"), each=2) 
final2$Medium<-rep(c("Limited","COMBO"), each=1) 
names(final2)<-c("letters","id","mu","Acclim_temp","Tx.label","Medium")

#plot 1
g1<-ggplot(data.frame(dat),aes(x=as.factor(Acclim_temp),y=mu,fill=Tx.label, linetype=Medium,color=Tx.label))+
  geom_boxplot()+
  scale_linetype_manual(values=c("solid","dashed"),labels=c("High","Low"))+
  labs(tag = "(a)")+
  scale_fill_manual(values=c("#3399FF","gray60","#FF6666"))+
  scale_color_manual(values=c("#004D9A","gray30","#CD0000"))+
  theme_bw()+
  theme(legend.position="none",axis.title.x = element_text(size = rel(1.2)),axis.title.y = element_text(size = rel(1.2)),
        axis.text.y = element_text(size = rel(1.2)),axis.text.x = element_text(size = rel(1.2)))+
  labs(linetype="Nutrient condition",color="Acclimation temperature", fill="Acclimation temperature")+
  geom_hline(yintercept=0,linetype="solid")+
  scale_x_discrete('Acclimation temperature (°C)')+
  scale_y_continuous('Growth rate (1/day)',limits=c(-0.2,0.8),breaks=c(-0.25,0,0.25,0.5,0.75))+
  geom_text(data = final2,aes(x=as.factor(Acclim_temp), y=mu,label=letters),vjust=-1.5, hjust=2)
g1

## Total MC (ppb) ##
#stats test & set-up for manual Tukey plot
mc<-lm(Mc_adj~Nutrients*factor(Acclim_temp),data=dat)
summary(mc)
summary(aov(Mc_adj~Nutrients*factor(Acclim_temp),data=dat))

m1<-aov(Mc_adj~as.factor(Acclim_temp),data=dat[dat$ Nutrients=='Low',])
m2<-aov(Mc_adj~as.factor(Acclim_temp),data=dat[dat$ Nutrients=='Normal',])
TukeyHSD(m1)
TukeyHSD(m2)

lab1<-data.frame("Letters"=c("α","a","β","b","β","c"),"id"=c("15.66_Low","15.66_Normal","28.63_Low","28.63_Normal","39.09_Low","39.09_Normal"))

names(lab1) <- c('Letters','id')
yvalue1 <- aggregate(dat$Mc_adj, list(dat$id), data=dat, quantile, probs=.75)  
final1 <- data.frame(lab1, yvalue1[,2])
final1$Acclim_temp<-rep(c(15.7, 28.6, 39.1), each=2) 
final1$Tx.label<-rep(c("Cold", "Fully-acclimated", "Hot"), each=2) 
final1$Medium<-rep(c("Limited","COMBO"), each=1) 
names(final1)<-c("letters","id","Mc_adj","Acclim_temp","Tx.label","Medium")

#plot 2
g2<-ggplot(data.frame(dat),aes(x=as.factor(Acclim_temp),y=Mc_adj,fill=Tx.label, linetype=Medium,color=Tx.label))+
  geom_boxplot()+
  scale_linetype_manual(values=c("solid","dashed"),labels=c("High","Low"))+
  labs(tag = "(b)")+
  scale_fill_manual(values=c("#3399FF","gray60","#FF6666"))+
  scale_color_manual(values=c("#004D9A","gray30","#CD0000"))+
  theme_bw()+
  theme(legend.position="bottom",legend.title=element_text(size = rel(1.5)),legend.text=element_text(size = rel(1.2)),axis.title.x = element_text(size = rel(1.2)),axis.title.y = element_text(size = rel(1.2)),
        axis.text.y = element_text(size = rel(1.2)),axis.text.x = element_text(size = rel(1.2)))+
  labs(linetype="Nutrient condition",color="Acclimation temperature", fill="Acclimation temperature")+
  geom_hline(yintercept=0,linetype="solid")+
  scale_x_discrete('Acclimation temperature (°C)')+
  scale_y_continuous('Total MC (ug/L)')+
  geom_text(data = final1,aes(x=as.factor(Acclim_temp), y=Mc_adj,label=letters),vjust=-1.5, hjust=2)
g2

## MC produced per cell (ppb/cell) ##
#stats test & set-up for manula Tukey plot
mcc<-lm(x0~Nutrients*factor(Acclim_temp),data=dat)
summary(mcc)
summary(aov(x0~Nutrients*factor(Acclim_temp),data=dat))

m5<-aov(x0~as.factor(Acclim_temp),data=dat[dat$ Nutrients=='Low',])
m6<-aov(x0~as.factor(Acclim_temp),data=dat[dat$ Nutrients=='Normal',])
TukeyHSD(m5)
TukeyHSD(m6)

lab3<-data.frame("Letters"=c("α","a","β","b","β","a"),"id"=c("15.66_Low","15.66_Normal","28.63_Low","28.63_Normal","39.09_Low","39.09_Normal"))

names(lab3) <- c('Letters','id')
yvalue3 <- aggregate(dat$x0, list(dat$id), data=dat, quantile, probs=.75)  
final3 <- data.frame(lab3, yvalue3[,2])
final3$Acclim_temp<-rep(c(15.7, 28.6, 39.1), each=2) 
final3$Tx.label<-rep(c("Cold", "Fully-acclimated", "Hot"), each=2) 
final3$Medium<-rep(c("Limited","COMBO"), each=1) 
names(final3)<-c("letters","id","x0","Acclim_temp","Tx.label","Medium")

# plot 3
g3<-ggplot(data.frame(dat),aes(x=as.factor(Acclim_temp),y=x0,fill=as.factor(Tx.label), linetype=Medium,color=as.factor(Tx.label)))+  #note y is different from (Mc_adj/counts_final)
  geom_boxplot()+
  scale_linetype_manual(values=c("solid","dashed"),labels=c("High","Low"))+
  scale_fill_manual(values=c("#3399FF","gray60","#FF6666"))+
  scale_color_manual(values=c("#004D9A","gray30","#CD0000"))+
  theme_bw()+
  theme(legend.position="none",axis.title.x = element_text(size = rel(1.2)),axis.title.y = element_text(size = rel(1.2)),
        axis.text.y = element_text(size = rel(1.2)),axis.text.x = element_text(size = rel(1.2)))+
  labs(linetype="Nutrient condition",color="Acclimation temperature", fill="Acclimation temperature", tag="(c)")+
  geom_hline(yintercept=0,linetype="solid")+
  scale_x_discrete('Acclimation temperature (°C)')+
  scale_y_continuous('MC (ug/L)/cell/day')+
  geom_text(data = final3,aes(x=as.factor(Acclim_temp), y=x0,label=letters),vjust=-1.5, hjust=2)
g3

#Figure 2

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend <- get_legend(g2)

g4<-grid.arrange(g1 + theme(legend.position = "none"),g2+ theme(legend.position = "none"),g3+ theme(legend.position = "none"),bottom=legend,nrow=1)
g4
ggsave(filename = '~/Desktop/Toxin_figs.png',g4,width = 12, height = 5, units = "in") 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#### Effect size differences ####

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# check out means
sum<-dat%>%
  group_by(Starting_Cond,Nutrients) %>%
  summarise(gr=mean(mu),
            gr_sd=sd(mu),
            gr_se=sd(mu)/sqrt(6)
            #mc=mean(x1), 
            #mc_sd=sd(x1),
            #mc_se=sd(x1)/sqrt(6)
            )
sum

sum<-dat%>%
  group_by(Starting_Cond,Nutrients) %>%
  summarise(tx=mean(x0),
            tx_sd=sd(x0),
            tx_se=sd(x0)/sqrt(6),
  )
sum

# average growth rate
iavg <- plyr::ddply(dat, c("Nutrients","Acclim_temp","Starting_Cond"), plyr::summarise, RFU_avg = mean(RFU_raw_initial), RFU_se = sd(RFU_raw_initial, na.rm=TRUE)/  
                      sqrt(length(RFU_raw_initial[!is.na(RFU_raw_initial)])))

# over performance function
opfunc <- function(x,y) {
  ((((x-y))/abs(y))*100)
}

# average under/overperformance in MC production
tavg <- plyr::ddply(dat, c("Nutrients","Acclim_temp","Starting_Cond"), plyr::summarise, mc_avg = mean(Mc_adj), mc_se = sd(Mc_adj, na.rm=TRUE)/  
                      sqrt(length(Mc_adj[!is.na(Mc_adj)])))
tavgb <- plyr::ddply(dat, c("Acclim_temp","Starting_Cond"), plyr::summarise, mc_avg = mean(Mc_adj), mc_se = sd(Mc_adj, na.rm=TRUE)/  
                       sqrt(length(Mc_adj[!is.na(Mc_adj)])))

opfunc(tavg$mc_avg[4],tavg$mc_avg[5]) #nutrient-replete only; cold, acc = 51.4% increase
opfunc(tavg$mc_avg[1],tavg$mc_avg[2]) #nutrient-limited only; cold, acc = 79.81% increase
opfunc(tavg$mc_avg[6],tavg$mc_avg[5]) #nutrient-replete only; hot, acc = 37.3% decrease
opfunc(tavg$mc_avg[3],tavg$mc_avg[2]) #nutrient-limited only; hot, acc = 24.3% decrease

# average under/overperformance in MC production PER CELL
tavg2 <- plyr::ddply(dat, c("Nutrients","Acclim_temp","Starting_Cond"), plyr::summarise, mc_avg = mean(x0), mc_se = sd(x0, na.rm=TRUE)/  
                      sqrt(length(x0[!is.na(x0)])))
opfunc(tavg2$mc_avg[4],tavg2$mc_avg[6]) #nutrient-replete; cold, hot = 9% increase
opfunc(tavg2$mc_avg[4],tavg2$mc_avg[5]) #nutrient-replete; cold, acc = 31% increase
opfunc(tavg2$mc_avg[1],tavg2$mc_avg[2]) #nutrient-limited; cold, acc = 472.6% increase
opfunc(tavg2$mc_avg[5],tavg2$mc_avg[6]) #nutrient-replete; acc, hot = 17% decrease
opfunc(tavg2$mc_avg[3],tavg2$mc_avg[2]) #nutrient-limited; hot, acc = 6% decrease

