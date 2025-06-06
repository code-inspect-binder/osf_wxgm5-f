#==== Experiment 3 Analyses ====

library(readxl)
library(lme4)
library(car)
library(emmeans)
library(ggplot2)
library(gridExtra)

data <- read.table('Exp3_data.csv', header=T, sep=',')

#Exclude trials with contextual precursors for main analyses 
dat <- subset(data, Distort %in% c(-12,-6,0,6,12))

dat$Orig_SD_Z <- scale(dat$Orig_SD)    #Original tapped IBI SD
dat$Prec_SD_Z <- scale(dat$Ctxt_SD)    #Precursor tapped IBI SD
sum(abs(dat$Orig_SD_Z) > 3 | abs(dat$Prec_SD_Z) > 3)
dat1 <- subset(dat, abs(Orig_SD_Z) < 3 & abs(Prec_SD_Z) < 3) 

dat1$Repr_Orig_BPM <- 60000/dat1$Orig_IBI                        #Reproduced original tempo (bpm)
dat1$Repr_Prec_BPM <- 60000/dat1$Ctxt_IBI                        #Reproduced precursor tempo (bpm)
dat1$Orig_BPM <- 60000/dat1$Actual                               #Actual original tempo (bpm)
dat1$Prec_BPM <- dat1$Orig_BPM + .01*dat1$Distort*dat1$Orig_BPM  #Actual precursor tempo (bpm)
dat1$Orig_ratio <- dat1$Repr_Orig_BPM/dat1$Orig_BPM
dat1$Prec_ratio <- dat1$Repr_Prec_BPM/dat1$Prec_BPM
sum(dat1$Orig_ratio > 1.5 | dat1$Orig_ratio < 0.75 |
      dat1$Prec_ratio > 1.5 | dat1$Prec_ratio < 0.75) 
dat2 <- subset(dat1, Orig_ratio <= 1.5 & Orig_ratio >= 0.75 & Prec_ratio <= 1.5 & Prec_ratio >= 0.75)

#Signed bias and absolute error scores (% of actual)
dat2$Prct_Orig_Bias <- ((dat2$Repr_Orig_BPM - dat2$Orig_BPM)/dat2$Orig_BPM)*100
dat2$Prct_Abs_Orig_Bias <- (abs(dat2$Repr_Orig_BPM - dat2$Orig_BPM)/dat2$Orig_BPM)*100
dat2$Prct_Prec_Bias <- ((dat2$Repr_Prec_BPM - dat2$Prec_BPM)/dat2$Prec_BPM)*100
dat2$Prct_Abs_Prec_Bias <- (abs(dat2$Repr_Prec_BPM - dat2$Prec_BPM)/dat2$Prec_BPM)*100

by_subj <- aggregate(cbind(Prct_Orig_Bias, Prct_Abs_Orig_Bias, 
                           Prct_Prec_Bias, Prct_Abs_Prec_Bias,
                           Orig_SD, Ctxt_SD) ~ 
                       AE*PA*MT*SA*E*Gen*FamOrderC*Context*Subject,
                     data=dat2,
                     FUN=mean)
by_subj$Prct_Orig_Bias_Z <- scale(by_subj$Prct_Orig_Bias)
by_subj$Prct_Abs_Orig_Bias_Z <- scale(by_subj$Prct_Abs_Orig_Bias)
by_subj$Orig_SD_Z <- scale(by_subj$Orig_SD)
by_subj$Prct_Prec_Bias_Z <- scale(by_subj$Prct_Prec_Bias)
by_subj$Prct_Abs_Prec_Bias_Z <- scale(by_subj$Prct_Abs_Prec_Bias)
by_subj$Ctxt_SD_Z <- scale(by_subj$Ctxt_SD)
#Participant #35 was an outlier!

### Correlations between GMSI subscales and performance measures
round(cor(subset(by_subj, Subject != 35)[,c(1:6,10:15)]), 2)

dat2 <- subset(dat2, Subject != 35)

### Analyzing number of tries during the practice routine 
prac <- aggregate(cbind(T1_Tries, T2_Tries) ~ Song*Subject, data=dat2, FUN=mean)
prac$total_tries <- prac$T1_Tries + prac$T2_Tries
prac$fail <- prac$T1_Tries==6 | prac$T2_Tries==6
prac1 <- aggregate(fail ~ Subject, data=prac, FUN=sum)
table(prac1$fail)
prac2 <- aggregate(total_tries ~ Subject, data=prac, FUN=mean)
by_subj1 <- subset(by_subj, Subject != 35)
by_subj1$total_tries <- prac2$total_tries
round(cor(by_subj1[,c(1:6,10:15,22)]), 2)
cor.test(by_subj1$total_tries, by_subj1$Prct_Abs_Prec_Bias)
cor.test(by_subj1$total_tries, by_subj1$Prct_Abs_Orig_Bias)
cor.test(by_subj1$total_tries, by_subj1$Gen)


#Reproduced precursor tempo expressed as a deviation from original tempo (rescale to aid model convergence)  
dat2$Repr_Prec_Dev <- ((dat2$Repr_Prec_BPM-dat2$Orig_BPM)/dat2$Orig_BPM)*100
dat2$Repr_Prec_Dev1 <- scale(dat2$Repr_Prec_Dev, center=F)   

#Need to also rescale Distort to aid model convergence
dat2$Distort1 <- scale(dat2$Distort, center=F)

### Which is more predictive of original reproduction bias? The precursor on each trial or the reproduced precursor tempo?
mod1a <- lmer(Prct_Orig_Bias ~ Distort1 + (Distort1|Subject), data=dat2)
summary(mod1a)
BIC(mod1a)
Anova(mod1a, type=3, test='Chisq')

mod1b <- lmer(Prct_Orig_Bias ~ Repr_Prec_Dev1 + (Repr_Prec_Dev1|Subject), data=dat2)
summary(mod1b)
BIC(mod1b)
Anova(mod1b, type=3, test='Chisq')

### Null model (intercepts only)
mod0 <- lmer(Prct_Orig_Bias ~ 1 + (1|Subject), data=dat2)

### Linear mixed effects model for original reproduction bias
contrasts(dat2$Context) <- 'contr.sum'
mod_full <- lmer(Prct_Orig_Bias ~ Repr_Prec_Dev1 + Context*(GMSI_Gen_Z + FamiliarC*FamOrderC) +
                   (Repr_Prec_Dev1 + FamiliarC||Subject), data=dat2,
                 control=lmerControl(optCtrl=list(maxfun=2000)))
Anova(mod_full, type=3, test='Chisq')
summary(mod_full)
BIC(mod_full)  

### Improvement of full model over null model
anova(mod0, mod_full)

#Try letting reproduced PT interact with Context and GMSI
mod_full1 <- lmer(Prct_Orig_Bias ~ Context*(GMSI_Gen_Z*Repr_Prec_Dev1 + FamiliarC*FamOrderC) +
                    (Repr_Prec_Dev1 + FamiliarC||Subject), data=dat2,
                  control=lmerControl(optCtrl=list(maxfun=2000)))
BIC(mod_full1) #higher BIC than original model

#Try letting reproduced PT interact with design variables
mod_full2 <- lmer(Prct_Orig_Bias ~ Context*(GMSI_Gen_Z + Repr_Prec_Dev1*FamiliarC*FamOrderC) +
                    (Repr_Prec_Dev1 + FamiliarC||Subject), data=dat2,
                  control=lmerControl(optCtrl=list(maxfun=2000)))
BIC(mod_full2) #higher BIC than original model

### Results were not affected by exluding songs that were unsuccesfully practiced 
dat3 <- subset(dat2, T1_Tries < 6 & T2_Tries < 6)
mod_full3 <- lmer(Prct_Orig_Bias ~ Repr_Prec_Dev1 + Context*(GMSI_Gen_Z + FamiliarC*FamOrderC) +
                    (Repr_Prec_Dev1 + FamiliarC||Subject), data=dat3,
                  control=lmerControl(optCtrl=list(maxfun=2000)))
summary(mod_full3)
Anova(mod_full3, type=3, test='Chisq')

emmeans(mod_full, "Context", mode='asymp')

###Get predicted means from model for GMSI x Context plot
means1 <- as.data.frame(emmeans(mod_full, 
                                "Context", 
                                by="GMSI_Gen_Z", 
                                at=list(GMSI_Gen_Z=c(-1, 0, 1)), 
                                mode='asymp'))


names(means1) <- c('Context','GMSI-General','Mean','SE','df','LL','UL')
means1$Context <- factor(means1$Context, levels=c('slow','fast'))
means1$`GMSI-General` <- factor(means1$`GMSI-General`, levels=c(-1,0,1), labels=c('Low (-1 SD)',
                                                                                  'Average',
                                                                                  'High (+1 SD)'))
means1$UL.SE <- means1$Mean + means1$SE
means1$LL.SE <- means1$Mean - means1$SE

### Testing models for each song separately

ELC <- subset(dat2, Song=='Electric Chapel')
mod_ELC <- lmer(Prct_Orig_Bias ~ Repr_Prec_Dev1 + Context*(GMSI_Gen_Z + FamOrderC) +
                  (Repr_Prec_Dev1 | Subject), data=ELC,
                control=lmerControl(optCtrl=list(maxfun=2000)))
summary(mod_ELC)
Anova(mod_ELC, type=3, test='Chisq')  
#No significant effects

RKB <- subset(dat2, Song=='Rock Bottom')
mod_RKB <- lmer(Prct_Orig_Bias ~ Repr_Prec_Dev1 + Context*(GMSI_Gen_Z + FamOrderC) +
                  (Repr_Prec_Dev1 | Subject), data=RKB,
                control=lmerControl(optCtrl=list(maxfun=2000)))
summary(mod_RKB)
Anova(mod_RKB, type=3, test='Chisq')  
#Sig. effect of reproduced PT and marginal effect of GMSI-General

BTB <- subset(dat2, Song=='The Boogie That Be')
mod_BTB <- lmer(Prct_Orig_Bias ~ Repr_Prec_Dev1 + Context*(GMSI_Gen_Z + FamOrderC) +
                  (Repr_Prec_Dev1 | Subject), data=BTB,
                control=lmerControl(optCtrl=list(maxfun=2000)))
summary(mod_BTB)
Anova(mod_BTB, type=3, test='Chisq')  
#Marginal effect of reproduced PT 

PKF <- subset(dat2, Song=='Poker Face')
mod_PKF <- lmer(Prct_Orig_Bias ~ Repr_Prec_Dev1 + Context*(GMSI_Gen_Z + FamOrderC) +
                  (Repr_Prec_Dev1 | Subject), data=PKF,
                control=lmerControl(optCtrl=list(maxfun=2000)))
summary(mod_PKF)
Anova(mod_PKF, type=3, test='Chisq')  
#Significant effect of Context (assim.)

MNI <- subset(dat2, Song=='My Name Is')
mod_MNI <- lmer(Prct_Orig_Bias ~ Repr_Prec_Dev1 + Context*(GMSI_Gen_Z + FamOrderC) +
                  (Repr_Prec_Dev1 | Subject), data=MNI,
                control=lmerControl(optCtrl=list(maxfun=2000)))
summary(mod_MNI)
Anova(mod_MNI, type=3, test='Chisq')  
#Significant effect of reproduced PT
#Marginal effect of Context (assim.)
#Significant effects of GMSI-General and Familiarity Order

LGS <- subset(dat2, Song=='Lets Get It Started')
mod_LGS <- lmer(Prct_Orig_Bias ~ Repr_Prec_Dev1 + Context*(GMSI_Gen_Z + FamOrderC) +
                  (Repr_Prec_Dev1 | Subject), data=LGS,
                control=lmerControl(optCtrl=list(maxfun=2000)))
summary(mod_LGS)
Anova(mod_LGS, type=3, test='Chisq')
#Significant effects of Context (assim.), reproduced PT, and GMSI-General

### Plotting song-specific effects

#Need to first clean the full dataset with contextual precursors included
dat <- data
dat$Orig_SD_Z <- scale(dat$Orig_SD)    
dat$Prec_SD_Z <- scale(dat$Ctxt_SD)    
dat1 <- subset(dat, abs(Orig_SD_Z) < 3 & abs(Prec_SD_Z) < 3) 
dat1$Repr_Orig_BPM <- 60000/dat1$Orig_IBI                        
dat1$Repr_Prec_BPM <- 60000/dat1$Ctxt_IBI                        
dat1$Orig_BPM <- 60000/dat1$Actual                               
dat1$Prec_BPM <- dat1$Orig_BPM + .01*dat1$Distort*dat1$Orig_BPM  
dat1$Orig_ratio <- dat1$Repr_Orig_BPM/dat1$Orig_BPM
dat1$Prec_ratio <- dat1$Repr_Prec_BPM/dat1$Prec_BPM
dat2 <- subset(dat1, Orig_ratio <= 1.5 & Orig_ratio >= 0.75 & Prec_ratio <= 1.5 & Prec_ratio >= 0.75)
dat2$Prct_Orig_Bias <- ((dat2$Repr_Orig_BPM - dat2$Orig_BPM)/dat2$Orig_BPM)*100
dat2$Prct_Abs_Orig_Bias <- (abs(dat2$Repr_Orig_BPM - dat2$Orig_BPM)/dat2$Orig_BPM)*100
dat2$Prct_Prec_Bias <- ((dat2$Repr_Prec_BPM - dat2$Prec_BPM)/dat2$Prec_BPM)*100
dat2$Prct_Abs_Prec_Bias <- (abs(dat2$Repr_Prec_BPM - dat2$Prec_BPM)/dat2$Prec_BPM)*100
dat2 <- subset(dat2, Subject != 35)

dat2_1 <- aggregate(Prct_Orig_Bias ~ Distort*Song*Context*Subject, data=dat2, FUN=mean)

dat2_2 <- aggregate(Prct_Orig_Bias ~ Distort*Song*Context, data=dat2_1, FUN=mean)
dat2_2_se <- aggregate(Prct_Orig_Bias ~ Distort*Song*Context, data=dat2_1, 
                       FUN=function(x) sd(x)/sqrt(length(x)))
dat2_2$UL <- dat2_2$Prct_Orig_Bias + dat2_2_se$Prct_Orig_Bias
dat2_2$LL <- dat2_2$Prct_Orig_Bias - dat2_2_se$Prct_Orig_Bias

dat2_2$Context <- factor(dat2_2$Context)

dat2_2$Type <- factor(ifelse(dat2_2$Distort %in% c(-12,-6,0,6,12), 'Target', 'Contextual'))

dat2_2$Ctxt_Type <- character(nrow(dat2_2))
dat2_2$Ctxt_Type[dat2_2$Context=='fast' & 
                   dat2_2$Type=='Contextual'] <- 'Fast, Context'
dat2_2$Ctxt_Type[dat2_2$Context=='fast' & 
                   dat2_2$Type=='Target'] <- 'Fast, Target'
dat2_2$Ctxt_Type[dat2_2$Context=='slow' & 
                   dat2_2$Type=='Contextual'] <- 'Slow, Context'
dat2_2$Ctxt_Type[dat2_2$Context=='slow' & 
                   dat2_2$Type=='Target'] <- 'Slow, Target'
dat2_2$Ctxt_Type <- factor(dat2_2$Ctxt_Type,
                           levels=c('Slow, Target','Slow, Context',
                                    'Fast, Target','Fast, Context'))

dat2_2$Song <- factor(dat2_2$Song, levels=c('Poker Face','Lets Get It Started','My Name Is',
                                            'Electric Chapel','The Boogie That Be','Rock Bottom'))

p1 <- ggplot(data=dat2_2, aes(x=Distort, y=Prct_Orig_Bias, group=Context)) +
  geom_point(aes(x=Distort, y=Prct_Orig_Bias, shape=Ctxt_Type), col='black', cex=3) +
  geom_errorbar(aes(x=Distort, ymin=LL, ymax=UL), width=1.5) +
  geom_line(aes(x=Distort, y=Prct_Orig_Bias)) +
  theme_bw(base_size = 14) +
  theme(panel.grid=element_blank(), legend.title = element_blank(),
        axis.text.x = element_text(angle=45, size=8),
        legend.text=element_text(size=8)) +
  geom_abline(aes(slope=0, intercept=0), lty=2) +
  scale_x_continuous(breaks=c(-30,-24,-18,-12,-6,0,6,12,18,24,30), expand=c(0,0)) +
  ylab('Reproduction Bias (%)') +
  xlab('Percent Tempo Distortion') +
  scale_shape_manual(values=c(1,2,16,17)) +
  facet_wrap(~Song, nrow=2, ncol=3)
p1



### Plotting mean bias across all tempo precursor levels and GMSI x Context plot

#Need to first clean the full dataset with contextual precursors included
dat <- data
dat$Orig_SD_Z <- scale(dat$Orig_SD)    
dat$Prec_SD_Z <- scale(dat$Ctxt_SD)    
dat1 <- subset(dat, abs(Orig_SD_Z) < 3 & abs(Prec_SD_Z) < 3) 
dat1$Repr_Orig_BPM <- 60000/dat1$Orig_IBI                        
dat1$Repr_Prec_BPM <- 60000/dat1$Ctxt_IBI                        
dat1$Orig_BPM <- 60000/dat1$Actual                               
dat1$Prec_BPM <- dat1$Orig_BPM + .01*dat1$Distort*dat1$Orig_BPM  
dat1$Orig_ratio <- dat1$Repr_Orig_BPM/dat1$Orig_BPM
dat1$Prec_ratio <- dat1$Repr_Prec_BPM/dat1$Prec_BPM
dat2 <- subset(dat1, Orig_ratio <= 1.5 & Orig_ratio >= 0.75 & Prec_ratio <= 1.5 & Prec_ratio >= 0.75)
dat2$Prct_Orig_Bias <- ((dat2$Repr_Orig_BPM - dat2$Orig_BPM)/dat2$Orig_BPM)*100
dat2$Prct_Abs_Orig_Bias <- (abs(dat2$Repr_Orig_BPM - dat2$Orig_BPM)/dat2$Orig_BPM)*100
dat2$Prct_Prec_Bias <- ((dat2$Repr_Prec_BPM - dat2$Prec_BPM)/dat2$Prec_BPM)*100
dat2$Prct_Abs_Prec_Bias <- (abs(dat2$Repr_Prec_BPM - dat2$Prec_BPM)/dat2$Prec_BPM)*100
dat2 <- subset(dat2, Subject != 35)

dat2_1 <- aggregate(Prct_Orig_Bias ~ Distort*Context*Subject, data=dat2, FUN=mean)

dat2_2 <- aggregate(Prct_Orig_Bias ~ Distort*Context, data=dat2_1, FUN=mean)
dat2_2_se <- aggregate(Prct_Orig_Bias ~ Distort*Context, data=dat2_1, 
                       FUN=function(x) sd(x)/sqrt(length(x)))
dat2_2$UL <- dat2_2$Prct_Orig_Bias + dat2_2_se$Prct_Orig_Bias
dat2_2$LL <- dat2_2$Prct_Orig_Bias - dat2_2_se$Prct_Orig_Bias

dat2_2$Context <- factor(dat2_2$Context)

dat2_2$Type <- factor(ifelse(dat2_2$Distort %in% c(-12,-6,0,6,12), 'Target', 'Contextual'))

dat2_2$Ctxt_Type <- character(nrow(dat2_2))
dat2_2$Ctxt_Type[dat2_2$Context=='fast' & 
                   dat2_2$Type=='Contextual'] <- 'Fast, Context'
dat2_2$Ctxt_Type[dat2_2$Context=='fast' & 
                   dat2_2$Type=='Target'] <- 'Fast, Target'
dat2_2$Ctxt_Type[dat2_2$Context=='slow' & 
                   dat2_2$Type=='Contextual'] <- 'Slow, Context'
dat2_2$Ctxt_Type[dat2_2$Context=='slow' & 
                   dat2_2$Type=='Target'] <- 'Slow, Target'
dat2_2$Ctxt_Type <- factor(dat2_2$Ctxt_Type,
                           levels=c('Slow, Target','Slow, Context',
                                    'Fast, Target','Fast, Context'))

p1 <- ggplot(data=dat2_2, aes(x=Distort, y=Prct_Orig_Bias, group=Context)) +
  geom_point(aes(x=Distort, y=Prct_Orig_Bias, shape=Ctxt_Type), col='black', cex=3) +
  geom_errorbar(aes(x=Distort, ymin=LL, ymax=UL), width=1.5) +
  geom_line(aes(x=Distort, y=Prct_Orig_Bias)) +
  theme_bw(base_size = 14) +
  theme(panel.grid=element_blank(), legend.title = element_blank(),
        axis.text.x = element_text(angle=45, size=8),
        legend.position=c(0.2,.825), legend.text=element_text(size=8)) +
  geom_abline(aes(slope=0, intercept=0), lty=2) +
  scale_x_continuous(breaks=c(-30,-24,-18,-12,-6,0,6,12,18,24,30), expand=c(0,0)) +
  scale_y_continuous(expand=c(0,.75), limits=c(-5,15)) +
  ylab('Reproduction Bias (%)') +
  xlab('Percent Tempo Distortion') +
  scale_shape_manual(values=c(1,2,16,17)) +
  ggtitle('A')

p2 <- ggplot(data=means1, aes(x=Context, y=Mean, group=`GMSI-General`, colour=`GMSI-General`, shape=`GMSI-General`)) +
  geom_point(cex=3) +
  geom_errorbar(aes(x=Context, ymin=LL.SE, ymax=UL.SE), width=.05) +
  geom_line() +
  theme_bw(base_size=14) +
  theme(panel.grid=element_blank(), legend.position=c(.25,.825),
        legend.text=element_text(size=8)) +
  geom_abline(aes(slope=0, intercept=0), lty=2) +
  scale_y_continuous(expand=c(0,.75), limits=c(-5,15)) +
  ylab('Reproduction Bias (%)') +
  xlab('Context') +
  scale_color_manual(name='GMSI-General', values=c('red','blue','black'),
                     labels=c('Low (-1 SD)', 'Average', 'High (+1 SD)')) +
  scale_shape_manual(name='GMSI-General', values=c(15,16,17),
                     labels=c('Low (-1 SD)', 'Average', 'High (+1 SD)')) +
  ggtitle('B')

grid.arrange(p1, p2, nrow=1, ncol=2)
