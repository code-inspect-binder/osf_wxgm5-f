#==== Experiment 1 Analyses ====

library(readxl)
library(lme4)
library(car)
library(emmeans)
library(ggplot2)
library(gridExtra)

data <- read.table('Exp1_data.csv', header=T, sep=',')

### Exclude contextual precursor trials for main analyses
dat <- subset(data, Distort %in% c(-12,-6,0,6,12))

### Data cleaning
dat$SD_Z <- scale(dat$Orig_SD)
sum(abs(dat$SD_Z) > 3)   

dat1 <- subset(dat, abs(SD_Z) < 3)   
dat1$Repr_BPM <- 60000/dat1$Orig_IBI      #Reproduced tempo (bpm)
dat1$Act_BPM <- 60000/dat1$Actual         #Actual tempo (bpm)
dat1$ratio <- dat1$Repr_BPM/dat1$Act_BPM  #Ratio of reproduced to actual tempo
sum(dat1$ratio > 1.5 | dat1$ratio < 0.75)

dat2 <- subset(dat1, ratio <= 1.5 & ratio >= 0.75)

dat2$Prct_Bias <- ((dat2$Repr_BPM - dat2$Act_BPM)/dat2$Act_BPM)*100         #Reproduction bias (% of actual)
dat2$Prct_Abs_Bias <- (abs(dat2$Repr_BPM - dat2$Act_BPM)/dat2$Act_BPM)*100  #Absolute error (% of actual)

by_subj <- aggregate(cbind(Prct_Bias, Prct_Abs_Bias, Orig_SD) ~ 
                       GMSI_AE*GMSI_PA*GMSI_MT*GMSI_SA*GMSI_E*GMSI_Gen*GroupC*Subject,
                     data=dat2,
                     FUN=mean)
by_subj$Prct_Bias_Z <- scale(by_subj$Prct_Bias)
by_subj$Prct_Abs_Bias_Z <- scale(by_subj$Prct_Abs_Bias)
by_subj$Orig_SD_Z <- scale(by_subj$Orig_SD)
#Participant #44 was an outlier!

#Exclude outlier participant
dat2 <- subset(dat2, Subject != 44)

### Analyzing number of tries during the practice routine 
prac <- aggregate(cbind(T1_Tries, T2_Tries) ~ Song*Subject, data=dat2, FUN=mean)
prac$total_tries <- prac$T1_Tries + prac$T2_Tries
prac$fail <- prac$T1_Tries==6 | prac$T2_Tries==6
prac1 <- aggregate(fail ~ Subject, data=prac, FUN=sum)
table(prac1$fail)
prac2 <- aggregate(total_tries ~ Subject, data=prac, FUN=mean)
by_subj1 <- subset(by_subj, Subject != 44)
by_subj1$total_tries <- prac2$total_tries
round(cor(by_subj1[,c(1:6,9:11,15)]), 2)
cor.test(by_subj1$total_tries, by_subj1$Prct_Abs_Bias)
cor.test(by_subj1$total_tries, by_subj1$GMSI_Gen)

### Correlations between GMSI subscales and performance measures
round(cor(subset(by_subj, Subject != 44)[,c(1:6,9:11)]), 2)

### Linear mixed-effects model for reproduction bias 
mod_full <- lmer(Prct_Bias ~ ContextC*(GMSI_Gen_Z + GroupC*OrderC*Distort) +
                   (ContextC | Subject), data=dat2)
summary(mod_full)
Anova(mod_full, type=3, test='Chisq')   # Wald tests

### Null model (random intercepts only)
mod0 <- lmer(Prct_Bias ~ 1 + (1|Subject), data=dat2)

### Test improvement of full model over null model
anova(mod0, mod_full)

### Simple slopes: Model for testing effect of Context at high GMSI
mod_1 <- lmer(Prct_Bias ~ ContextC*(GMSI_Gen_hi + GroupC*OrderC*Distort) +
                (ContextC | Subject), data=dat2)
summary(mod_1)
Anova(mod_1, type=3, test='Chisq')

### Simple slopes: Model for testing effect of Context at low GMSI
mod_2 <- lmer(Prct_Bias ~ ContextC*(GMSI_Gen_lo + GroupC*OrderC*Distort) +
                (ContextC | Subject), data=dat2)
summary(mod_2)
Anova(mod_2, type=3, test='Chisq')

### Results after exluding songs that were unsuccesfully practiced 
dat3 <- subset(dat2, T1_Tries < 6 & T2_Tries < 6)
mod_full1 <- lmer(Prct_Bias ~ ContextC*(GMSI_Gen_Z + GroupC*OrderC*Distort) +
                    (ContextC | Subject), data=dat3)
summary(mod_full1)
Anova(mod_full1, type=3, test='Chisq')

### Get predicted means from model for plotting interaction
means1 <- as.data.frame(emmeans(mod_full, 
                                "ContextC", 
                                by="GMSI_Gen_Z", 
                                at=list(ContextC=c(-1,1),
                                        GMSI_Gen_Z=c(-1, 0, 1)), 
                                mode='asymp'))

names(means1) <- c('Context','GMSI-General','Mean','SE','df','LL','UL')
means1$Context <- factor(means1$Context, levels=c(-1,1), labels=c('Slow','Fast'))
means1$`GMSI-General` <- factor(means1$`GMSI-General`, levels=c(-1,0,1), labels=c('Low (-1 SD)',
                                                                                  'Average',
                                                                                  'High (+1 SD)'))
means1$UL.SE <- means1$Mean + means1$SE
means1$LL.SE <- means1$Mean - means1$SE

### Plotting mean bias across all tempo precursor levels and GMSI x Context plot

#Need to first clean the full dataset with contextual precursors included
dat <- data
dat$SD_Z <- scale(dat$Orig_SD)
dat1 <- subset(dat, abs(SD_Z) < 3)   
dat1$Repr_BPM <- 60000/dat1$Orig_IBI      
dat1$Act_BPM <- 60000/dat1$Actual         
dat1$ratio <- dat1$Repr_BPM/dat1$Act_BPM  
dat2 <- subset(dat1, ratio <= 1.5 & ratio >= 0.75)
dat2$Prct_Bias <- ((dat2$Repr_BPM - dat2$Act_BPM)/dat2$Act_BPM)*100         
dat2$Prct_Abs_Bias <- (abs(dat2$Repr_BPM - dat2$Act_BPM)/dat2$Act_BPM)*100  
dat2 <- subset(dat2, Subject != 44)

dat2_1 <- aggregate(Prct_Bias ~ Distort*Context*Subject, data=dat2, FUN=mean)
dat2_2 <- aggregate(Prct_Bias ~ Distort*Context, data=dat2_1, FUN=mean)
dat2_2_se <- aggregate(Prct_Bias ~ Distort*Context, data=dat2_1, 
                       FUN=function(x) sd(x)/sqrt(length(x)))
dat2_2$UL <- dat2_2$Prct_Bias + dat2_2_se$Prct_Bias
dat2_2$LL <- dat2_2$Prct_Bias - dat2_2_se$Prct_Bias

dat2_2$Context <- factor(dat2_2$Context)

dat2_2$Type <- factor(ifelse(dat2_2$Distort %in% c(-12,-6,0,6,12), 'Target', 'Contextual'))

dat2_2$Ctxt_Type <- character(nrow(dat2_2))
dat2_2$Ctxt_Type[dat2_2$Context=='Fast' & 
                   dat2_2$Type=='Contextual'] <- 'Fast, Context'
dat2_2$Ctxt_Type[dat2_2$Context=='Fast' & 
                   dat2_2$Type=='Target'] <- 'Fast, Target'
dat2_2$Ctxt_Type[dat2_2$Context=='Slow' & 
                   dat2_2$Type=='Contextual'] <- 'Slow, Context'
dat2_2$Ctxt_Type[dat2_2$Context=='Slow' & 
                   dat2_2$Type=='Target'] <- 'Slow, Target'
dat2_2$Ctxt_Type <- factor(dat2_2$Ctxt_Type,
                           levels=c('Slow, Target','Slow, Context',
                                    'Fast, Target','Fast, Context'))

p1 <- ggplot(data=dat2_2, aes(x=Distort, y=Prct_Bias, group=Context)) +
  geom_point(aes(x=Distort, y=Prct_Bias, shape=Ctxt_Type), col='black', cex=3) +
  geom_errorbar(aes(x=Distort, ymin=LL, ymax=UL), width=1.5) +
  geom_line(aes(x=Distort, y=Prct_Bias)) +
  theme_bw(base_size = 14) +
  theme(panel.grid=element_blank(), legend.title = element_blank(),
        axis.text.x = element_text(angle=45, size=8),
        legend.position=c(0.8,.825), legend.text=element_text(size=8)) +
  geom_abline(aes(slope=0, intercept=0), lty=2) +
  scale_x_continuous(breaks=c(-30,-24,-18,-12,-6,0,6,12,18,24,30), expand=c(0,0)) +
  scale_y_continuous(expand=c(0,.75), limits=c(-5,10)) +
  ylab('Reproduction Bias (%)') +
  xlab('Percent Tempo Distortion') +
  scale_shape_manual(values=c(1,2,16,17)) +
  ggtitle('A')

p2 <- ggplot(data=means1, aes(x=Context, y=Mean, group=`GMSI-General`, colour=`GMSI-General`, shape=`GMSI-General`)) +
  geom_point(cex=3) +
  geom_errorbar(aes(x=Context, ymin=LL.SE, ymax=UL.SE), width=.05) +
  geom_line() +
  theme_bw(base_size=14) +
  theme(panel.grid=element_blank(), legend.position=c(.73,.825),
        legend.text=element_text(size=8)) +
  geom_abline(aes(slope=0, intercept=0), lty=2) +
  scale_y_continuous(expand=c(0,.75), limits=c(-5,10)) +
  ylab('Reproduction Bias (%)') +
  xlab('Context') +
  scale_color_manual(name='GMSI-General', values=c('red','blue','black'),
                     labels=c('Low (-1 SD)', 'Average', 'High (+1 SD)')) +
  scale_shape_manual(name='GMSI-General', values=c(15,16,17),
                     labels=c('Low (-1 SD)', 'Average', 'High (+1 SD)')) +
  ggtitle('B')

grid.arrange(p1, p2, nrow=1, ncol=2)

### Testing models for each song separately

ELC <- subset(dat2, Song=='Electric Chapel')
mod_ELC <- lmer(Prct_Bias ~ ContextC*(GMSI_Gen_Z + OrderC*Distort) +
                   (ContextC | Subject), data=ELC)
summary(mod_ELC)
Anova(mod_ELC, type=3, test='Chisq')  
#Sig. effect of Precursor

RKB <- subset(dat2, Song=='Rock Bottom')
mod_RKB <- lmer(Prct_Bias ~ ContextC*(GMSI_Gen_Z + OrderC*Distort) +
                  (ContextC | Subject), data=RKB)
summary(mod_RKB)
Anova(mod_RKB, type=3, test='Chisq')  
#Sig. effect of GMSI-General (marginal Context x GMSI interaction)

BTB <- subset(dat2, Song=='The Boogie That Be')
mod_BTB <- lmer(Prct_Bias ~ ContextC*(GMSI_Gen_Z + OrderC*Distort) +
                  (ContextC | Subject), data=BTB)
summary(mod_BTB)
Anova(mod_BTB, type=3, test='Chisq')  
#***Model failed to converge***
#Sig. effect of Precursor 

PKF <- subset(dat2, Song=='Poker Face')
mod_PKF <- lmer(Prct_Bias ~ ContextC*(GMSI_Gen_Z + OrderC*Distort) +
                  (ContextC | Subject), data=PKF)
summary(mod_PKF)
Anova(mod_PKF, type=3, test='Chisq')  
#No significant effects

MNI <- subset(dat2, Song=='My Name Is')
mod_MNI <- lmer(Prct_Bias ~ ContextC*(GMSI_Gen_Z + OrderC*Distort) +
                  (ContextC | Subject), data=MNI)
summary(mod_MNI)
Anova(mod_MNI, type=3, test='Chisq')  
#***Identifiability issue***
#Marginal effects of Context and Precursor 

LGS <- subset(dat2, Song=='Lets Get It Started')
mod_LGS <- lmer(Prct_Bias ~ ContextC*(GMSI_Gen_Z + OrderC*Distort) +
                  (ContextC | Subject), data=LGS)
summary(mod_LGS)
Anova(mod_LGS, type=3, test='Chisq')  
#No significant effects

### Plotting song-specific effects

#Need to first clean the full dataset with contextual precursors included
dat <- data
dat$SD_Z <- scale(dat$Orig_SD)
dat1 <- subset(dat, abs(SD_Z) < 3)   
dat1$Repr_BPM <- 60000/dat1$Orig_IBI      
dat1$Act_BPM <- 60000/dat1$Actual         
dat1$ratio <- dat1$Repr_BPM/dat1$Act_BPM  
dat2 <- subset(dat1, ratio <= 1.5 & ratio >= 0.75)
dat2$Prct_Bias <- ((dat2$Repr_BPM - dat2$Act_BPM)/dat2$Act_BPM)*100         
dat2$Prct_Abs_Bias <- (abs(dat2$Repr_BPM - dat2$Act_BPM)/dat2$Act_BPM)*100  
dat2 <- subset(dat2, Subject != 44)

dat2_1 <- aggregate(Prct_Bias ~ Distort*Song*Context*Subject, data=dat2, FUN=mean)
dat2_2 <- aggregate(Prct_Bias ~ Distort*Song*Context, data=dat2_1, FUN=mean)
dat2_2_se <- aggregate(Prct_Bias ~ Distort*Song*Context, data=dat2_1, 
                       FUN=function(x) sd(x)/sqrt(length(x)))
dat2_2$UL <- dat2_2$Prct_Bias + dat2_2_se$Prct_Bias
dat2_2$LL <- dat2_2$Prct_Bias - dat2_2_se$Prct_Bias

dat2_2$Context <- factor(dat2_2$Context)

dat2_2$Type <- factor(ifelse(dat2_2$Distort %in% c(-12,-6,0,6,12), 'Target', 'Contextual'))

dat2_2$Ctxt_Type <- character(nrow(dat2_2))
dat2_2$Ctxt_Type[dat2_2$Context=='Fast' & 
                   dat2_2$Type=='Contextual'] <- 'Fast, Context'
dat2_2$Ctxt_Type[dat2_2$Context=='Fast' & 
                   dat2_2$Type=='Target'] <- 'Fast, Target'
dat2_2$Ctxt_Type[dat2_2$Context=='Slow' & 
                   dat2_2$Type=='Contextual'] <- 'Slow, Context'
dat2_2$Ctxt_Type[dat2_2$Context=='Slow' & 
                   dat2_2$Type=='Target'] <- 'Slow, Target'
dat2_2$Ctxt_Type <- factor(dat2_2$Ctxt_Type,
                           levels=c('Slow, Target','Slow, Context',
                                    'Fast, Target','Fast, Context'))

dat2_2$Song <- factor(dat2_2$Song, 
                      levels=c('Poker Face','Lets Get It Started','My Name Is',
                               'Electric Chapel','The Boogie That Be','Rock Bottom'))

p1 <- ggplot(data=dat2_2, aes(x=Distort, y=Prct_Bias, group=Context)) +
  geom_point(aes(x=Distort, y=Prct_Bias, shape=Ctxt_Type), col='black', cex=3) +
  geom_errorbar(aes(x=Distort, ymin=LL, ymax=UL), width=1.5) +
  geom_line(aes(x=Distort, y=Prct_Bias)) +
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


