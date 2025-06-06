#==== Experiment 2: Reproduction Task ====

library(readxl)
library(lme4)
library(car)
library(emmeans)
library(ggplot2)
library(gridExtra)

data <- read.table('Exp2_tapping.csv', header=T, sep=',')

#Convert to bpm
data$Repr_BPM <- 60000/data$IBI        #Reproduced tempo
data$Actual_BPM <- 60000/data$Actual   #Actual tempo

data$ContextC <- ifelse(data$Context=='fast', 1, -1)
data$OrderC <- ifelse(data$Order=='J1T2', 1, -1)

#Exclude contextual precursor trials for main analyses
dat <- subset(data, Distort %in% c(-12,-6,0,6,12))

### Data cleaning 
dat$SD_Z <- scale(dat$SD_IBI)
sum(abs(dat$SD_Z) > 3)
dat1 <- subset(dat, abs(SD_Z) < 3) 

dat1$ratio <- dat1$Repr_BPM/dat1$Actual_BPM
sum(dat1$ratio > 1.5 | dat1$ratio < 0.75) 
dat2 <- subset(dat1, ratio <= 1.5 & ratio >= 0.75)

dat2$Prct_Bias <- ((dat2$Repr_BPM - dat2$Actual_BPM)/dat2$Actual_BPM)*100
dat2$Prct_Abs_Bias <- (abs(dat2$Repr_BPM - dat2$Actual_BPM)/dat2$Actual_BPM)*100

by_subj <- aggregate(cbind(Prct_Bias, Prct_Abs_Bias, SD_IBI) ~ 
                       GMSI_AE*GMSI_PA*GMSI_MT*GMSI_SA*GMSI_E*GMSI_Gen*Context*Order*Subject,
                     data=dat2,
                     FUN=mean)
by_subj$Prct_Bias_Z <- scale(by_subj$Prct_Bias)
by_subj$Prct_Abs_Bias_Z <- scale(by_subj$Prct_Abs_Bias)
by_subj$Orig_SD_Z <- scale(by_subj$SD_IBI)

### Correlations between GMSI subscales and performance measures
round(cor(by_subj[,c(1:6,10:12)]), 2)

### Analyzing number of tries during the practice routine 
prac <- aggregate(Tries ~ Subject, data=dat2, FUN=mean)
sum(prac$Tries==6)
by_subj$Tries <- prac$Tries
round(cor(by_subj[,c(1:6,10:12,16)]), 2)
cor.test(by_subj$Tries, by_subj$Prct_Abs_Bias)

### Linear mixed effects model for reproduction bias
mod_full <- lmer(Prct_Bias ~ ContextC*(GMSI_Gen_Z + OrderC*Distort) + (1|Subject),
                 data=dat2)
summary(mod_full)
Anova(mod_full, type=3, test='Chisq')

### Null model (random intercepts only)
mod0 <- lmer(Prct_Bias ~ 1 + (1|Subject), data=dat2)

### Test improvement of full model over null model
anova(mod0, mod_full)

### Get predicted means from the model for GMSI x Context plot 
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

#Re-fit model excluding participants who failed to meet practice criterion
dat3 <- subset(dat2, !(Subject %in% c(1,2,15,23,32,54)))

mod_full1 <- lmer(Prct_Bias ~ ContextC*(GMSI_Gen_Z + OrderC*Distort) + (1|Subject),
                  data=dat3)
summary(mod_full1)
Anova(mod_full1, type=3, test='Chisq')

### Testing models for each song separately

ELC <- subset(dat2, Song=='Electric Chapel')
mod_ELC <- lmer(Prct_Bias ~ ContextC*(GMSI_Gen_Z + OrderC*Distort) + (1|Subject),
                data=ELC)
summary(mod_ELC)
Anova(mod_ELC, type=3, test='Chisq')  
#Sig. effect of Precursor, Context (assim.), and Task Order

RKB <- subset(dat2, Song=='Rock Bottom')
mod_RKB <- lmer(Prct_Bias ~ ContextC*(GMSI_Gen_Z + OrderC*Distort) + (1|Subject),
                data=RKB)
summary(mod_RKB)
Anova(mod_RKB, type=3, test='Chisq')  
#No sig. effects

BTB <- subset(dat2, Song=='The Boogie That Be')
mod_BTB <- lmer(Prct_Bias ~ ContextC*(GMSI_Gen_Z + OrderC*Distort) + (1|Subject),
                data=BTB)
summary(mod_BTB)
Anova(mod_BTB, type=3, test='Chisq')  
#Sig. effect of Precursor 

### Plotting song-specific effects

#Need to first clean the full dataset with contextual precursors included
dat <- data
dat$SD_Z <- scale(dat$SD_IBI)
dat1 <- subset(dat, abs(SD_Z) < 3) 
dat1$ratio <- dat1$Repr_BPM/dat1$Actual_BPM
dat2 <- subset(dat1, ratio <= 1.5 & ratio >= 0.75)
dat2$Prct_Bias <- ((dat2$Repr_BPM - dat2$Actual_BPM)/dat2$Actual_BPM)*100
dat2$Prct_Abs_Bias <- (abs(dat2$Repr_BPM - dat2$Actual_BPM)/dat2$Actual_BPM)*100

dat2_1 <- aggregate(Prct_Bias ~ Distort*Song*Context*Subject, data=dat2, FUN=mean)
dat2_2 <- aggregate(Prct_Bias ~ Distort*Song*Context, data=dat2_1, FUN=mean)
dat2_2_se <- aggregate(Prct_Bias ~ Distort*Song*Context, data=dat2_1, 
                       FUN=function(x) sd(x)/sqrt(length(x)))
dat2_2$UL <- dat2_2$Prct_Bias + dat2_2_se$Prct_Bias
dat2_2$LL <- dat2_2$Prct_Bias - dat2_2_se$Prct_Bias

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
  scale_y_continuous(limits=c(-10,20)) +
  ylab('Reproduction Bias (%)') +
  xlab('Percent Tempo Distortion') +
  scale_shape_manual(values=c(1,2,16,17)) +
  facet_wrap(~Song, nrow=1, ncol=3)
p1

### Plotting mean bias across all tempo precursor levels and GMSI x Context plot

#Need to first clean the full dataset with contextual precursors included
dat <- data
dat$SD_Z <- scale(dat$SD_IBI)
dat1 <- subset(dat, abs(SD_Z) < 3) 
dat1$ratio <- dat1$Repr_BPM/dat1$Actual_BPM
dat2 <- subset(dat1, ratio <= 1.5 & ratio >= 0.75)
dat2$Prct_Bias <- ((dat2$Repr_BPM - dat2$Actual_BPM)/dat2$Actual_BPM)*100
dat2$Prct_Abs_Bias <- (abs(dat2$Repr_BPM - dat2$Actual_BPM)/dat2$Actual_BPM)*100

dat2_1 <- aggregate(Prct_Bias ~ Distort*Context*Subject, data=dat2, FUN=mean)
dat2_2 <- aggregate(Prct_Bias ~ Distort*Context, data=dat2_1, FUN=mean)
dat2_2_se <- aggregate(Prct_Bias ~ Distort*Context, data=dat2_1, 
                       FUN=function(x) sd(x)/sqrt(length(x)))
dat2_2$UL <- dat2_2$Prct_Bias + dat2_2_se$Prct_Bias
dat2_2$LL <- dat2_2$Prct_Bias - dat2_2_se$Prct_Bias

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

p1 <- ggplot(data=dat2_2, aes(x=Distort, y=Prct_Bias, group=Context)) +
  geom_point(aes(x=Distort, y=Prct_Bias, shape=Ctxt_Type), col='black', cex=3) +
  geom_errorbar(aes(x=Distort, ymin=LL, ymax=UL), width=1.5) +
  geom_line(aes(x=Distort, y=Prct_Bias)) +
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
  theme(panel.grid=element_blank(), legend.position=c(.73,.825),
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

