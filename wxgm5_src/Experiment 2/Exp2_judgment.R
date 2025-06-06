#==== Experiment 2: Perceptual Judgment Task ====

library(readxl)
library(lme4)
library(car)
library(emmeans)
library(ggplot2)
library(gridExtra)

#==== Analyzing binary judgments ====

jdata <- read.table('Exp2_binaryjudge.csv', sep=',', header=T)

jdat <- subset(jdata, Distort %in% c(-12, -6, 0, 6, 12))
jdat$Subject <- factor(jdat$Subject)

#Participants 17 and 53 rated all targets as faster than original
jdat_all_1s <- aggregate(faster ~ Subject, data=jdat, FUN=prod)

#Participants 8 and 23 rated all targets as slower than original
jdat$slower <- 1-jdat$faster
jdat_all_0s <- aggregate(slower ~ Subject, data=jdat, FUN=prod)

jdat$cannot_fit <- numeric(nrow(jdat))
jdat$cannot_fit[jdat$Subject %in% c(17,53,8,23)] <- 1

#Fit logistic model to each individual's data to estimate PSEs
jdat$PSE <- 500
for(i in as.numeric(levels(jdat$Subject))){
  dat <- subset(jdat, Subject==i)
  if(dat$cannot_fit[1]==0){
    fit.glm <- glm(faster ~ Distort, family=binomial, data=dat)
    jdat$PSE[jdat$Subject==i] <- -coef(fit.glm)[1]/coef(fit.glm)[2]  #PSE is -b0/b1
  }
}

#Make sure all PSEs fall within [-12,12]
jdat$PSE[jdat$Subject==17 | jdat$Subject==53] <- -12.00
jdat$PSE[jdat$Subject==8 | jdat$Subject==23] <- 12.00
jdat$PSE[jdat$PSE > 12.00] <- 12.00
jdat$PSE[jdat$PSE < -12.00] <- -12.00

#One row per individual
jdat1 <- aggregate(PSE ~ Context*Order*GMSI_AE_Z*GMSI_PA_Z*GMSI_MT_Z*GMSI_SA_Z*GMSI_E_Z*GMSI_Gen_Z*Subject, 
                   data=jdat, FUN=mean)
jdat1$Context <- factor(jdat1$Context)
jdat1$Order <- factor(jdat1$Order)
contrasts(jdat1$Context) <- 'contr.sum'; contrasts(jdat1$Order) <- 'contr.sum'

### Regression model for PSEs 
mod.PSE <- lm(PSE ~ Context*(GMSI_Gen_Z + Order), data=jdat1)
summary(mod.PSE)

emmeans(mod.PSE, "Context", at=list(GMSI_Gen_Z=0))

cor.test(jdat1$GMSI_Gen_Z, abs(jdat1$PSE))

#==== Analyzing graded judgments ====

jdata2 <- read.table('Exp2_AllJudge.csv', sep=',', header=T)

jdat2 <- subset(jdata2, Distort %in% c(-12,-6,0,6,12))

jdat2$Context <- factor(jdat2$Context)
jdat2$Order <- factor(jdat2$Order)

### Linear mixed effects model for graded judgments
mod_full <- lmer(GradedJudge ~ Context*(GMSI_Gen_Z + Order*Distort) + (1|Subject), 
                 data=jdat2)
summary(mod_full)
Anova(mod_full, type=3, test='Chisq')

### Null (intercepts-only) model
mod0 <- lmer(GradedJudge ~ 1 + (1|Subject), 
             data=jdat2)
anova(mod0, mod_full)

### Refit the model with each participant's ratings on a common scale 
jdat2$Graded_rescaled <- jdat2$GradedJudge/ave(jdat2$GradedJudge, jdat2$Subject, FUN=sd)
mod_full1 <- lmer(Graded_rescaled ~ Context*(GMSI_Gen_Z + Order*Distort) + (1|Subject), 
                  data=jdat2)
summary(mod_full1)
Anova(mod_full1, type=3, test='Chisq')


#==== Plotting graded and binary judgments ====

dat1 <- aggregate(faster ~ Distort*Context*Subject, data=jdata, FUN=mean)
dat1_1 <- aggregate(faster ~ Distort*Context, data=dat1, FUN=mean)

dat1_1_se <- aggregate(faster ~ Distort*Context, data=dat1, 
                       FUN=function(x) sd(x)/sqrt(length(x)))
dat1_1$UL <- dat1_1$faster + dat1_1_se$faster
dat1_1$LL <- dat1_1$faster - dat1_1_se$faster

dat1_1$Context <- factor(dat1_1$Context, levels=c('slow','fast'), 
                         labels=c('Slow','Fast'))

dat1_1$Type <- factor(ifelse(dat1_1$Distort %in% c(-12,-6,0,6,12), 'Target', 'Contextual'))

dat1_1$Ctxt_Type <- character(nrow(dat1_1))
dat1_1$Ctxt_Type[dat1_1$Context=='Fast' & 
                   dat1_1$Type=='Contextual'] <- 'Fast, Context'
dat1_1$Ctxt_Type[dat1_1$Context=='Fast' & 
                   dat1_1$Type=='Target'] <- 'Fast, Target'
dat1_1$Ctxt_Type[dat1_1$Context=='Slow' & 
                   dat1_1$Type=='Contextual'] <- 'Slow, Context'
dat1_1$Ctxt_Type[dat1_1$Context=='Slow' & 
                   dat1_1$Type=='Target'] <- 'Slow, Target'
dat1_1$Ctxt_Type <- factor(dat1_1$Ctxt_Type,
                           levels=c('Slow, Target','Slow, Context',
                                    'Fast, Target','Fast, Context'))

dat2 <- aggregate(GradedJudge ~ Distort*Context*Subject, data=jdata2, FUN=mean)

dat2_1 <- aggregate(GradedJudge ~ Distort*Context, data=dat2, FUN=mean)
dat2_1_se <- aggregate(GradedJudge ~ Distort*Context, data=dat2, 
                       FUN=function(x) sd(x)/sqrt(length(x)))
dat2_1$UL <- dat2_1$GradedJudge + dat2_1_se$GradedJudge
dat2_1$LL <- dat2_1$GradedJudge - dat2_1_se$GradedJudge

names(dat1_1)[3] <- 'jdgmt'
names(dat2_1)[3] <- 'jdgmt'

dat2_1$Context <- factor(dat2_1$Context, levels=c('slow','fast'), 
                         labels=c('Slow','Fast'))

dat2_1$Type <- factor(ifelse(dat2_1$Distort %in% c(-12,-6,0,6,12), 'Target', 'Contextual'))

dat2_1$Ctxt_Type <- character(nrow(dat2_1))
dat2_1$Ctxt_Type[dat2_1$Context=='Fast' & 
                   dat2_1$Type=='Contextual'] <- 'Fast, Context'
dat2_1$Ctxt_Type[dat2_1$Context=='Fast' & 
                   dat2_1$Type=='Target'] <- 'Fast, Target'
dat2_1$Ctxt_Type[dat2_1$Context=='Slow' & 
                   dat2_1$Type=='Contextual'] <- 'Slow, Context'
dat2_1$Ctxt_Type[dat2_1$Context=='Slow' & 
                   dat2_1$Type=='Target'] <- 'Slow, Target'
dat2_1$Ctxt_Type <- factor(dat2_1$Ctxt_Type,
                           levels=c('Slow, Target','Slow, Context',
                                    'Fast, Target','Fast, Context'))

dat1_1$jdg_type <- 'Binary Judgments'
dat2_1$jdg_type <- 'Graded Judgments'

dat3 <- rbind(dat1_1, dat2_1)
dat3$jdg_type <- factor(dat3$jdg_type)

p1 <- ggplot(data=subset(dat3, jdg_type=='Binary Judgments'), 
             aes(x=Distort, y=jdgmt, group=Context)) +
  geom_point(aes(x=Distort, y=jdgmt, shape=Ctxt_Type), col='black', cex=3) +
  geom_errorbar(aes(x=Distort, ymin=LL, ymax=UL), width=1.5) +
  geom_line(aes(x=Distort, y=jdgmt)) +
  theme_bw(base_size = 14) +
  theme(panel.grid=element_blank(), legend.title = element_blank(),
        axis.text.x = element_text(angle=45, size=8),
        legend.position='none') +
  geom_abline(aes(slope=0, intercept=.50), lty=2) +
  geom_vline(aes(xintercept=0), lty=2) +
  scale_x_continuous(breaks=c(-30,-24,-18,-12,-6,0,6,12,18,24,30), expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0)) +
  ylab('Proportion Judged Faster') +
  xlab('Percent Tempo Distortion') +
  scale_shape_manual(values=c(1,2,16,17)) +
  ggtitle('A')

p2 <- ggplot(data=subset(dat3, jdg_type=='Graded Judgments'), 
             aes(x=Distort, y=jdgmt, group=Context)) +
  geom_point(aes(x=Distort, y=jdgmt, shape=Ctxt_Type), col='black', cex=3) +
  geom_errorbar(aes(x=Distort, ymin=LL, ymax=UL), width=1.5) +
  geom_line(aes(x=Distort, y=jdgmt)) +
  theme_bw(base_size = 14) +
  theme(panel.grid=element_blank(), legend.title = element_blank(),
        strip.background = element_blank(), axis.text.x = element_text(angle=45, size=8),
        legend.position=c(.8,.2)) +
  geom_abline(aes(slope=0, intercept=0), lty=2) +
  geom_vline(aes(xintercept=0), lty=2) +
  scale_x_continuous(breaks=c(-30,-24,-18,-12,-6,0,6,12,18,24,30), expand=c(0,0)) +
  scale_y_continuous(expand=c(0,0), limits=c(-3,3)) +
  ylab('Graded Judgment') +
  xlab('Percent Tempo Distortion') +
  scale_shape_manual(values=c(1,2,16,17)) +
  ggtitle('B')

grid.arrange(p1, p2, nrow=1, ncol=2)

