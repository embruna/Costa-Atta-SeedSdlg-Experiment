# library(ggbiplot)
# library(grid)
# library(gridExtra)
# library(reshape2)
# library(MuMIn)
# library(arm)
# library(popbio)

library(broom)
library(lme4)
library(ggplot2)
library(dplyr)
library(tidyr)

#Clear out everything from the environment
rm(list=ls())

######################################################
######################################################
### DATA ENTRY AND CLEANUP
######################################################
######################################################

#Step 1: load the individual CSV files and save them as dataframes

seeds<-read.csv("./Data/seeds.csv", dec=".", header = TRUE, sep = ",", check.names=FALSE )
seedlings<-read.csv("./Data/seedlings.csv", dec=".", header = TRUE, sep = ",", check.names=FALSE )

# Change columns to factors as needed
str(seedlings)
# This will make it easier to remember what 0 & 1 are and for making plots later
seedlings["fate"] <- NA
seedlings = within(seedlings, {fate = ifelse(surv == 0, "alive", "dead")})
# there were three summary_1 of each species in each plot; this adds a column with the replicate
seedlings["sdlg_replicate"] <- NA 
seedlings$sdlg_replicate<-rep(seq(1,3, by=1), times=480)
seedlings<-rename(seedlings, pair=trail)

# Change columns to factors or integers as needed
seedlings$nest<-as.factor(seedlings$nest)
seedlings$pair<-as.factor(seedlings$pair)
seedlings$plot<-as.factor(seedlings$plot)
seedlings$fate<-as.factor(seedlings$fate)
seedlings$sdlg_replicate<-as.integer(seedlings$sdlg_replicate)



str(seeds)
seeds$nest<-as.factor(seeds$nest)
seeds$plate<-as.factor(seeds$plate)




#################################################################################
# GLMM SEEDLINGS
# Otimo Exemplo: http://www.simonqueenborough.info/R/specialist/mixed-models.html
#################################################################################

# Some summary tables
# How many Survived? 
survival_summary<-table(seedlings$fate)
# How many Survived by Species? 
survival_summary<-table(seedlings$fate, seedlings$trt)
survival_summary
# Now the GLMMs
options(na.action = "na.fail") #for calc of QAIC see page 42: https://cran.r-project.org/web/packages/MuMIn/MuMIn.pdf

# Simplest models: effect of trt vs trt plus random effect of nest
# sdlg.glm <- glm(surv ~ trt, data = seedlings, family = binomial) #effect of TRT
sdlg.glmer1 <- glmer(surv ~ (1|nest), family = binomial, data = seedlings) #Random effect of nest
sdlg.glmer2 <- glmer(surv ~ trt+(1|nest), family = binomial, data = seedlings) #TRT + Random effect of nest
sdlg.glmer3 <- glmer(surv ~ trt+(1|pair/nest), family = binomial, data = seedlings) #TRT + Random effect of pair nested in nest
sdlg.glmer4 <- glmer(surv ~ trt+spp+(1|nest), family = binomial, data = seedlings) #TRT + SPP Random effect of nest
sdlg.glmer5 <- glmer(surv ~ trt+spp+(1|pair/nest), family = binomial, data = seedlings) #TRT + SPP Random effect of pair nested in nest
sdlg.glmer6 <- glmer(surv ~ trt*spp+(1|nest), family = binomial, data = seedlings) #TRT + SPP + TRT* SPP + Random effect of nest
sdlg.glmer7 <- glmer(surv ~ trt*spp+(1|pair/nest), family = binomial, data = seedlings) #TRT + SPP + TRT* SPP + Random effect of oair nested in nest

summary(sdlg.glmer1)
summary(sdlg.glmer2)
summary(sdlg.glmer3)
summary(sdlg.glmer4)
summary(sdlg.glmer5)
summary(sdlg.glmer6)
summary(sdlg.glmer7)
# Is it necessary to include PLOT as a random effect?
# Test this using likelihood ratio test, as recommended by Bolker et al. 2009
LL1 <- logLik(sdlg.glmer5)
LL2 <- logLik(sdlg.glmer7)
D <- as.numeric(-2 * (LL1 - LL2))
pchisq(D, df = 1, lower = FALSE)[1]
# A: YES

# Could also do like this?
anova(sdlg.glmer1, sdlg.glmer2, test = "Chisq") #Add TRT improves fit
anova(sdlg.glmer2, sdlg.glmer3, test = "Chisq") #Adding Plot nested in Nest improves fit
anova(sdlg.glmer3, sdlg.glmer5, test = "Chisq") #Adding species ID improves fit
anova(sdlg.glmer5, sdlg.glmer7, test = "Chisq") #Adding TRT by Species Interaction doesn't improve fit (NB: worried about model convergence)
anova(sdlg.glmer4, sdlg.glmer6, test = "Chisq") #Adding Interaction also doesn't improve fit when not nesting plot in nest 

AIC(sdlg.glmer1,sdlg.glmer2,sdlg.glmer3,sdlg.glmer4,sdlg.glmer5,sdlg.glmer6,sdlg.glmer7)
#Best model includes treatment and species but not treatment by species interaction!!!

# TO REPORT: http://www.ashander.info/posts/2015/10/model-selection-glms-aic-what-to-report/
model.names.sdlgs <- c("1 Nest as a Random Effect", 
                    "2 Treatment + Random Effect of Nest", 
                    "3 Treatment + replicate nested in the random effect of Nest", 
                    "4 Treatment + Species + Random Effect of Nest",
                    "5 Treatment + Species + replicate nested in the random effect of Nest",
                    "6 Treatment + Species + Treatment * Species Interaction + Random Effect of Nest",
                    "7 Treatment + Species + Treatment * Species Interaction + replicate nested in the random effect of Nest")


summ.table.sdgls <- do.call(rbind, lapply(list(sdlg.glmer1,sdlg.glmer2,sdlg.glmer3,sdlg.glmer4,
                                             sdlg.glmer5,sdlg.glmer6,sdlg.glmer7), broom::glance))
summ.table.sdgls
table.sdlgs <- c("df.residual", "deviance", "AIC")
reported.table.sdlgs <- summ.table.sdgls[table.sdlgs]
names(reported.table.sdlgs) <- c("Resid. Df", "Resid. Dev", "AIC")

reported.table.sdlgs[['dAIC']] <-  with(reported.table.sdlgs, AIC - min(AIC))
reported.table.sdlgs[['weight']] <- with(reported.table.sdlgs, exp(- 0.5 * dAIC) / sum(exp(- 0.5 * dAIC)))
reported.table.sdlgs$AIC <- NULL
reported.table.sdlgs$weight <- round(reported.table.sdlgs$weight, 2)
reported.table.sdlgs$dAIC <- round(reported.table.sdlgs$dAIC, 1)
row.names(reported.table.sdlgs) <- model.names.sdlgs
reported.table.sdlgs

write.table(reported.table.sdlgs, file = "Sdlgs_GLMM.csv", , sep = ",", na = "NA", row.names = TRUE)

# Overall survival %
survival_summary1<-table(seedlings$fate)
survival_summary1.prop<-prop.table(as.matrix(survival_summary1))
survival_summary1.prop #Proportions in each category

# Surv % by trt
survival_summary2<-table(seedlings$fate, seedlings$trt)
survival_summary2.prop<-prop.table(as.matrix(survival_summary2))
survival_summary2.prop #Proportions in each category

# Surv % by species
survival_summary3<-table(seedlings$fate, seedlings$spp)
survival_summary3.prop<-prop.table(as.matrix(survival_summary3))
survival_summary3.prop #Proportions in each category

##########################
###### SDLG FIGURES ######
##########################
summary(seedlings)
str(seedlings)
sdlg_summary<-group_by(seedlings,trt,spp)

sdlg_summary <- summarize(sdlg_summary, dead_sdlgs=sum(surv), init_sdlgs=60)
sdlg_summary<-droplevels(sdlg_summary)
sdlg_summary<-as.data.frame(sdlg_summary)
sdlg_summary["mort.rate"] <- NA 
sdlg_summary$mort.rate<-sdlg_summary$dead_sdlgs/sdlg_summary$init_sdlgs






sdlg_summary<-spread(sdlg_summary, trt, surv_sdlgs)
sdlg_summary["prop.dead.closed"] <- NA 
sdlg_summary["prop.dead.open"] <- NA 

sdlg_summary$prop.dead.closed<-sdlg_summary$close/sdlg_summary$init_sdlgs
sdlg_summary$prop.dead.open<-sdlg_summary$open/sdlg_summary$init_sdlgs

sdlg_summary["increase_surv"] <- NA 
sdlg_summary$increase_surv<-(sdlg_summary$prop.dead.open-sdlg_summary$prop.dead.closed)

levels(sdlg_summary$spp)
levels(sdlg_summary$spp) <- c("Alibertia", "Brosimum ","Coussarea","Eugenia","Guapira", "Maprounea", "Matayba", "Miconia", "Myrcia", "Siparuna", "Tapirira", "Virola")

sdlg.fig.1 <- ggplot(sdlg_summary, aes(x=spp, y=mort.rate, fill=trt))+
  geom_bar(stat="identity", position=position_dodge(), colour="black")+
  scale_fill_grey() + 
  theme_classic()+
  scale_y_continuous(breaks=seq(0,0.8,0.1))+
  ylab("Mortaltiy Rate")+
  xlab("Species")+
  ggtitle("Mortality of seedlings (protected vs. exposed)")
sdlg.fig.1



#################################################################################
# GLMM SEEDS
# Otimo Exemplo: http://www.simonqueenborough.info/R/specialist/mixed-models.html
#################################################################################

# This is "per tray"
seeds["seeds_remaining"] <- NA
seeds$seeds_remaining<-seeds$seeds_initial-seeds$seeds_removed
seeds["seeds_prop_remaining"] <- NA
seeds$seeds_prop_remaining<-(seeds$seeds_remaining/seeds$seeds_initial)
summary(seeds)


seeds.glmer0 <- glmer(cbind(seeds$seeds_remaining, seeds$seeds_removed) ~ (1|nest), family = binomial, data = seeds) #Random effect of nest
seeds.glmer1 <- glmer(cbind(seeds$seeds_remaining, seeds$seeds_removed) ~ (1|plate/nest), family = binomial, data = seeds) #Random effect of nest
seeds.glmer2 <- glmer(cbind(seeds$seeds_remaining, seeds$seeds_removed) ~ trt+(1|nest), family = binomial, data = seeds) #TRT + Random effect of nest
seeds.glmer3 <- glmer(cbind(seeds$seeds_remaining, seeds$seeds_removed) ~ spp+(1|nest), family = binomial, data = seeds) #TRT + Random effect of nest
seeds.glmer4 <- glmer(cbind(seeds$seeds_remaining, seeds$seeds_removed) ~ trt+spp+(1|nest), family = binomial, data = seeds) #TRT + SPP Random effect of nest
seeds.glmer5 <- glmer(cbind(seeds$seeds_remaining, seeds$seeds_removed) ~ trt*spp+(1|nest), family = binomial, data = seeds) #TRT + SPP + TRT* SPP + Random effect of nest
seeds.glmer6 <- glmer(cbind(seeds$seeds_remaining, seeds$seeds_removed) ~ trt+(1|plate/nest), family = binomial, data = seeds) #TRT + Random effect of pair nested in nest
seeds.glmer7 <- glmer(cbind(seeds$seeds_remaining, seeds$seeds_removed) ~ trt+spp+(1|plate/nest), family = binomial, data = seeds) #TRT + SPP Random effect of pair nested in nest
#seeds.glmer8 <- glmer(cbind(seeds$seeds_remaining, seeds$seeds_removed) ~ trt*spp+(1|plate/nest), family = binomial, data = seeds) #TRT + SPP + TRT* SPP + Random effect of oair nested in nest
# Model 8 doesn't converge
anova(seeds.glmer0, seeds.glmer1, test = "Chisq") #No improved fit to adding plate nested in nest

summary(seeds.glmer0)
summary(seeds.glmer1)
summary(seeds.glmer3)
summary(seeds.glmer4)
summary(seeds.glmer5)

AIC(seeds.glmer1,seeds.glmer2,seeds.glmer3,seeds.glmer4,seeds.glmer5,sdlg.glmer6,sdlg.glmer7)


###NO ADVANTAGE TO NESTING PLATES IN NEST, SO POOL!!!!

#######################################################################################################
### THIS POOLS ALL THE DATA FROM DIFFERENT PLATES AT A NEST AND DOES THE ANALYSIS ON THE POOLED DATA###
#######################################################################################################
# For analyses, do the proportion, no. remaining per nest in each treatment (pool all trays for a nest) 
# this divides the dataframe into two dataframes, one for the species for which there was 1 seed per tray...
seeds_one<-filter(seeds, seeds_initial==1)
summary(seeds_one)

# ...and another for the seeds for whihc there were three seeds per tray
seeds_three<-filter(seeds, seeds_initial==3)
summary(seeds_three)

# Pool the data from the different plates and add a column with the proportion of seeds remaining
summary_3<-group_by(seeds_three,nest, trt,spp)
summary_3 <- summarize(summary_3, seeds_initial=sum(seeds_initial),seeds_removed=sum(seeds_removed),seeds_remaining=sum(seeds_remaining))
summary_3<-mutate(summary_3, prop_remaining=seeds_remaining/seeds_initial)
summary_3<-droplevels(summary_3)
summary_3<-as.data.frame(summary_3)

summary_1<-group_by(seeds_one,nest, trt,spp)
summary_1 <- summarize(summary_1, seeds_initial=sum(seeds_initial),seeds_removed=sum(seeds_removed),seeds_remaining=sum(seeds_remaining))
summary_1<-mutate(summary_1, prop_remaining=seeds_remaining/seeds_initial)
summary_1<-droplevels(summary_1)
summary_1<-as.data.frame(summary_1)

#ALL TOGETHER
pooled<-group_by(seeds,nest, trt,spp)
pooled <- summarize(pooled, seeds_initial=sum(seeds_initial),seeds_removed=sum(seeds_removed),seeds_remaining=sum(seeds_remaining))
pooled<-mutate(pooled, prop_remaining=seeds_remaining/seeds_initial)
pooled<-droplevels(pooled)
pooled<-as.data.frame(pooled)


# Now the GLMMs

############################################
###### ALL SPP TOGETHER #######
############################################

options(na.action = "na.fail") #for calc of QAIC see page 42: https://cran.r-project.org/web/packages/MuMIn/MuMIn.pdf

# Simplest models: effect of trt vs trt plus random effect of nest
# seeds.glm <- glm(surv ~ trt, data = summary_1, family = binomial) #effect of TRT
seeds.glmer1 <- glmer(cbind(pooled$seeds_remaining, pooled$seeds_removed) ~ (1|nest), family = binomial, data = pooled) #Random effect of nest
seeds.glmer2 <- glmer(cbind(pooled$seeds_remaining, pooled$seeds_removed) ~ trt+(1|nest), family = binomial, data = pooled) #TRT + Random effect of nest
seeds.glmer3 <- glmer(cbind(pooled$seeds_remaining, pooled$seeds_removed) ~ spp+(1|nest), family = binomial, data = pooled) #TRT + Random effect of nest
seeds.glmer4 <- glmer(cbind(pooled$seeds_remaining, pooled$seeds_removed) ~ trt+spp+(1|nest), family = binomial, data = pooled) #TRT + SPP Random effect of nest
seeds.glmer5 <- glmer(cbind(pooled$seeds_remaining, pooled$seeds_removed) ~ trt*spp+(1|nest), family = binomial, data = pooled) #TRT + SPP + TRT* SPP + Random effect of nest

summary(seeds.glmer1)
summary(seeds.glmer2)
summary(seeds.glmer3)
summary(seeds.glmer4)
summary(seeds.glmer5)


# Could also do like this?
anova(seeds.glmer1, seeds.glmer2, test = "Chisq") #Add TRT improves fit
anova(seeds.glmer1, seeds.glmer3, test = "Chisq") #Adding spp improves fit
anova(seeds.glmer1, seeds.glmer4, test = "Chisq") #Adding both  improves fit
anova(seeds.glmer2, seeds.glmer4, test = "Chisq") #Adding both improves fit over trt
anova(seeds.glmer3, seeds.glmer4, test = "Chisq") #Adding both improves fit over spp
anova(seeds.glmer4, seeds.glmer5, test = "Chisq") #Adding trt by Species Interaction improve fit over just main effects

AIC(seeds.glmer1,seeds.glmer2,seeds.glmer3,seeds.glmer4,seeds.glmer5)
#Best model includes treatment and species and treatment by species interaction!!!

# TO REPORT: http://www.ashander.info/posts/2015/10/model-selection-glms-aic-what-to-report/
model.names.seeds <- c("1 Nest as a Random Effect", 
                       "2 Treatment + Random Effect of Nest", 
                       "3 Species   + Random Effect of Nest",
                       "4 Treatment + Species + Random effect of Nest",
                       "5 Treatment + Species + Treatment * Species Interaction + Random effect of Nest")


summ.table.seeds <- do.call(rbind, lapply(list(seeds.glmer1,seeds.glmer2,seeds.glmer3,seeds.glmer4, seeds.glmer5), broom::glance))
summ.table.seeds
table.seeds <- c("df.residual", "deviance", "AIC")
reported.table.seeds <- summ.table.seeds[table.seeds]
names(reported.table.seeds) <- c("Resid. Df", "Resid. Dev", "AIC")

reported.table.seeds[['dAIC']] <-  with(reported.table.seeds, AIC - min(AIC))
reported.table.seeds[['weight']] <- with(reported.table.seeds, exp(- 0.5 * dAIC) / sum(exp(- 0.5 * dAIC)))
reported.table.seeds$AIC <- NULL
reported.table.seeds$weight <- round(reported.table.seeds$weight, 2)
reported.table.seeds$dAIC <- round(reported.table.seeds$dAIC, 1)
row.names(reported.table.seeds) <- model.names.seeds
reported.table.seeds

write.table(reported.table.seeds, file = "pooled.csv", , sep = ",", na = "NA", row.names = TRUE)









############################################
###### SPP WITH THREE SEEDS PER TRAY #######
############################################
options(na.action = "na.fail") #for calc of QAIC see page 42: https://cran.r-project.org/web/packages/MuMIn/MuMIn.pdf

# Simplest models: effect of trt vs trt plus random effect of nest
# seeds.glm <- glm(surv ~ trt, data = summary_1, family = binomial) #effect of TRT
seeds.glmer1 <- glmer(cbind(summary_3$seeds_remaining, summary_3$seeds_removed) ~ (1|nest), family = binomial, data = summary_3) #Random effect of nest
seeds.glmer2 <- glmer(cbind(summary_3$seeds_remaining, summary_3$seeds_removed) ~ trt+(1|nest), family = binomial, data = summary_3) #TRT + Random effect of nest
seeds.glmer3 <- glmer(cbind(summary_3$seeds_remaining, summary_3$seeds_removed) ~ spp+(1|nest), family = binomial, data = summary_3) #TRT + Random effect of nest
seeds.glmer4 <- glmer(cbind(summary_3$seeds_remaining, summary_3$seeds_removed) ~ trt+spp+(1|nest), family = binomial, data = summary_3) #TRT + SPP Random effect of nest
seeds.glmer5 <- glmer(cbind(summary_3$seeds_remaining, summary_3$seeds_removed) ~ trt*spp+(1|nest), family = binomial, data = summary_3) #TRT + SPP + TRT* SPP + Random effect of nest

summary(seeds.glmer1)
summary(seeds.glmer2)
summary(seeds.glmer3)
summary(seeds.glmer4)
summary(seeds.glmer5)


# Could also do like this?
anova(seeds.glmer1, seeds.glmer2, test = "Chisq") #Add TRT improves fit
anova(seeds.glmer1, seeds.glmer3, test = "Chisq") #Adding spp improves fit
anova(seeds.glmer1, seeds.glmer4, test = "Chisq") #Adding both  improves fit
anova(seeds.glmer2, seeds.glmer4, test = "Chisq") #Adding both improves fit over trt
anova(seeds.glmer3, seeds.glmer4, test = "Chisq") #Adding both improves fit over spp
anova(seeds.glmer4, seeds.glmer5, test = "Chisq") #Adding trt by Species Interaction improve fit over just main effects

AIC(seeds.glmer1,seeds.glmer2,seeds.glmer3,seeds.glmer4,seeds.glmer5)
#Best model includes treatment and species and treatment by species interaction!!!

# TO REPORT: http://www.ashander.info/posts/2015/10/model-selection-glms-aic-what-to-report/
model.names.seeds <- c("1 Nest as a Random Effect", 
                       "2 Treatment + Random Effect of Nest", 
                       "3 Species   + Random Effect of Nest",
                       "4 Treatment + Species + Random effect of Nest",
                       "5 Treatment + Species + Treatment * Species Interaction + Random effect of Nest")


summ.table.seeds <- do.call(rbind, lapply(list(seeds.glmer1,seeds.glmer2,seeds.glmer3,seeds.glmer4, seeds.glmer5), broom::glance))
summ.table.seeds
table.seeds <- c("df.residual", "deviance", "AIC")
reported.table.seeds <- summ.table.seeds[table.seeds]
names(reported.table.seeds) <- c("Resid. Df", "Resid. Dev", "AIC")

reported.table.seeds[['dAIC']] <-  with(reported.table.seeds, AIC - min(AIC))
reported.table.seeds[['weight']] <- with(reported.table.seeds, exp(- 0.5 * dAIC) / sum(exp(- 0.5 * dAIC)))
reported.table.seeds$AIC <- NULL
reported.table.seeds$weight <- round(reported.table.seeds$weight, 2)
reported.table.seeds$dAIC <- round(reported.table.seeds$dAIC, 1)
row.names(reported.table.seeds) <- model.names.seeds
reported.table.seeds

write.table(reported.table.seeds, file = "seeds_3_GLMM.csv", , sep = ",", na = "NA", row.names = TRUE)




############################################
###### SPP WITH ONE SEED PER TRAY #######
############################################

options(na.action = "na.fail") #for calc of QAIC see page 42: https://cran.r-project.org/web/packages/MuMIn/MuMIn.pdf

# Simplest models: effect of trt vs trt plus random effect of nest
seeds.glmer1.1 <- glmer(cbind(summary_1$seeds_remaining, summary_1$seeds_removed) ~ (1|nest), family = binomial, data = summary_1) #Random effect of nest
seeds.glmer1.2 <- glmer(cbind(summary_1$seeds_remaining, summary_1$seeds_removed) ~ trt+(1|nest), family = binomial, data = summary_1) #TRT + Random effect of nest
seeds.glmer1.3 <- glmer(cbind(summary_1$seeds_remaining, summary_1$seeds_removed) ~ spp+(1|nest), family = binomial, data = summary_1) #TRT + Random effect of nest
seeds.glmer1.4 <- glmer(cbind(summary_1$seeds_remaining, summary_1$seeds_removed) ~ trt+spp+(1|nest), family = binomial, data = summary_1) #TRT + SPP Random effect of nest
seeds.glmer1.5 <- glmer(cbind(summary_1$seeds_remaining, summary_1$seeds_removed) ~ trt*spp+(1|nest), family = binomial, data = summary_1) #TRT + SPP + TRT* SPP + Random effect of nest

summary(seeds.glmer1.1)
summary(seeds.glmer1.2)
summary(seeds.glmer1.3)
summary(seeds.glmer1.4)
summary(seeds.glmer1.5)


# Could also do like this?
anova(seeds.glmer1.1, seeds.glmer1.2, test = "Chisq") #Add TRT improves fit
anova(seeds.glmer1.1, seeds.glmer1.3, test = "Chisq") #Adding spp improves fit
anova(seeds.glmer1.1, seeds.glmer1.4, test = "Chisq") #Adding both  improves fit
anova(seeds.glmer1.2, seeds.glmer1.4, test = "Chisq") #Adding both improves fit over trt
anova(seeds.glmer1.3, seeds.glmer1.4, test = "Chisq") #Adding both improves fit over spp
anova(seeds.glmer1.4, seeds.glmer1.5, test = "Chisq") #Adding trt by Species Interaction improve fit over just main effects

AIC(seeds.glmer1.1,seeds.glmer1.2,seeds.glmer1.3,seeds.glmer1.4,seeds.glmer1.5)
#Best model includes treatment and species and treatment by species interaction!!!

# TO REPORT: http://www.ashander.info/posts/2015/10/model-selection-glms-aic-what-to-report/
model.names.seeds.1 <- c("1 Nest as a Random Effect", 
                       "2 Treatment + Random Effect of Nest", 
                       "3 Species   + Random Effect of Nest",
                       "4 Treatment + Species + Random effect of Nest",
                       "5 Treatment + Species + Treatment * Species Interaction + Random effect of Nest")


summ.table.seeds.1 <- do.call(rbind, lapply(list(seeds.glmer1.1,seeds.glmer1.2,seeds.glmer1.3,seeds.glmer1.4, seeds.glmer1.5), broom::glance))
summ.table.seeds.1
table.seeds.1 <- c("df.residual", "deviance", "AIC")
reported.table.seeds.1 <- summ.table.seeds[table.seeds.1]
names(reported.table.seeds.1) <- c("Resid. Df", "Resid. Dev", "AIC")

reported.table.seeds.1[['dAIC']] <-  with(reported.table.seeds.1, AIC - min(AIC))
reported.table.seeds.1[['weight']] <- with(reported.table.seeds.1, exp(- 0.5 * dAIC) / sum(exp(- 0.5 * dAIC)))
reported.table.seeds.1$AIC <- NULL
reported.table.seeds.1$weight <- round(reported.table.seeds.1$weight, 2)
reported.table.seeds.1$dAIC <- round(reported.table.seeds.1$dAIC, 1)
row.names(reported.table.seeds.1) <- model.names.seeds.1
reported.table.seeds.1

write.table(reported.table.seeds.1, file = "seeds_1_GLMM.csv", , sep = ",", na = "NA", row.names = TRUE)


#########################
######## FIGURES ########
#########################

#All plant species

levels(pooled$spp)
levels(pooled$spp) <- c("Alibertia", "Brosimum ","Coussarea","Eugenia","Guapira", "Maprounea", "Matayba", "Miconia", "Myrcia", "Siparuna", "Tapirira", "Virola")

seed_fig_pooled.1 <- ggplot(pooled, aes(x=trt, y=seeds_remaining))+geom_boxplot()+ theme_classic() # EFFECT OF TRT
seed_fig_pooled.2 <- ggplot(pooled, aes(x=spp, y=seeds_remaining))+geom_boxplot()+ theme_classic() # EFFECT OF SPP

seed_fig_pooled <- ggplot(pooled, aes(x=trt, y=seeds_remaining))+geom_boxplot()+ theme_classic()
seed_fig_pooled <- seed_fig_pooled+facet_wrap(~ spp, ncol=3)
seed_fig_pooled<-seed_fig_pooled+labs(y="Seeds Remaining", x="Treatment", title="Fig. 1")
seed_fig_pooled

# Seeds with 3 per plate
levels(summary_3$spp)
levels(summary_3$spp) <- c("Alibertia", "Coussarea","Guapira", "Maprounea", "Matayba", "Miconia", "Myrcia", "Siparuna", "Tapirira", "Virola")

seed_fig_1.1 <- ggplot(summary_3, aes(x=trt, y=seeds_remaining))+geom_boxplot()+ theme_classic() # EFFECT OF TRT
seed_fig_1.2 <- ggplot(summary_3, aes(x=spp, y=seeds_remaining))+geom_boxplot()+ theme_classic() # EFFECT OF SPP

seed_fig_1 <- ggplot(summary_3, aes(x=trt, y=seeds_remaining))+geom_boxplot()+ theme_classic()
seed_fig_1 <- seed_fig_1+facet_wrap(~ spp, ncol=2)
seed_fig_1<-seed_fig_1+labs(y="Seeds Remaining", x="Treatment", title="Fig. 1")
seed_fig_1

# Seeds with 1 per plate
levels(summary_1$spp)
levels(summary_1$spp) <- c("Brosimum ", "Eugenia")

seed_fig_2.1 <- ggplot(summary_1, aes(x=trt, y=seeds_remaining))+geom_boxplot()+ theme_classic() # EFFECT OF TRT
seed_fig_2.2 <- ggplot(summary_1, aes(x=spp, y=seeds_remaining))+geom_boxplot()+ theme_classic() # EFFECT OF SPP

seed_fig_2 <- ggplot(summary_1, aes(x=trt, y=seeds_remaining))+geom_boxplot()+ theme_classic() #TRT x SPP
seed_fig_2 <- seed_fig_2+facet_wrap(~ spp, ncol=2)
seed_fig_2<-seed_fig_2+labs(y="Seeds Remaining", x="Treatment", title="Fig. 2")
seed_fig_2






