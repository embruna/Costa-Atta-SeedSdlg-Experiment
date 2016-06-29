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
# there were three seedlings of each species in each plot; this adds a column with the replicate
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





#################################################################################
# GLMM SEEDS
# Otimo Exemplo: http://www.simonqueenborough.info/R/specialist/mixed-models.html
#################################################################################



# Some summary tables
# How many Survived? 

# This is "per tray"

seeds["seeds_remaining"] <- NA
seeds$seeds_remaining<-seeds$seeds_initial-seeds$seeds_removed
seeds["seeds_prop_remaining"] <- NA
seeds$seeds_prop_remaining<-(seeds$seeds_remaining/seeds$seeds_initial)
summary(seeds)

# For analyses, do the proportion, no. remaining per nest in each treatment (pool all trays for a nest) 

seeds_removed_summary<-table(seeds$seeds_prop_remaining)


# How many Survived by Species? 
survival_summary<-table(seedlings$fate, seedlings$trt)
survival_summary
# Now the GLMMs
options(na.action = "na.fail") #for calc of QAIC see page 42: https://cran.r-project.org/web/packages/MuMIn/MuMIn.pdf
