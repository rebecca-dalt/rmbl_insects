## daltonetal-malaise-analysis-2022.R
## analysis for Dalton, RM., Underwood, NC., Inouye, DW., Soule M., and B. Inouye. ## long-term declines in insect abundance and biomass in a subalpine habitat. 
## last updated: 13 July 2023

## for more information about this data set, 
## see "daltonetal-malaise-metadata-2022.pdf" 

## Files needed for this analysis:
  # "daltonetal-malaise-week-2020-final.csv": weekly data from 1984 - 2020
  # "daltonetal-malaise-year-2020-final.csv": yearly data from 1984 - 2020

## additional information about data used in this analysis: 
  # snowmelt date and total snowfall data collected by billy barr at RMBL:
  # http://www.gothicwx.org 
  # floral abundance data were collected by David Inouye and colleagues at RMBL:
  # https://www.bio.fsu.edu/~nunderwood/homepage/RMBLphenologyproject.html 
  # temperature and summer precipitation data from Crested Butte Weather Station
  # Station USC00051959 
  # https://www.ncdc.noaa.gov/cdo-web/datasets/GHCND/stations/GHCND:USC00051959/detail

rm(list=ls())
wd <- ("")
#add your working directory here
setwd(wd)

# packages for analysis
library(car) 		     # for vif() function
library(lmerTest) 	 # for lmer and associated bits
library(merTools) 	 # for prediction intervals from lmer model fits
library(tidyverse)	 # for summarizing
library(lmerTest)
library(MuMIn)       #R.squaredGLMM

##################################################################################
###################### Open data and scale variables #############################
##################################################################################

# first, read in the weekly data frame
data.week <- read.csv("daltonetal-malaise-week-2020-final.csv", header = T, sep = ",")

# next create a reduced data set without 1984 and 1985
data.week.sub <- subset(data.week, year != 1984 & year != 1985)

# scale the variables

# full dataset (1984 - 2020):
data.week$s.year <- c(scale(data.week$year, scale = T, center = T))
data.week$s.dps <- c(scale(data.week$dps, scale = T, center = T))
data.week$s.mean.temp.ave.cb <- c(scale(data.week$mean.temp.ave.cb, scale = T, center = T))
data.week$s.sum.precip.cb <- c(scale(data.week$sum.precip.cb, scale = T, center = T))
data.week$s.mean.flw.day <- c(scale(data.week$mean.flw.day, scale = T, center = T))
data.week$s.snowpack <- c(scale(data.week$snowpack, scale = T, center = T))
# add a column for year as a random effect
data.week$f.year <- as.factor(data.week$year)

# reduced dataset (1986 - 2020):
data.week.sub$s.year <- c(scale(data.week.sub$year, scale = T, center = T))
data.week.sub$s.dps <- c(scale(data.week.sub$dps, scale = T, center = T))
data.week.sub$s.mean.temp.ave.cb <- c(scale(data.week.sub$mean.temp.ave.cb, scale = T, center = T))
data.week.sub$s.sum.precip.cb <- c(scale(data.week.sub$sum.precip.cb, scale = T, center = T))
data.week.sub$s.mean.flw.day <- c(scale(data.week.sub$mean.flw.day, scale = T, center = T))
data.week.sub$s.snowpack <- c(scale(data.week.sub$snowpack, scale = T, center = T))
# add a column for year as a random effect
data.week.sub$f.year <- as.factor(data.week.sub$year)

# second, read in the yearly data frame
data.year <- read.csv("daltonetal-malaise-year-2020-final.csv", header = T, sep = ",")

# next create a reduced data set without 1984 and 1985
data.year.sub <- subset(data.year, year != 1984 & year != 1985) 

# scale the variables

# full dataset (1984 - 2020):
data.year$s.year <- c(scale(data.year$year, center = T, scale = T))
data.year$s.summer.precip.cb <- c(scale(data.year$summer.precip.cb, center = T, scale = T))
data.year$s.summer.temp.cb <- c(scale(data.year$summer.temp.cb, center = T, scale = T))
data.year$s.snowmelt <- c(scale(data.year$snowmelt, center = T, scale = T))
data.year$s.lag.summer.precip.cb <- c(scale(data.year$lag.summer.precip.cb, center = T, scale = T))
data.year$s.lag.summer.temp.cb <- c(scale(data.year$lag.summer.temp.cb, center = T, scale = T))
data.year$s.snowpack <- c(scale(data.year$snowpack, center = T, scale = T))
data.year$s.lag.snowpack <- c(scale(data.year$lag.snowpack, center = T, scale = T))
data.year$s.yr.floralcount <- c(scale(data.year$yr.floralcount, center = T, scale = T))
data.year$s.lag.yr.floralcount <- c(scale(data.year$lag.yr.floralcount, center = T, scale = T))
data.year$f.year <- as.factor(data.year$year)

# reduced dataset (1986 - 2020):
data.year.sub$s.year <- c(scale(data.year.sub$year, center = T, scale = T))
data.year.sub$s.summer.precip.cb <- c(scale(data.year.sub$summer.precip.cb, center = T, scale = T))
data.year.sub$s.summer.temp.cb <- c(scale(data.year.sub$summer.temp.cb, center = T, scale = T))
data.year.sub$s.snowmelt <- c(scale(data.year.sub$snowmelt, center = T, scale = T))
data.year.sub$s.lag.summer.precip.cb <- c(scale(data.year.sub$lag.summer.precip.cb, center = T, scale = T))
data.year.sub$s.lag.summer.temp.cb <- c(scale(data.year.sub$lag.summer.temp.cb, center = T, scale = T))
data.year.sub$s.snowpack <- c(scale(data.year.sub$snowpack, center = T, scale = T))
data.year.sub$s.lag.snowpack <- c(scale(data.year.sub$lag.snowpack, center = T, scale = T))
data.year.sub$s.yr.floralcount <- c(scale(data.year.sub$yr.floralcount, center = T, scale = T))
data.year.sub$s.lag.yr.floralcount <- c(scale(data.year.sub$lag.yr.floralcount, center = T, scale = T))
data.year.sub$f.year <- as.factor(data.year.sub$year)

##################################################################################
############################ Statistics for manuscript ###########################
##################################################################################

# 1983 is included in these data to allow for lagged environmental variables
# subset to the sampling period

sampling.period.yr <- subset(data.year, year >= 1984) 
sampling.period.wk <- subset(data.week, year >= 1984) 

# report mean annual snowfall in the methods section
mean(sampling.period.yr$snowpack) #1045.459 mean total snowfall (cm)
(sd(sampling.period.yr$snowpack))/sqrt(length(sampling.period.yr$snowpack)) #44.29765 cm standard error total snowfall
sd(sampling.period.yr$snowpack) #269 cm standard deviation total snowfall

# report mean summer temperature in the methods section
mean(sampling.period.yr$summer.temp.cb) # 11.27071 mean summer temp (C)
(sd(sampling.period.yr$summer.temp.cb))/sqrt(length(sampling.period.yr$summer.temp.cb)) #0.1135633 cm standard error summer temp
sd(sampling.period.yr$summer.temp.cb) #0.691 C standard deviation summer temp

# report mean summer precipitation in the methods section
mean(sampling.period.yr$summer.precip.cb) # 16.61649 mean total summer precip (cm)
(sd(sampling.period.yr$summer.precip.cb))/sqrt(length(sampling.period.yr$summer.precip.cb)) #0.8947016 cm standard error total snowfall
sd(sampling.period.yr$summer.precip.cb) #5.44 cm standard deviation summer precip

# report average number of days malaise trap deployed
mean(sampling.period.wk$total.days, na.rm = T) # 2.064039 mean total days
(sd(sampling.period.wk$total.days))/sqrt(length(sampling.period.wk$total.days)) #0.02008541 days standard error
sd(sampling.period.wk$total.days) #0.41 days SD

# report average number of weekly samples
mean(sampling.period.yr$yr.samples, na.rm = T) # 10.97297 mean total weekly samples
(sd(sampling.period.yr$yr.samples))/sqrt(length(sampling.period.yr$yr.samples)) #0.584441 standard error
sd(sampling.period.yr$yr.samples) #3.6 sd days

# report total samples from 1984 to 2020
nrow(data.week) # 406 samples 
length(na.omit(data.week$total.g)) # 379 biomass samples for full data set
length(na.omit(data.week$total.no)) # 373 abundance samples for full data set
length(na.omit(data.week.sub$total.g)) # 363 biomass samples for reduced data set
length(na.omit(data.week.sub$total.no)) # 357 abundance samples for reduced data set

# report floral abundance and range with and without 2015
subset(sampling.period.yr, year == 2015)$yr.floralcount # 261272 flowers
min(subset(sampling.period.yr, year != 2015)$yr.floralcount, na.rm = T) # 21509 min total flowers
max(subset(sampling.period.yr, year != 2015)$yr.floralcount, na.rm = T) # 119590 max total flowers

##################################################################################
########## Analysis 1: Predict grams/day or number/day from weekly data ##########
##########  Results are presented in Figure 1 and Tables S1 and S2 ###############
##################################################################################

vars <-c("year", "dps", "mean.temp.ave.cb",
         "sum.precip.cb", "mean.flw.day")
variables <- data.week[vars]
cormat <- cor(variables, use = "complete.obs") #1990 has missing flower data
# all cor values below threshold of 0.7

###################################################################################
# Q1A. Total insect biomass (total.g) for both the full dataset (1984 - 2020)
# and the reduced dataset (1986 - 2020)
###################################################################################

# weekly total biomass: reduced dataset (1986 - 2020): lm
q1.86.total.g.lm <- lm(sqrt(total.g/total.days) ~  
                      s.year + s.dps + I(s.dps^2) + s.dps * s.year + 
                      s.mean.temp.ave.cb + s.mean.temp.ave.cb * s.year + 
                      s.sum.precip.cb + s.sum.precip.cb*s.year +
                      s.mean.flw.day + s.mean.flw.day*s.year, 
                      data = data.week.sub, na.action = na.omit)

summary(q1.86.total.g.lm)
vif(q1.86.total.g.lm)
plot(q1.86.total.g.lm)
r.squaredGLMM(q1.86.total.g.lm) # 0.2387614

# Weekly total biomass: reduced dataset (1986 - 2020): lmer
q1.86.total.g.lmer <- lmer(sqrt(total.g/total.days) ~
                      s.year + s.dps + I(s.dps^2) + s.dps*s.year + 
                      s.mean.temp.ave.cb + s.mean.temp.ave.cb*s.year + 
                      s.sum.precip.cb + s.sum.precip.cb*s.year +
                      s.mean.flw.day + s.mean.flw.day*s.year + (1|f.year), 
                      data = data.week.sub, na.action = na.omit,
                      REML = F) #REML = F to compare models

summary(q1.86.total.g.lmer)
vif(q1.86.total.g.lmer)
plot(q1.86.total.g.lmer)
r.squaredGLMM(q1.86.total.g.lmer) #  M = 0.2451431, C = 0.3724755

# compare fixed effects to mixed effects model
AIC(q1.86.total.g.lm, q1.86.total.g.lmer) #AIC is lower for lmer

# change REML = T
q1.86.total.g.lmer <- lmer(sqrt(total.g/total.days) ~
                      s.year + s.dps + I(s.dps^2) + s.dps*s.year + 
                      s.mean.temp.ave.cb + s.mean.temp.ave.cb*s.year + 
                      s.sum.precip.cb + s.sum.precip.cb*s.year +
                      s.mean.flw.day + s.mean.flw.day*s.year + (1|f.year), 
                      data = data.week.sub, na.action = na.omit,
                      REML = T) 

summary(q1.86.total.g.lmer) # reported in TABLE S1
vif(q1.86.total.g.lmer)
plot(q1.86.total.g.lmer)
r.squaredGLMM(q1.86.total.g.lmer) # reported in TABLE S1
CI.g.86 <-confint(q1.86.total.g.lmer) # marginal CIs for slopes
AICc(q1.86.total.g.lmer) # reported in TABLE S1: -85.07313

# next, calculate the predicted values for FIGURE 1

new.week.86 <- as.data.frame(matrix(0, ncol = 4 , nrow = length(data.week.sub$s.year)))
new.week.86 <- cbind(data.week.sub$s.year, new.week.86, data.week.sub$f.year)
names(new.week.86) <- c("s.year", "s.dps", "s.mean.temp.ave.cb", 
                        "s.sum.precip.cb", "s.mean.flw.day", "f.year")  

pred.g.86 <- as.data.frame(predictInterval(q1.86.total.g.lmer, new.week.86, 
                        which = "fixed", level = .90 )) #90% prediction interval
pred.g.86 <- cbind(data.week.sub$s.year, pred.g.86)
names(pred.g.86) <- c("s.year", "fit", "upr", "lwr") 


### backtransform, because fit was done on sqrt values
pred.g.86$fit <- (pred.g.86$fit)^2
pred.g.86$upr <- (pred.g.86$upr)^2
pred.g.86$lwr <- (pred.g.86$lwr)^2

# calculate percentage of decline
pred.g.1986 <- mean(subset(pred.g.86[,2], pred.g.86$s.year == min(pred.g.86$s.year)) )
pred.g.2020 <- mean(subset(pred.g.86[,2], pred.g.86$s.year == max(pred.g.86$s.year)) )
g.percent <- ((pred.g.2020 - pred.g.1986)/ (pred.g.1986)) * 100 
g.percent 
# note: each time code is run it produces a slight different value based on predictions, 
# ranging from ~ 46.9 to 49.0. 
# ~47 % decline reported in paper as a conservative estimate

# weekly total biomass: full dataset (1984 - 2020)

q1.total.g.lm <- lm(sqrt(total.g/total.days) ~  
                      s.year + s.dps + I(s.dps^2) + s.dps*s.year + 
                      s.mean.temp.ave.cb + s.mean.temp.ave.cb*s.year + 
                      s.sum.precip.cb + s.sum.precip.cb*s.year +
                      s.mean.flw.day + s.mean.flw.day*s.year, 
                      data = data.week, na.action = na.omit)

summary(q1.total.g.lm)
vif(q1.total.g.lm)
plot(q1.total.g.lm)
# note the fitted vs residuals plot shows pattern and truncation

q1.total.g.lmer <- lmer(sqrt(total.g/total.days) ~
                      s.year + s.dps + I(s.dps^2) + s.dps*s.year + 
                      s.mean.temp.ave.cb + s.mean.temp.ave.cb*s.year + 
                      s.sum.precip.cb + s.sum.precip.cb*s.year +
                      s.mean.flw.day + s.mean.flw.day*s.year + (1|f.year), 
                      data = data.week, na.action = na.omit,
                      REML = F)  #REML = F to compare models

plot(q1.total.g.lmer)	# residuals vs fitted looks better than lm
summary(q1.total.g.lmer)
anova(q1.total.g.lmer) 	
vif(q1.total.g.lmer)	# all about 2 or less, keep in model
AIC(q1.total.g.lm, q1.total.g.lmer) # model with random effect of year is much better fit

q1.total.g.lmer <- lmer(sqrt(total.g/total.days) ~
                      s.year + s.dps + I(s.dps^2) + s.dps*s.year + 
                      s.mean.temp.ave.cb + s.mean.temp.ave.cb*s.year + 
                      s.sum.precip.cb + s.sum.precip.cb*s.year +
                      s.mean.flw.day + s.mean.flw.day*s.year + (1|f.year), 
                      data = data.week, na.action = na.omit,
                      REML = T) 

plot(q1.total.g.lmer)
summary(q1.total.g.lmer) # reported in TABLE S2
anova(q1.total.g.lmer) 	
r.squaredGLMM(q1.total.g.lmer) # reported in TABLE S2
CI.g <-confint(q1.total.g.lmer) # marginal CIs for slopes
AICc(q1.total.g.lmer) # reported in TABLE S2  -79.1007

###################################################################################
# Q1B. Total insect abundance (total.no) for both the full dataset (1984 - 2020)
# and the reduced dataset (1986 - 2020)
###################################################################################

# Weekly total abundance: reduced dataset (1986 - 2020)  # lmer and with sqrt

q1.86.total.no.lmer <-  lmer(sqrt(total.no/total.days) ~  
                        s.year + s.dps + I(s.dps^2) + s.dps * s.year + 
                        s.mean.temp.ave.cb + s.mean.temp.ave.cb * s.year + 
                        s.sum.precip.cb + s.sum.precip.cb*s.year +
                        s.mean.flw.day + s.mean.flw.day*s.year + (1|f.year), 
                        data = data.week.sub, na.action = na.omit)

vif(q1.86.total.no.lmer)
summary(q1.86.total.no.lmer) # reported in TABLE S1
r.squaredGLMM(q1.86.total.no.lmer) # reported in TABLE S1
CI.no.86 <-confint(q1.86.total.no.lmer) # marginal CIs for slopes
AICc(q1.86.total.no.lmer) # reported in TABLE S1 - 1747.7

# prediction intervals for creating a figure

pred.no.86 <- as.data.frame(predictInterval(q1.86.total.no.lmer,
              new.week.86, which = "fixed", level = .90 )) #90% prediction interval
pred.no.86 <- cbind(data.week.sub$s.year, pred.no.86)
names(pred.no.86) <- c("s.year", "fit", "upr", "lwr") 

# backtransform, because fit was done on sqrt values

pred.no.86$fit <- (pred.no.86$fit)^2
pred.no.86$upr <- (pred.no.86$upr)^2
pred.no.86$lwr <- (pred.no.86$lwr)^2

# calculate percentage of decline

pred.no.1986 <- mean(subset(pred.no.86[,2], pred.no.86$s.year == min(pred.no.86$s.year)) )
pred.no.2020 <- mean(subset(pred.no.86[,2], pred.no.86$s.year == max(pred.no.86$s.year)) )
no.percent <- ((pred.no.2020 - pred.no.1986)/ (pred.no.1986)) * 100 
no.percent  
# note: each time code is run it produces a slight different value based on predictions, 
# ranging from ~ 61.5 to ~61.9% 
# ~61.5 % decline reported in paper as a conservative estimate

# Weekly total biomass: full data set (1984 - 2020) 

q1.total.no.lm <-  lm(sqrt(total.no/total.days) ~ 
                      s.year + s.dps + I(s.dps^2) + s.dps * s.year + 
                      s.mean.temp.ave.cb + s.mean.temp.ave.cb * s.year + 
                      s.sum.precip.cb + s.sum.precip.cb*s.year +
                      s.mean.flw.day + s.mean.flw.day*s.year , 
                      data = data.week, na.action = na.omit)

q1.total.no.lmer <-  lmer(sqrt(total.no/total.days) ~ 
                      s.year + s.dps + I(s.dps^2) + s.dps * s.year + 
                      s.mean.temp.ave.cb + s.mean.temp.ave.cb * s.year + 
                      s.sum.precip.cb + s.sum.precip.cb*s.year +
                      s.mean.flw.day + s.mean.flw.day*s.year + (1|f.year), 
                      data = data.week, na.action = na.omit, REML = F)

AIC(q1.total.no.lmer, q1.total.no.lm)  # AIC suggests lmer is better model

#### Set REML back to the default T, after AIC comparison is made

q1.total.no.lmer <-  lmer(sqrt(total.no/total.days) ~ 
                      s.year + s.dps + I(s.dps^2) + s.dps * s.year + 
                      s.mean.temp.ave.cb + s.mean.temp.ave.cb * s.year + 
                      s.sum.precip.cb + s.sum.precip.cb*s.year +
                      s.mean.flw.day + s.mean.flw.day*s.year + (1|f.year), 
                      data = data.week, na.action = na.omit, REML = T)

summary(q1.total.no.lmer) # reported in TABLE S2
vif(q1.total.no.lmer)     # all VIFs below 3
r.squaredGLMM(q1.total.no.lmer) # reported in TABLE S2
AICc(q1.total.no.lm) #[1] 1957.13

###################################################################################
# Q1 Supplement. Dipteran, Hymenopteran, and other insect abundance (total.no) and 
# biomass (total.g) for both the full dataset (1984 - 2020)
# and the reduced dataset (1986 - 2020)
# Results presented in Fig. S2, Tables S4 - S8
###################################################################################

new.week <- as.data.frame(matrix(0, ncol = 4 , nrow = length(data.week$s.year)))
new.week <- cbind(data.week$s.year, new.week, data.week$f.year)
names(new.week) <- c("s.year", "s.dps", "s.mean.temp.ave.cb", 
                     "s.sum.precip.cb", "s.mean.flw.day", "f.year")  

# Reduced dataset: 1986 - 2020
# All insects (reported in Table S1)
# Biomass (g.total)

summary(s1.86.g.total <- lmer(sqrt(total.g/total.days) ~  
                            s.year + s.dps + I(s.dps^2) + s.dps * s.year + 
                            s.mean.temp.ave.cb + s.mean.temp.ave.cb * s.year + 
                            s.sum.precip.cb + s.sum.precip.cb*s.year +
                            s.mean.flw.day + s.mean.flw.day*s.year + (1|f.year), 
                            data = data.week.sub, na.action = na.omit))
AICc(s1.86.g.total)

s1.86.pred.g.total <- predictInterval(s1.86.g.total, new.week, which = "fixed", level = .9)
s1.86.pred.g.total <- cbind(data.week$s.year, s1.86.pred.g.total) 
names(s1.86.pred.g.total) <- c("s.year", "fit", "upr", "lwr")

s1.86.pred.g.total$fit <- (s1.86.pred.g.total$fit)^2
s1.86.pred.g.total$upr <- (s1.86.pred.g.total$upr)^2
s1.86.pred.g.total$lwr <- (s1.86.pred.g.total$lwr)^2
s1.86.pred.g.total$taxa <- "all"

# Reduced dataset: 1986 - 2020
# Hymenoptera (reported in Table S3)
# Biomass (hymenoptera.g)

summary(s1.86.g.hym <- lmer(sqrt(hymenoptera.g/total.days) ~ 
                            s.year + s.dps + I(s.dps^2) + s.dps * s.year + 
                            s.mean.temp.ave.cb + s.mean.temp.ave.cb * s.year + 
                            s.sum.precip.cb + s.sum.precip.cb*s.year +
                            s.mean.flw.day + s.mean.flw.day*s.year + (1|f.year), 
                            data = data.week.sub, na.action = na.omit))
r.squaredGLMM(s1.86.g.hym)
AICc(s1.86.g.hym)

s1.86.pred.g.hym <- predictInterval(s1.86.g.hym, new.week, which = "fixed", level = .9)
s1.86.pred.g.hym <- cbind(data.week$s.year, s1.86.pred.g.hym) 
names(s1.86.pred.g.hym) <- c("s.year", "fit", "upr", "lwr")

s1.86.pred.g.hym$fit <- (s1.86.pred.g.hym$fit)^2
s1.86.pred.g.hym$upr <- (s1.86.pred.g.hym$upr)^2
s1.86.pred.g.hym$lwr <- (s1.86.pred.g.hym$lwr)^2
s1.86.pred.g.hym$taxa <- "hym"

# Reduced dataset: 1986 - 2020
# Diptera (reported in Table S5)
# Biomass (diptera.g)

summary(s1.86.g.dip <- lmer(sqrt(diptera.g/total.days) ~  
                            s.year + s.dps + I(s.dps^2) + s.dps * s.year + 
                            s.mean.temp.ave.cb + s.mean.temp.ave.cb * s.year + 
                            s.sum.precip.cb + s.sum.precip.cb*s.year +
                            s.mean.flw.day + s.mean.flw.day*s.year + (1|f.year), 
                            data = data.week.sub, na.action = na.omit))

r.squaredGLMM(s1.86.g.dip)
AICc(s1.86.g.dip)
s1.86.pred.g.dip <- predictInterval(s1.86.g.dip, new.week, which = "fixed", level = .9)
s1.86.pred.g.dip <- cbind(data.week$s.year, s1.86.pred.g.dip) 
names(s1.86.pred.g.dip) <- c("s.year", "fit", "upr", "lwr")

s1.86.pred.g.dip$fit <- (s1.86.pred.g.dip$fit)^2
s1.86.pred.g.dip$upr <- (s1.86.pred.g.dip$upr)^2
s1.86.pred.g.dip$lwr <- (s1.86.pred.g.dip$lwr)^2
s1.86.pred.g.dip$taxa <- "dip"

# Reduced dataset: 1986 - 2020
# Other (reported in Table S7)
# Biomass (other.g)

summary(s1.86.g.other <- lmer(sqrt(other.g/total.days) ~ 
                            s.year + s.dps + I(s.dps^2) + s.dps * s.year + 
                            s.mean.temp.ave.cb + s.mean.temp.ave.cb * s.year + 
                            s.sum.precip.cb + s.sum.precip.cb*s.year +
                            s.mean.flw.day + s.mean.flw.day*s.year + (1|f.year), 
                            data = data.week.sub, na.action = na.omit))
r.squaredGLMM(s1.86.g.other)
AICc(s1.86.g.other)

s1.86.pred.g.other <- predictInterval(s1.86.g.other, new.week, which = "fixed", level = .9)
s1.86.pred.g.other <- cbind(data.week$s.year, s1.86.pred.g.other) 
names(s1.86.pred.g.other) <- c("s.year", "fit", "upr", "lwr")

s1.86.pred.g.other$fit <- (s1.86.pred.g.other$fit)^2
s1.86.pred.g.other$upr <- (s1.86.pred.g.other$upr)^2
s1.86.pred.g.other$lwr <- (s1.86.pred.g.other$lwr)^2
s1.86.pred.g.other$taxa <- "other"

# Full dataset: 1984 - 2020
# Hymenoptera (reported in Table S4)
# Biomass (hymenoptera.g)

summary(s1.g.hym <- lmer(sqrt(hymenoptera.g/total.days) ~  
                        s.year + s.dps + I(s.dps^2) + s.dps * s.year + 
                        s.mean.temp.ave.cb + s.mean.temp.ave.cb * s.year + 
                        s.sum.precip.cb + s.sum.precip.cb*s.year +
                        s.mean.flw.day + s.mean.flw.day*s.year + (1|f.year), 
                        data = data.week, na.action = na.omit))
r.squaredGLMM(s1.g.hym)
AICc(s1.g.hym)

s1.pred.g.hym <- predictInterval(s1.g.hym, new.week, which = "fixed", level = .9)
s1.pred.g.hym <- cbind(data.week$s.year, s1.pred.g.hym) 
names(s1.pred.g.hym) <- c("s.year", "fit", "upr", "lwr")

s1.pred.g.hym$fit <- (s1.pred.g.hym$fit)^2
s1.pred.g.hym$upr <- (s1.pred.g.hym$upr)^2
s1.pred.g.hym$lwr <- (s1.pred.g.hym$lwr)^2
s1.pred.g.hym$taxa <- "hym"

# Full dataset: 1984 - 2020
# Diptera (reported in Table S6)
# Biomass (diptera.g)

summary(s1.g.dip <- lmer(sqrt(diptera.g/total.days) ~ 
                          s.year + s.dps + I(s.dps^2) + s.dps * s.year + 
                          s.mean.temp.ave.cb + s.mean.temp.ave.cb * s.year + 
                          s.sum.precip.cb + s.sum.precip.cb*s.year +
                          s.mean.flw.day + s.mean.flw.day*s.year  + (1|f.year), 
                          data = data.week, na.action = na.omit))
r.squaredGLMM(s1.g.dip)
AICc(s1.g.dip)

s1.pred.g.dip <- predictInterval(s1.g.dip, new.week, which = "fixed", level = .9)
s1.pred.g.dip <- cbind(data.week$s.year, s1.pred.g.dip) 
names(s1.pred.g.dip) <- c("s.year", "fit", "upr", "lwr")

s1.pred.g.dip$fit <- (s1.pred.g.dip$fit)^2
s1.pred.g.dip$upr <- (s1.pred.g.dip$upr)^2
s1.pred.g.dip$lwr <- (s1.pred.g.dip$lwr)^2
s1.pred.g.dip$taxa <- "dip"

# Full dataset: 1984 - 2020
# Other (reported in Table S8)
# Biomass (other.g)

summary(s1.g.other <- lmer(sqrt(other.g/total.days) ~ 
                          s.year + s.dps + I(s.dps^2) + s.dps * s.year + 
                          s.mean.temp.ave.cb + s.mean.temp.ave.cb * s.year + 
                          s.sum.precip.cb + s.sum.precip.cb*s.year +
                          s.mean.flw.day + s.mean.flw.day*s.year + (1|f.year), 
                          data = data.week, na.action = na.omit))
r.squaredGLMM(s1.g.other)
AICc(s1.g.other)

s1.pred.g.other <- predictInterval(s1.g.other, new.week, which = "fixed", level = .9)
s1.pred.g.other <- cbind(data.week$s.year, s1.pred.g.other) 
names(s1.pred.g.other) <- c("s.year", "fit", "upr", "lwr")

s1.pred.g.other$fit <- (s1.pred.g.other$fit)^2
s1.pred.g.other$upr <- (s1.pred.g.other$upr)^2
s1.pred.g.other$lwr <- (s1.pred.g.other$lwr)^2
s1.pred.g.other$taxa <- "other"

s1.g.models <- rbind(s1.pred.g.other, 
                     s1.pred.g.hym, 
                     s1.pred.g.dip)

s1.g.models$type <- "full"


# Combine the reduced data frame biomass prediction data frames together

s1.86.g.models <- rbind(s1.86.pred.g.other, 
                        s1.86.pred.g.hym, 
                        s1.86.pred.g.dip)

s1.86.g.models$type <- "reduced"
s1.86.g.models$s.year <- ifelse( s1.86.g.models$s.year >= -1.63912620, as.numeric(s1.86.g.models$s.year), NA)

s1.models.g <- rbind(s1.86.g.models, s1.g.models)

# Reduced dataset: 1986 - 2020
# Hymenoptera (reported in Table S3)
# Abundance (hymenoptera.no)

summary(s1.86.no.hym <- lmer(sqrt(hymenoptera.no/total.days) ~ 
                            s.year + s.dps + I(s.dps^2) + s.dps * s.year + 
                            s.mean.temp.ave.cb + s.mean.temp.ave.cb * s.year + 
                            s.sum.precip.cb + s.sum.precip.cb*s.year +
                            s.mean.flw.day + s.mean.flw.day*s.year + (1|f.year), 
                            data = data.week.sub, na.action = na.omit))
r.squaredGLMM(s1.86.no.hym)
AICc(s1.86.no.hym)

s1.86.pred.no.hym <- predictInterval(s1.86.no.hym, new.week, which = "fixed", level = .9)
s1.86.pred.no.hym <- cbind(data.week$s.year, s1.86.pred.no.hym) 
names(s1.86.pred.no.hym) <- c("s.year", "fit", "upr", "lwr")

s1.86.pred.no.hym$fit <- (s1.86.pred.no.hym$fit)^2
s1.86.pred.no.hym$upr <- (s1.86.pred.no.hym$upr)^2
s1.86.pred.no.hym$lwr <- (s1.86.pred.no.hym$lwr)^2
s1.86.pred.no.hym$taxa <- "hym"

# Reduced dataset: 1986 - 2020
# Diptera (reported in Table S5)
# Abundance (diptera.no)

summary(s1.86.no.dip <- lmer(sqrt(diptera.no/total.days) ~  
                            s.year + s.dps + I(s.dps^2) + s.dps * s.year + 
                            s.mean.temp.ave.cb + s.mean.temp.ave.cb * s.year + 
                            s.sum.precip.cb + s.sum.precip.cb*s.year +
                            s.mean.flw.day + s.mean.flw.day*s.year + (1|f.year), 
                            data = data.week.sub, na.action = na.omit))

r.squaredGLMM(s1.86.no.dip)
AICc(s1.86.no.dip)
s1.86.pred.no.dip <- predictInterval(s1.86.no.dip, new.week, which = "fixed", level = .9)
s1.86.pred.no.dip <- cbind(data.week$s.year, s1.86.pred.no.dip) 
names(s1.86.pred.no.dip) <- c("s.year", "fit", "upr", "lwr")

s1.86.pred.no.dip$fit <- (s1.86.pred.no.dip$fit)^2
s1.86.pred.no.dip$upr <- (s1.86.pred.no.dip$upr)^2
s1.86.pred.no.dip$lwr <- (s1.86.pred.no.dip$lwr)^2
s1.86.pred.no.dip$taxa <- "dip"

# Reduced dataset: 1986 - 2020
# Other (reported in Table S7)
# Abundance (other.no)

summary(s1.86.no.other <- lmer(sqrt(other.no/total.days) ~ 
                              s.year + s.dps + I(s.dps^2) + s.dps * s.year + 
                              s.mean.temp.ave.cb + s.mean.temp.ave.cb * s.year + 
                              s.sum.precip.cb + s.sum.precip.cb*s.year +
                              s.mean.flw.day + s.mean.flw.day*s.year + (1|f.year),
                              data = data.week.sub, na.action = na.omit))
r.squaredGLMM(s1.86.no.other)
AICc(s1.86.no.other)

s1.86.pred.no.other <- predictInterval(s1.86.no.other, new.week, which = "fixed", level = .9)
s1.86.pred.no.other <- cbind(data.week$s.year, s1.86.pred.no.other) 
names(s1.86.pred.no.other) <- c("s.year", "fit", "upr", "lwr")

s1.86.pred.no.other$fit <- (s1.86.pred.no.other$fit)^2
s1.86.pred.no.other$upr <- (s1.86.pred.no.other$upr)^2
s1.86.pred.no.other$lwr <- (s1.86.pred.no.other$lwr)^2
s1.86.pred.no.other$taxa <- "other"

s1.86.no.models <- rbind(s1.86.pred.no.other, 
                         s1.86.pred.no.hym, 
                         s1.86.pred.no.dip)

s1.86.no.models$type <- "reduced"
s1.86.no.models$s.year <- ifelse( s1.86.no.models$s.year >= -1.63912620, as.numeric(s1.86.no.models$s.year), NA)

# Full dataset: 1984 - 2020
# Hymenoptera (reported in Table S4)
# Abundance (hymenoptera.no)

summary(s1.no.hym <- lmer(sqrt(hymenoptera.no/total.days) ~  
                          s.year + s.dps + I(s.dps^2) + s.dps * s.year + 
                          s.mean.temp.ave.cb + s.mean.temp.ave.cb * s.year + 
                          s.sum.precip.cb + s.sum.precip.cb*s.year +
                          s.mean.flw.day + s.mean.flw.day*s.year + (1|f.year), 
                          data = data.week, na.action = na.omit))
r.squaredGLMM(s1.no.hym)
AICc(s1.no.hym)

s1.pred.no.hym <- predictInterval(s1.no.hym, new.week, which = "fixed", level = .9)
s1.pred.no.hym <- cbind(data.week$s.year, s1.pred.no.hym) 
names(s1.pred.no.hym) <- c("s.year", "fit", "upr", "lwr")

s1.pred.no.hym$fit <- (s1.pred.no.hym$fit)^2
s1.pred.no.hym$upr <- (s1.pred.no.hym$upr)^2
s1.pred.no.hym$lwr <- (s1.pred.no.hym$lwr)^2
s1.pred.no.hym$taxa <- "hym"

# Full dataset: 1984 - 2020
# Diptera (reported in Table S6)
# Abundance (diptera.no)

summary(s1.no.dip <- lmer(sqrt(diptera.no/total.days) ~ 
                          s.year + s.dps + I(s.dps^2) + s.dps * s.year + 
                          s.mean.temp.ave.cb + s.mean.temp.ave.cb * s.year + 
                          s.sum.precip.cb + s.sum.precip.cb*s.year +
                          s.mean.flw.day + s.mean.flw.day*s.year  + (1|f.year), 
                          data = data.week, na.action = na.omit))
r.squaredGLMM(s1.no.dip)
AICc(s1.no.dip)

s1.pred.no.dip <- predictInterval(s1.no.dip, new.week, which = "fixed", level = .9)
s1.pred.no.dip <- cbind(data.week$s.year, s1.pred.no.dip) 
names(s1.pred.no.dip) <- c("s.year", "fit", "upr", "lwr")

s1.pred.no.dip$fit <- (s1.pred.no.dip$fit)^2
s1.pred.no.dip$upr <- (s1.pred.no.dip$upr)^2
s1.pred.no.dip$lwr <- (s1.pred.no.dip$lwr)^2
s1.pred.no.dip$taxa <- "dip"

# Full dataset: 1984 - 2020
# Other (reported in Table S8)
# Abundance (other.no)

summary(s1.no.other <- lmer(sqrt(other.no/total.days) ~ 
                          s.year + s.dps + I(s.dps^2) + s.dps * s.year + 
                          s.mean.temp.ave.cb + s.mean.temp.ave.cb * s.year + 
                          s.sum.precip.cb + s.sum.precip.cb*s.year +
                          s.mean.flw.day + s.mean.flw.day*s.year + (1|f.year), 
                          data = data.week, na.action = na.omit))
r.squaredGLMM(s1.no.other)
AICc(s1.no.other)

s1.pred.no.other <- predictInterval(s1.no.other, new.week, which = "fixed", level = .9)
s1.pred.no.other <- cbind(data.week$s.year, s1.pred.no.other) 
names(s1.pred.no.other) <- c("s.year", "fit", "upr", "lwr")

s1.pred.no.other$fit <- (s1.pred.no.other$fit)^2
s1.pred.no.other$upr <- (s1.pred.no.other$upr)^2
s1.pred.no.other$lwr <- (s1.pred.no.other$lwr)^2
s1.pred.no.other$taxa <- "other"

s1.no.models <- rbind(s1.pred.no.other, 
                      s1.pred.no.hym, 
                      s1.pred.no.dip)

s1.no.models$type <- "full"

s1.models.no <- rbind(s1.86.no.models, s1.no.models)

##################################################################################
########## Analysis 2: Predict grams/day from subsets of weekly data #############
################  Results are presented in text and Figure S5 ####################
##################################################################################

### FULL DATA SET: 1984 - 2020 ###

# First, run the model from above using all of the data (1984 - 2020)

q1.all <- lmer(sqrt(total.g/total.days) ~  
              s.year + s.dps + I(s.dps^2) + s.dps * s.year + 
              s.mean.temp.ave.cb + s.mean.temp.ave.cb * s.year + 
              s.sum.precip.cb + s.sum.precip.cb*s.year +
              s.mean.flw.day + s.mean.flw.day*s.year + (1|f.year),
              data = data.week, na.action = na.omit, REML = T)  

summary(q1.all)

# Next, create prediction intervals for FIGURE S5

new.week.all <- as.data.frame(matrix(0, ncol = 4 , 
                                     nrow = length(data.week$s.year)))
new.week.all <- cbind(data.week$s.year, new.week.all, data.week$f.year)
names(new.week.all) <- c("s.year", "s.dps", "s.mean.temp.ave.cb", 
                         "s.sum.precip.cb", "s.mean.flw.day", "f.year") 

int.fit.data.all <- predictInterval(q1.all, new.week.all, which = "fixed", level = .9)
pred.all <- cbind(data.week$s.year, data.week$year, int.fit.data.all) 
names(pred.all) <- c("s.year", "year", "fit", "upr", "lwr")

# backtransform, because fit was done on sqrt values

pred.all$fit <- (pred.all$fit)^2
pred.all$upr <- (pred.all$upr)^2
pred.all$lwr <- (pred.all$lwr)^2

### REDUCED DATA SET: 1986 - 2020 ###

# Next, repeat steps for the reduced model (1986 - 2020)

q1.reduced <- lmer(sqrt(total.g/total.days) ~  
                s.year + s.dps + I(s.dps^2) + s.dps * s.year + 
                s.mean.temp.ave.cb + s.mean.temp.ave.cb * s.year + 
                s.sum.precip.cb + s.sum.precip.cb*s.year +
                s.mean.flw.day + s.mean.flw.day*s.year + (1|f.year),
                data = data.week.sub, na.action = na.omit, REML = T)  

summary(q1.reduced)

new.week.reduced <- as.data.frame(matrix(0, ncol = 4 , 
                                         nrow = length(data.week.sub$s.year)))
new.week.reduced <- cbind(data.week.sub$s.year, new.week.reduced, data.week.sub$f.year)
names(new.week.reduced) <- c("s.year", "s.dps", "s.mean.temp.ave.cb", 
                             "s.sum.precip.cb", "s.mean.flw.day", "f.year") 

int.fit.data.reduced <- predictInterval(q1.reduced, new.week.reduced, which = "fixed", level = .9)
pred.reduced <- cbind(data.week.sub$s.year, data.week.sub$year, int.fit.data.reduced) 
names(pred.reduced) <- c("s.year", "year", "fit", "upr", "lwr")

### backtransform, because fit was done on sqrt values
pred.reduced$fit <- (pred.reduced$fit)^2
pred.reduced$upr <- (pred.reduced$upr)^2
pred.reduced$lwr <- (pred.reduced$lwr)^2

### Most recent 10 years of data: 2011 - 2020 ###

# First, we need to rescale the variables

data10 <- subset(data.week, data.week$year > 2010)
data10$s.dps <- c(scale(data10$dps, scale = T, center = T))
data10$s.year <- c(scale(data10$year, scale = T, center = T))
data10$s.mean.temp.ave.cb <- c(scale(data10$mean.temp.ave.cb, scale = T, center = T))
data10$s.sum.precip.cb <- c(scale(data10$sum.precip.cb, scale = T, center = T))
data10$s.mean.flw.day <- c(scale(data10$mean.flw.day, scale = T, center = T))
data10$s.snowpack <- c(scale(data10$snowpack, scale = T, center = T))

q1.10 <- lmer(sqrt(total.g/total.days) ~  
              s.year + s.dps + I(s.dps^2) + s.dps * s.year + 
              s.mean.temp.ave.cb + s.mean.temp.ave.cb * s.year + 
              s.sum.precip.cb + s.sum.precip.cb*s.year +
              s.mean.flw.day + s.mean.flw.day*s.year + (1|f.year),
              data = data10, na.action = na.omit, REML = T)  

summary(q1.10) 

# prediction lines

new.week.10 <- as.data.frame(matrix(0, ncol = 4 , 
                                    nrow = length(data10$s.year)))
new.week.10 <- cbind(data10$s.year, new.week.10, data10$f.year)
names(new.week.10) <- c("s.year", "s.dps", "s.mean.temp.ave.cb", 
                        "s.sum.precip.cb", "s.mean.flw.day", "f.year") 
int.fit.data10 <- predictInterval(q1.10, new.week.10, which = "fixed", level = .9)
pred.10 <- cbind(data10$s.year, data10$year, int.fit.data10) 
names(pred.10) <- c("s.year", "year", "fit", "upr", "lwr")

# backtransform, because fit was done on sqrt values

pred.10$fit <- (pred.10$fit)^2
pred.10$upr <- (pred.10$upr)^2
pred.10$lwr <- (pred.10$lwr)^2

### Most recent 15 years of data: 2006 - 2020 ###

data15 <- subset(data.week, data.week$year > 2005)
data15$s.dps <- c(scale(data15$dps, scale = T, center = T))
data15$s.year <- c(scale(data15$year, scale = T, center = T))
data15$s.mean.temp.ave.cb <- c(scale(data15$mean.temp.ave.cb, scale = T, center = T))
data15$s.sum.precip.cb <- c(scale(data15$sum.precip.cb, scale = T, center = T))
data15$s.mean.flw.day <- c(scale(data15$mean.flw.day, scale = T, center = T))
data15$s.snowpack <- c(scale(data15$snowpack, scale = T, center = T))

q1.15 <- lmer(sqrt(total.g/total.days) ~  
              s.year + s.dps + I(s.dps^2) + s.dps * s.year + 
              s.mean.temp.ave.cb + s.mean.temp.ave.cb * s.year + 
              s.sum.precip.cb + s.sum.precip.cb*s.year +
              s.mean.flw.day + s.mean.flw.day*s.year + (1|f.year),
              data = data15, na.action = na.omit, REML = T) 
summary(q1.15)

# prediction lines

new.week.15 <- as.data.frame(matrix(0, ncol = 4 , 
                                    nrow = length(data15$s.year)))
new.week.15 <- cbind(data15$s.year, new.week.15, data15$f.year)
names(new.week.15) <- c("s.year", "s.dps", "s.mean.temp.ave.cb", 
                        "s.sum.precip.cb", "s.mean.flw.day", "f.year") 
int.fit.data15 <- predictInterval(q1.15, new.week.15, which = "fixed", level = .9)
pred.15 <- cbind(data15$s.year, data15$year, int.fit.data15) 
names(pred.15) <- c("s.year", "year", "fit", "upr", "lwr")

# backtransform, because fit was done on sqrt values

pred.15$fit <- (pred.15$fit)^2
pred.15$upr <- (pred.15$upr)^2
pred.15$lwr <- (pred.15$lwr)^2

### Most recent 20 years of data: 2001 - 2020 ###

data20 <- subset(data.week, data.week$year > 2000)
data20$s.dps <- c(scale(data20$dps, scale = T, center = T))
data20$s.year <- c(scale(data20$year, scale = T, center = T))
data20$s.mean.temp.ave.cb <- c(scale(data20$mean.temp.ave.cb, scale = T, center = T))
data20$s.sum.precip.cb <- c(scale(data20$sum.precip.cb, scale = T, center = T))
data20$s.mean.flw.day <- c(scale(data20$mean.flw.day, scale = T, center = T))
data20$s.snowpack <- c(scale(data20$snowpack, scale = T, center = T))

q1.20 <- lmer(sqrt(total.g/total.days) ~  
              s.year + s.dps + I(s.dps^2) + s.dps * s.year + 
              s.mean.temp.ave.cb + s.mean.temp.ave.cb * s.year + 
              s.sum.precip.cb + s.sum.precip.cb*s.year +
              s.mean.flw.day + s.mean.flw.day*s.year + (1|f.year),
              data = data20, na.action = na.omit, REML = T) 

summary(q1.20) 

# prediction lines

new.week.20 <- as.data.frame(matrix(0, ncol = 4 , 
                                    nrow = length(data20$s.year)))
new.week.20 <- cbind(data20$s.year, new.week.20, data20$f.year)
names(new.week.20) <- c("s.year", "s.dps", "s.mean.temp.ave.cb", 
                        "s.sum.precip.cb", "s.mean.flw.day", "f.year") 
int.fit.data20 <- predictInterval(q1.20, new.week.20, which = "fixed", level = .9)
pred.20 <- cbind(data20$s.year, data20$year, int.fit.data20) 
names(pred.20) <- c("s.year", "year", "fit", "upr", "lwr")

# backtransform, because fit was done on sqrt values

pred.20$fit <- (pred.20$fit)^2
pred.20$upr <- (pred.20$upr)^2
pred.20$lwr <- (pred.20$lwr)^2

### Early 10 years of data: 1986 - 1995 ###

data10b <- subset(data.week, data.week$year > 1985 & data.week$year < 1996)
data10b$s.dps <- c(scale(data10b$dps, scale = T, center = T))
data10b$s.year <- c(scale(data10b$year, scale = T, center = T))
data10b$s.mean.temp.ave.cb <- c(scale(data10b$mean.temp.ave.cb, scale = T, center = T))
data10b$s.sum.precip.cb <- c(scale(data10b$sum.precip.cb, scale = T, center = T))
data10b$s.mean.flw.day <- c(scale(data10b$mean.flw.day, scale = T, center = T))
data10b$s.snowpack <- c(scale(data10b$snowpack, scale = T, center = T))

q1.10b <- lmer(sqrt(total.g/total.days) ~  
              s.year + s.dps + I(s.dps^2) + s.dps * s.year + 
              s.mean.temp.ave.cb + s.mean.temp.ave.cb * s.year + 
              s.sum.precip.cb + s.sum.precip.cb*s.year +
              s.mean.flw.day + s.mean.flw.day*s.year + (1|f.year),
              data = data10b, na.action = na.omit, REML = T) 
summary(q1.10b)

##################################################################################
######## Analysis 3: Predict grams/year or number/year from yearly data  #########
###########  Results are presented in Figs 2 & 3 and Tables S9 - S12 #############
##################################################################################

# Calculate the residuals for the reduced data set for biomass and abundance
# First, compare using total.g/total.days to sqrt transformation

mean.dps.total.g.86.lm <-  lm((total.g/total.days) ~ 
                              I(s.dps^2) + s.dps + s.mean.temp.ave.cb + 
                              s.sum.precip.cb + s.mean.flw.day, 
                              data = data.week.sub, na.action = na.omit)

summary(mean.dps.total.g.86.lm)
r.squaredGLMM(mean.dps.total.g.86.lm)

mean.dps.total.g.86.lmsqrt <- lm(sqrt(total.g/total.days) ~ 
                                I(s.dps^2) + s.dps + s.mean.temp.ave.cb + 
                                s.sum.precip.cb + s.mean.flw.day, 
                                data = data.week.sub, na.action = na.omit)

summary(mean.dps.total.g.86.lmsqrt)
r.squaredGLMM(mean.dps.total.g.86.lmsqrt)

AIC(mean.dps.total.g.86.lm, mean.dps.total.g.86.lmsqrt) #sqrt is better
mean.dps.total.g.86 <- mean.dps.total.g.86.lmsqrt

pred.dps.total.g.86 <- predict(mean.dps.total.g.86, data = data.week.sub)
dps.total.g.86 <- data.week.sub[!is.na(data.week.sub$total.g), ]
dps.total.g.86 <- dps.total.g.86[!is.na(dps.total.g.86$s.mean.flw.day), ]
dps.total.g.86$resid <- resid(mean.dps.total.g.86)

dps.total.g.86.summary <- as.data.frame(dps.total.g.86 %>% 
                                          group_by(year) %>% 
                                          summarize(sum.resid.total.g = sum(resid)))

mean.dps.total.no.86.lm <- lm((total.no/total.days) ~ 
                              I(s.dps^2) + s.dps + s.mean.temp.ave.cb + 
                              s.sum.precip.cb + s.mean.flw.day,
                              data = data.week.sub, na.action = na.omit)

summary(mean.dps.total.no.86.lm)
r.squaredGLMM(mean.dps.total.no.86.lm)

mean.dps.total.no.86.lmsqrt <- lm(sqrt(total.no/total.days) ~ 
                                I(s.dps^2) + s.dps + s.mean.temp.ave.cb +
                                s.sum.precip.cb + s.mean.flw.day,
                                data = data.week.sub, na.action = na.omit)

summary(mean.dps.total.no.86.lmsqrt)
r.squaredGLMM(mean.dps.total.no.86.lmsqrt)

AIC(mean.dps.total.no.86.lm, mean.dps.total.no.86.lmsqrt) #sqrt is better
mean.dps.total.no.86 <- mean.dps.total.no.86.lmsqrt

dps.total.no.86 <- data.week.sub[!is.na(data.week.sub$total.no), ]
dps.total.no.86 <- dps.total.no.86[!is.na(dps.total.no.86$s.mean.flw.day), ]
dps.total.no.86$resid <- resid(mean.dps.total.no.86)

dps.total.no.86.summary <- as.data.frame(dps.total.no.86 %>% 
                                           group_by(year) %>% 
                                           summarize(sum.resid.total.no = sum(resid)))

total.resid <- merge(dps.total.g.86.summary, dps.total.no.86.summary, by = c("year"))

data.year.sub <- merge(data.year.sub, total.resid, by = c("year"))
data.year.sub <- data.year.sub[!is.na(data.year.sub$s.lag.yr.floralcount), ]
data.year.sub <- data.year.sub[!is.na(data.year.sub$s.yr.floralcount), ]

# Run global models predicting residuals based on weather
# Abundance: reduced data (1986 - 2020)
# Results presented in Figs 2 and 3 and Table S9

q2.86.total.no <- lm(sum.resid.total.no ~ s.year  + 
                    s.summer.temp.cb + s.lag.summer.temp.cb +
                    s.summer.precip.cb + s.lag.summer.precip.cb +
                    s.snowpack + s.lag.snowpack + 
                    s.yr.floralcount + s.lag.yr.floralcount + 
                    s.year * s.summer.temp.cb + s.year * s.lag.summer.temp.cb +
                    s.year * s.summer.precip.cb + s.year * s.lag.summer.precip.cb +
                    s.year * s.snowpack + s.year * s.lag.snowpack +
                    s.year * s.yr.floralcount + s.year * s.lag.yr.floralcount,
                    weights = sqrt(yr.samples), 
                    data = data.year.sub, na.action = na.omit)

r.squaredGLMM(q2.86.total.no)
vif(q2.86.total.no) # first removed s.year:s.lag.yr.floralcount:    16.257401  

q2.86.total.no <- lm(sum.resid.total.no ~ s.year +
                    s.summer.temp.cb + s.lag.summer.temp.cb +
                    s.summer.precip.cb + s.lag.summer.precip.cb +
                    s.snowpack + s.lag.snowpack + 
                    s.yr.floralcount + s.lag.yr.floralcount + 
                    s.year * s.summer.temp.cb + s.year * s.lag.summer.temp.cb +
                    s.year * s.summer.precip.cb + s.year * s.lag.summer.precip.cb +
                    s.year * s.snowpack + s.year * s.lag.snowpack +
                    s.year * s.yr.floralcount,
                    weights = sqrt(yr.samples), 
                    data = data.year.sub, na.action = na.omit)

vif(q2.86.total.no) # next removed : s.year:s.yr.floralcount: vif =     5.633738   

q2.86.total.no <- lm(sum.resid.total.no ~ s.year +
                    s.summer.temp.cb + s.lag.summer.temp.cb +
                    s.summer.precip.cb + s.lag.summer.precip.cb +
                    s.snowpack + s.lag.snowpack + 
                    s.yr.floralcount + s.lag.yr.floralcount + 
                    s.year * s.summer.temp.cb + s.year * s.lag.summer.temp.cb +
                    s.year * s.summer.precip.cb + s.year * s.lag.summer.precip.cb +
                    s.year * s.snowpack + s.year * s.lag.snowpack,
                    weights = sqrt(yr.samples), 
                    data = data.year.sub, na.action = na.omit)

vif(q2.86.total.no) # next removed : s.year:s.lag.summer.temp.cb  vif = 4.268112    

q2.86.total.no <- lm(sum.resid.total.no ~ s.year +
                    s.summer.temp.cb + s.lag.summer.temp.cb +
                    s.summer.precip.cb + s.lag.summer.precip.cb +
                    s.snowpack + s.lag.snowpack + 
                    s.yr.floralcount + s.lag.yr.floralcount + 
                    s.year * s.summer.temp.cb + 
                    s.year * s.summer.precip.cb + s.year * s.lag.summer.precip.cb +
                    s.year * s.snowpack + s.year * s.lag.snowpack,
                    weights = sqrt(yr.samples), 
                    data = data.year.sub, na.action = na.omit)

vif(q2.86.total.no) # next removed : s.year:s.lag.snowpack   vif =    3.488182   

q2.86.total.no <- lm(sum.resid.total.no ~ s.year +
                    s.summer.temp.cb + s.lag.summer.temp.cb +
                    s.summer.precip.cb + s.lag.summer.precip.cb +
                    s.snowpack + s.lag.snowpack + 
                    s.yr.floralcount + s.lag.yr.floralcount + 
                    s.year * s.summer.precip.cb + s.year * s.lag.summer.precip.cb +
                    s.year * s.snowpack,
                    weights = sqrt(yr.samples), 
                    data = data.year.sub, na.action = na.omit)

vif(q2.86.total.no) # VIF are all less than 5
summary(q2.86.total.no)

# Use a step approach to select the best fit model

summary(q2.86.selected.total.no <- step(q2.86.total.no, direction = "both",
                                     scope = list(lower = ~s.year)))
confint(q2.86.selected.total.no)
# used extractAIC() here instead of AICc because of AIC used in step() function
extractAIC(q2.86.selected.total.no, scale = 0) 

# Run global models predicting residuals based on weather
# Biomass: reduced data (1986 - 2020)
# Results presented in Fig 2 and Table S9

q2.86.total.g <- lm(sum.resid.total.g ~ s.year +
                    s.summer.temp.cb + s.lag.summer.temp.cb +
                    s.summer.precip.cb + s.lag.summer.precip.cb +
                    s.snowpack + s.lag.snowpack + 
                    s.yr.floralcount + s.lag.yr.floralcount + 
                    s.year * s.summer.temp.cb + s.year * s.lag.summer.temp.cb +
                    s.year * s.summer.precip.cb + s.year * s.lag.summer.precip.cb +
                    s.year * s.snowpack + s.year * s.lag.snowpack +
                    s.year * s.yr.floralcount + s.year * s.lag.yr.floralcount,
                    weights = sqrt(yr.samples), na.action = na.omit, 
                    data = data.year.sub)

vif(q2.86.total.g) # removed same factors as we did for abundance  

q2.86.total.g <- lm(sum.resid.total.g ~ s.year +
                    s.summer.temp.cb + s.lag.summer.temp.cb +
                    s.summer.precip.cb + s.lag.summer.precip.cb +
                    s.snowpack + s.lag.snowpack + 
                    s.yr.floralcount + s.lag.yr.floralcount + 
                    s.year * s.summer.precip.cb + s.year * s.lag.summer.precip.cb +
                    s.year * s.snowpack + s.year * s.lag.snowpack,
                    weights = sqrt(yr.samples), 
                    data = data.year.sub, na.action = na.omit)
vif(q2.86.total.g) 

summary(q2.86.selected.total.g <- step(q2.86.total.g, direction = "both"))              
confint(q2.86.selected.total.g)  
extractAIC(q2.86.selected.total.g, scale = 0) 

# Calculate the residuals for the full data set for biomass and abundance

mean.dps.total.g <-  lm(sqrt(total.g/total.days) ~ I(s.dps^2) + s.dps + 
                        s.mean.temp.ave.cb + s.sum.precip.cb + s.mean.flw.day,
                        data = data.week, na.action = na.omit)

summary(mean.dps.total.g)

pred.dps.total.g <- predict(mean.dps.total.g, data=data.week)
dps.total.g <- data.week[!is.na(data.week$total.g), ]
dps.total.g <- dps.total.g[!is.na(dps.total.g$mean.flw.day), ]
dps.total.g$resid <- resid(mean.dps.total.g)

dps.total.g.summary <- as.data.frame(dps.total.g %>% 
                                     group_by(year) %>% 
                                     summarize(sum.resid.total.g = sum(resid)))


mean.dps.total.no <- lm(sqrt(total.no/total.days) ~ I(s.dps^2) + s.dps + 
                        s.mean.temp.ave.cb + s.sum.precip.cb + s.mean.flw.day,
                        data = data.week, na.action = na.omit)

summary(mean.dps.total.no)
dps.total.no <- data.week[!is.na(data.week$total.no), ]
dps.total.no <- dps.total.no[!is.na(dps.total.no$mean.flw.day), ]
dps.total.no$resid <- resid(mean.dps.total.no)

dps.total.no.summary <- as.data.frame(dps.total.no %>% 
                                      group_by(year) %>% 
                                      summarize(sum.resid.total.no = sum(resid)))

total.resid <- merge(dps.total.g.summary, dps.total.no.summary, by = c("year"))

data.year <- merge(data.year, total.resid, by = c("year"))
data.year <- data.year[!is.na(data.year$s.lag.yr.floralcount), ]
data.year <- data.year[!is.na(data.year$s.yr.floralcount), ]

# Run global models predicting residuals based on weather
# Abundance: full data (1984 - 2020)
# Results presented in Table S10

q2.total.no <- lm(sum.resid.total.no ~ s.year  + 
                  s.summer.temp.cb + s.lag.summer.temp.cb +
                  s.summer.precip.cb + s.lag.summer.precip.cb +
                  s.snowpack + s.lag.snowpack + 
                  s.yr.floralcount + s.lag.yr.floralcount + 
                  s.year * s.summer.temp.cb + s.year * s.lag.summer.temp.cb +
                  s.year * s.summer.precip.cb + s.year * s.lag.summer.precip.cb +
                  s.year * s.snowpack + s.year * s.lag.snowpack +
                  s.year * s.yr.floralcount + s.year * s.lag.yr.floralcount,
                  weights = sqrt(yr.samples), 
                  data = data.year, na.action = na.omit)

vif(q2.total.no) # first removed s.year:s.lag.yr.floralcount, vif = 18.733455 

q2.total.no <- lm(sum.resid.total.no ~ s.year +
                  s.summer.temp.cb + s.lag.summer.temp.cb +
                  s.summer.precip.cb + s.lag.summer.precip.cb +
                  s.snowpack + s.lag.snowpack + 
                  s.yr.floralcount + s.lag.yr.floralcount + 
                  s.year * s.summer.temp.cb + s.year * s.lag.summer.temp.cb +
                  s.year * s.summer.precip.cb + s.year * s.lag.summer.precip.cb +
                  s.year * s.snowpack + s.year * s.lag.snowpack +
                  s.year * s.yr.floralcount,
                  weights = sqrt(yr.samples), 
                  data = data.year, na.action = na.omit)

vif(q2.total.no) # next, removed s.year : s.yr.floralcount, vif =     5.793368  

q2.total.no <- lm(sum.resid.total.no ~ s.year +
                  s.summer.temp.cb + s.lag.summer.temp.cb +
                  s.summer.precip.cb + s.lag.summer.precip.cb +
                  s.snowpack + s.lag.snowpack + 
                  s.yr.floralcount + s.lag.yr.floralcount + 
                  s.year * s.summer.temp.cb + s.year * s.lag.summer.temp.cb +
                  s.year * s.summer.precip.cb + s.year * s.lag.summer.precip.cb +
                  s.year * s.snowpack + s.year * s.lag.snowpack,
                  weights = sqrt(yr.samples), 
                  data = data.year, na.action = na.omit)

vif(q2.total.no) # next, removed s.year:s.lag.snowpack:      2.450458     

q2.total.no <- lm(sum.resid.total.no ~ s.year +
                  s.summer.temp.cb + s.lag.summer.temp.cb +
                  s.summer.precip.cb + s.lag.summer.precip.cb +
                  s.snowpack + s.lag.snowpack + 
                  s.yr.floralcount + s.lag.yr.floralcount + 
                  s.year * s.summer.temp.cb  + s.year * s.lag.summer.temp.cb +
                  s.year * s.summer.precip.cb + s.year * s.lag.summer.precip.cb +
                  s.year * s.snowpack,
                  weights = sqrt(yr.samples), 
                  data = data.year, na.action = na.omit)

vif(q2.total.no) #  next, removed s.year:s.lag.summer.temp.cb:     1.737203  

q2.total.no <- lm(sum.resid.total.no ~ s.year +
                  s.summer.temp.cb + s.lag.summer.temp.cb +
                  s.summer.precip.cb + s.lag.summer.precip.cb +
                  s.snowpack + s.lag.snowpack + 
                  s.yr.floralcount + s.lag.yr.floralcount + 
                  s.year * s.summer.temp.cb  +
                  s.year * s.summer.precip.cb + s.year * s.lag.summer.precip.cb +
                  s.year * s.snowpack,
                  weights = sqrt(yr.samples), data = data.year, na.action = na.omit)

vif(q2.total.no) # all below VIF of 5

# Use a step approach to select the best fit model
summary(q2.total.no) 
summary(q2.selected.total.no <- step(q2.total.no, direction = "both", scope = list(lower = ~s.year)))
confint(q2.selected.total.no)
extractAIC(q2.selected.total.no, scale = 0)

# Run global models predicting residuals based on weather
# Biomass: full data (1984 - 2020)
# Results presented in Table S10

q2.total.g <- lm(sum.resid.total.g ~ s.year +
                  s.summer.temp.cb + s.lag.summer.temp.cb +
                  s.summer.precip.cb + s.lag.summer.precip.cb +
                  s.snowpack + s.lag.snowpack + 
                  s.yr.floralcount + s.lag.yr.floralcount + 
                  s.year * s.summer.temp.cb + s.year * s.lag.summer.temp.cb +
                  s.year * s.summer.precip.cb + s.year * s.lag.summer.precip.cb +
                  s.year * s.snowpack + s.year * s.lag.snowpack +
                  s.year * s.yr.floralcount + s.year * s.lag.yr.floralcount,
                  weights = sqrt(yr.samples), 
                  na.action = na.omit, data = data.year)

# Removed same covariates as total.no for full data set

q2.total.g <- lm(sum.resid.total.g ~ s.year +
                  s.summer.temp.cb + s.lag.summer.temp.cb +
                  s.summer.precip.cb + s.lag.summer.precip.cb +
                  s.snowpack + s.lag.snowpack + 
                  s.yr.floralcount + s.lag.yr.floralcount + 
                  s.year * s.summer.temp.cb  +
                  s.year * s.summer.precip.cb + s.year * s.lag.summer.precip.cb +
                  s.year * s.snowpack, na.action = na.omit,
                  weights = sqrt(yr.samples), data = data.year)

vif(q2.total.g) #  all VIFs below 5

summary(q2.selected.total.g <- step(q2.total.g, direction = "both", scope = list(lower = ~s.year)))                        
vif(q2.selected.total.g)
confint(q2.selected.total.g)
extractAIC(q2.selected.total.g, scale = 0)

# Calculate the residuals for the yearly data set for biomass and abundance
# without 2015 and 2016 - high floral abundance in 2015 and 2016 has lagged
# floral abundance

# First, subset the df and rescale variables without 2015 and 2016
data.year.out <- subset(data.year, year != 2015 & year != 2016)
data.year.sub.out <- subset(data.year.out, year != 1984 & year != 1985) 

data.year.out$s.year <- c(scale(data.year.out$year, center = T, scale = T))
data.year.out$s.summer.precip.cb <- c(scale(data.year.out$summer.precip.cb, center = T, scale = T))
data.year.out$s.summer.temp.cb <- c(scale(data.year.out$summer.temp.cb, center = T, scale = T))
data.year.out$s.snowmelt <- c(scale(data.year.out$snowmelt, center = T, scale = T))
data.year.out$s.lag.summer.precip.cb <- c(scale(data.year.out$lag.summer.precip.cb, center = T, scale = T))
data.year.out$s.lag.summer.temp.cb <- c(scale(data.year.out$lag.summer.temp.cb, center = T, scale = T))
data.year.out$s.snowpack <- c(scale(data.year.out$snowpack, center = T, scale = T))
data.year.out$s.lag.snowpack <- c(scale(data.year.out$lag.snowpack, center = T, scale = T))
data.year.out$s.yr.floralcount <- c(scale(data.year.out$yr.floralcount, center = T, scale = T))
data.year.out$s.lag.yr.floralcount <- c(scale(data.year.out$lag.yr.floralcount, center = T, scale = T))

data.year.sub.out$s.year <- c(scale(data.year.sub.out$year, center = T, scale = T))
data.year.sub.out$s.summer.precip.cb <- c(scale(data.year.sub.out$summer.precip.cb, center = T, scale = T))
data.year.sub.out$s.summer.temp.cb <- c(scale(data.year.sub.out$summer.temp.cb, center = T, scale = T))
data.year.sub.out$s.snowmelt <- c(scale(data.year.sub.out$snowmelt, center = T, scale = T))
data.year.sub.out$s.lag.summer.precip.cb <- c(scale(data.year.sub.out$lag.summer.precip.cb, center = T, scale = T))
data.year.sub.out$s.lag.summer.temp.cb <- c(scale(data.year.sub.out$lag.summer.temp.cb, center = T, scale = T))
data.year.sub.out$s.snowpack <- c(scale(data.year.sub.out$snowpack, center = T, scale = T))
data.year.sub.out$s.lag.snowpack <- c(scale(data.year.sub.out$lag.snowpack, center = T, scale = T))
data.year.sub.out$s.yr.floralcount <- c(scale(data.year.sub.out$yr.floralcount, center = T, scale = T))
data.year.sub.out$s.lag.yr.floralcount <- c(scale(data.year.sub.out$lag.yr.floralcount, center = T, scale = T))

## 1. Calculate the residuals for the reduced data set (1986 - 2020)
data.week.out <- subset(data.week, year != 2015 & year != 2016)
data.week.sub.out <- subset(data.week.out, year != 1984  & year != 1985)

# scale the variables
# full dataset (1984 - 2020):
data.week.out$s.year <- c(scale(data.week.out$year, scale = T, center = T))
data.week.out$s.dps <- c(scale(data.week.out$dps, scale = T, center = T))
data.week.out$s.mean.temp.ave.cb <- c(scale(data.week.out$mean.temp.ave.cb, scale = T, center = T))
data.week.out$s.sum.precip.cb <- c(scale(data.week.out$sum.precip.cb, scale = T, center = T))
data.week.out$s.mean.flw.day <- c(scale(data.week.out$mean.flw.day, scale = T, center = T))
data.week.out$s.snowpack <- c(scale(data.week.out$snowpack, scale = T, center = T))

# reduced dataset (1986 - 2020):
data.week.sub.out$s.year <- c(scale(data.week.sub.out$year, scale = T, center = T))
data.week.sub.out$s.dps <- c(scale(data.week.sub.out$dps, scale = T, center = T))
data.week.sub.out$s.mean.temp.ave.cb <- c(scale(data.week.sub.out$mean.temp.ave.cb, scale = T, center = T))
data.week.sub.out$s.sum.precip.cb <- c(scale(data.week.sub.out$sum.precip.cb, scale = T, center = T))
data.week.sub.out$s.mean.flw.day <- c(scale(data.week.sub.out$mean.flw.day, scale = T, center = T))
data.week.sub.out$s.snowpack <- c(scale(data.week.sub.out$snowpack, scale = T, center = T))

mean.dps.total.g.86 <-  lm(sqrt(total.g/total.days) ~ I(s.dps^2) + s.dps + 
                             s.mean.temp.ave.cb + s.sum.precip.cb + s.mean.flw.day,
                           data = data.week.sub.out, na.action = na.omit)

#with(data.week.sub.out, plot(sqrt(total.g/total.days) ~ s.dps))  # w, w/o sqrt 
summary(mean.dps.total.g.86)

pred.dps.total.g.86 <- predict(mean.dps.total.g.86, data = data.week.sub.out)
dps.total.g.86 <- data.week.sub.out[!is.na(data.week.sub.out$total.g), ]
dps.total.g.86 <- dps.total.g.86[!is.na(dps.total.g.86$s.mean.flw.day), ]
dps.total.g.86$resid <- resid(mean.dps.total.g.86)

dps.total.g.86.summary <- as.data.frame(dps.total.g.86 %>% 
                                          group_by(year) %>% 
                                          summarize(sum.resid.total.g = sum(resid)))


mean.dps.total.no.86 <- lm(sqrt(total.no/total.days) ~ I(s.dps^2) + s.dps + 
                             s.mean.temp.ave.cb + s.sum.precip.cb + s.mean.flw.day,
                           data = data.week.sub.out, na.action = na.omit)

summary(mean.dps.total.no.86)
dps.total.no.86 <- data.week.sub.out[!is.na(data.week.sub.out$total.no), ]
dps.total.no.86 <- dps.total.no.86[!is.na(dps.total.no.86$s.mean.flw.day), ]
dps.total.no.86$resid <- resid(mean.dps.total.no.86)

dps.total.no.86.summary <- as.data.frame(dps.total.no.86 %>% 
                                           group_by(year) %>% 
                                           summarize(sum.resid.total.no = sum(resid)))

total.resid <- merge(dps.total.g.86.summary, dps.total.no.86.summary, by = c("year"))

data.year.sub.out <- merge(data.year.sub.out, total.resid, by = c("year"))
data.year.sub.out <- data.year.sub.out[!is.na(data.year.sub.out$s.lag.yr.floralcount), ]
data.year.sub.out <- data.year.sub.out[!is.na(data.year.sub.out$s.yr.floralcount), ]

# 2. Run global models predicting residuals based on weather

q2.86.total.no <- lm(sum.resid.total.no ~ s.year  + 
                       s.summer.temp.cb + s.lag.summer.temp.cb +
                       s.summer.precip.cb + s.lag.summer.precip.cb +
                       s.snowpack + s.lag.snowpack + 
                       s.yr.floralcount + s.lag.yr.floralcount + 
                       s.year * s.summer.temp.cb + s.year * s.lag.summer.temp.cb +
                       s.year * s.summer.precip.cb + s.year * s.lag.summer.precip.cb +
                       s.year * s.snowpack + s.year * s.lag.snowpack +
                       s.year * s.yr.floralcount + s.year * s.lag.yr.floralcount,
                     weights = sqrt(yr.samples), 
                     data = data.year.sub.out, na.action = na.omit)

vif(q2.86.total.no) # first removed s.year:s.lag.yr.floralcount:     5.187155 
vif(q2.86.total.no) # next removed : s.year:s.yr.floralcount: vif =    4.627527   

q2.86.total.no <- lm(sum.resid.total.no ~ s.year +
                       s.summer.temp.cb + s.lag.summer.temp.cb +
                       s.summer.precip.cb + s.lag.summer.precip.cb +
                       s.snowpack + s.lag.snowpack + 
                       s.yr.floralcount + s.lag.yr.floralcount + 
                       s.year * s.summer.temp.cb + s.year * s.lag.summer.temp.cb +
                       s.year * s.summer.precip.cb + s.year * s.lag.summer.precip.cb +
                       s.year * s.snowpack + s.year * s.lag.snowpack,
                     weights = sqrt(yr.samples), data = data.year.sub.out, na.action = na.omit)
vif(q2.86.total.no) # next removed : s.year:s.lag.summer.temp.cb  vif =    4.704864    

q2.86.total.no <- lm(sum.resid.total.no ~ s.year +
                       s.summer.temp.cb + s.lag.summer.temp.cb +
                       s.summer.precip.cb + s.lag.summer.precip.cb +
                       s.snowpack + s.lag.snowpack + 
                       s.yr.floralcount + s.lag.yr.floralcount + 
                       s.year * s.summer.temp.cb + 
                       s.year * s.summer.precip.cb + s.year * s.lag.summer.precip.cb +
                       s.year * s.snowpack + s.year * s.lag.snowpack,
                     weights = sqrt(yr.samples), data = data.year.sub.out, na.action = na.omit)
vif(q2.86.total.no) #  next removed : s.year:s.summer.temp.cb    vif =           3.795924        

q2.86.total.no <- lm(sum.resid.total.no ~ s.year +
                       s.summer.temp.cb + 
                       s.summer.precip.cb + s.lag.summer.precip.cb +
                       s.snowpack + s.lag.snowpack + 
                       s.yr.floralcount + s.lag.yr.floralcount + 
                       s.year * s.summer.temp.cb + 
                       s.year * s.summer.precip.cb + s.year * s.lag.summer.precip.cb +
                       s.year * s.snowpack + s.year * s.lag.snowpack,
                     weights = sqrt(yr.samples), data = data.year.sub.out, na.action = na.omit)
vif(q2.86.total.no) #  next removed : s.year:s.summer.temp.cb   vif =  3.463058          

q2.86.total.no <- lm((sum.resid.total.no/yr.samples) ~ s.year + 
                       s.summer.temp.cb + 
                       s.summer.precip.cb + s.lag.summer.precip.cb +
                       s.snowpack + s.lag.snowpack + 
                       s.yr.floralcount + s.lag.yr.floralcount + 
                       s.year * s.summer.precip.cb + s.year * s.lag.summer.precip.cb +
                       s.year * s.snowpack + s.year * s.lag.snowpack,
                     weights = sqrt(yr.samples), data = data.year.sub.out, na.action = na.omit)

vif(q2.86.total.no) # VIF are all less than 3, except year and summer temp forcing into model

summary(q2.86.total.no) 

# 3. Use a step approach to select the best fit model
summary(q2.selected.total.no <- step(q2.86.total.no, direction = "both",
                                     scope = list(lower = ~s.year)))   #temp, precip and lag snowpack
extractAIC(q2.selected.total.no, scale = 0)

# 4. Repeat same process for total insect biomass (global model, stepwise)
# vif(q2.86.total.g) # removed same as for number  

q2.86.total.g <- lm((sum.resid.total.g/yr.samples) ~ s.year +
                      s.summer.temp.cb + s.lag.summer.temp.cb +
                      s.summer.precip.cb + s.lag.summer.precip.cb +
                      s.snowpack + s.lag.snowpack + 
                      s.yr.floralcount + s.lag.yr.floralcount + 
                      s.year * s.summer.precip.cb + s.year * s.lag.summer.precip.cb +
                      s.year * s.snowpack + s.year * s.lag.snowpack,
                    weights = sqrt(yr.samples), data = data.year.sub.out, na.action = na.omit)
vif(q2.86.total.g) 

summary(q2.selected.total.g <- step(q2.86.total.g, direction = "both"))  # precip, lag precip, floralcount, marginal lag floralcount
vif(q2.selected.total.g)
extractAIC(q2.selected.total.g)

################################################################################

## 1. Calculate the residuals for the full data set (1984 - 2020) without
## 2015 which had a high floral abundance and 2016 which had a high lag 
## floral abundance

mean.dps.total.g <-  lm(sqrt(total.g/total.days) ~ I(s.dps^2) + s.dps + 
                          s.mean.temp.ave.cb + s.sum.precip.cb + s.mean.flw.day,
                        data = data.week.out, na.action = na.omit)

summary(mean.dps.total.g)

pred.dps.total.g <- predict(mean.dps.total.g, data=data.week.out)
dps.total.g <- data.week.out[!is.na(data.week.out$total.g), ]
dps.total.g <- dps.total.g[!is.na(dps.total.g$mean.flw.day), ]
dps.total.g$resid <- resid(mean.dps.total.g)

dps.total.g.summary <- as.data.frame(dps.total.g %>% 
                                       group_by(year) %>% 
                                       summarize(sum.resid.total.g = sum(resid)))


mean.dps.total.no <- lm(sqrt(total.no/total.days) ~ I(s.dps^2) + s.dps + 
                          s.mean.temp.ave.cb + s.sum.precip.cb + s.mean.flw.day,
                        data = data.week.out, na.action = na.omit)

summary(mean.dps.total.no)
dps.total.no <- data.week.out[!is.na(data.week.out$total.no), ]
dps.total.no <- dps.total.no[!is.na(dps.total.no$mean.flw.day), ]
dps.total.no$resid <- resid(mean.dps.total.no)

dps.total.no.summary <- as.data.frame(dps.total.no %>% 
                                        group_by(year) %>% 
                                        summarize(sum.resid.total.no = sum(resid)))

total.resid <- merge(dps.total.g.summary, dps.total.no.summary, by = c("year"))

data.year.out <- merge(data.year.out, total.resid, by = c("year"))
data.year.out <- data.year.out[!is.na(data.year.out$s.lag.yr.floralcount), ]
data.year.out <- data.year.out[!is.na(data.year.out$s.yr.floralcount), ]

# 2. Run global models predicting residuals based on weather
q2.total.no <- lm(sum.resid.total.no ~ s.year  + 
                    s.summer.temp.cb + s.lag.summer.temp.cb +
                    s.summer.precip.cb + s.lag.summer.precip.cb +
                    s.snowpack + s.lag.snowpack + 
                    s.yr.floralcount + s.lag.yr.floralcount + 
                    s.year * s.summer.temp.cb + s.year * s.lag.summer.temp.cb +
                    s.year * s.summer.precip.cb + s.year * s.lag.summer.precip.cb +
                    s.year * s.snowpack + s.year * s.lag.snowpack +
                    s.year * s.yr.floralcount + s.year * s.lag.yr.floralcount,
                  weights = sqrt(yr.samples), 
                  data = data.year.out, na.action = na.omit)

vif(q2.total.no) # first removed s.year:s.snowpack, vif =  6.694982   

q2.total.no <- lm(sum.resid.total.no ~ s.year  + 
                    s.summer.temp.cb + s.lag.summer.temp.cb +
                    s.summer.precip.cb + s.lag.summer.precip.cb +
                    s.snowpack + s.lag.snowpack + 
                    s.yr.floralcount + s.lag.yr.floralcount + 
                    s.year * s.summer.temp.cb + s.year * s.lag.summer.temp.cb +
                    s.year * s.summer.precip.cb + s.year * s.lag.summer.precip.cb +
                    s.year * s.lag.snowpack +
                    s.year * s.yr.floralcount + s.year * s.lag.yr.floralcount,
                  weights = sqrt(yr.samples), 
                  data = data.year.out, na.action = na.omit)
vif(q2.total.no) # next, removed s.year : s.yr.floralcount, vif =  3.328974  

q2.total.no <- lm(sum.resid.total.no ~ s.year  + 
                    s.summer.temp.cb + s.lag.summer.temp.cb +
                    s.summer.precip.cb + s.lag.summer.precip.cb +
                    s.snowpack + s.lag.snowpack + 
                    s.yr.floralcount + s.lag.yr.floralcount + 
                    s.year * s.summer.temp.cb + s.year * s.lag.summer.temp.cb +
                    s.year * s.summer.precip.cb + s.year * s.lag.summer.precip.cb +
                    s.year * s.lag.snowpack +
                    s.year * s.lag.yr.floralcount,
                  weights = sqrt(yr.samples), 
                  data = data.year.out, na.action = na.omit)

vif(q2.total.no) # next, removed s.floralcount:     3.100850     

q2.total.no <- lm(sum.resid.total.no ~ s.year  + 
                    s.summer.temp.cb + s.lag.summer.temp.cb +
                    s.summer.precip.cb + s.lag.summer.precip.cb +
                    s.snowpack + s.lag.snowpack + 
                    s.lag.yr.floralcount + 
                    s.year * s.summer.temp.cb + s.year * s.lag.summer.temp.cb +
                    s.year * s.summer.precip.cb + s.year * s.lag.summer.precip.cb +
                    s.year * s.lag.snowpack +
                    s.year * s.lag.yr.floralcount,
                  weights = sqrt(yr.samples), 
                  data = data.year.out, na.action = na.omit)
vif(q2.total.no) # next, removed s.lag.summer.temp.cb and year:lag.temp:     3.110661    

q2.total.no <- lm((sum.resid.total.no/yr.samples) ~ s.year  + 
                    s.summer.temp.cb + 
                    s.summer.precip.cb + s.lag.summer.precip.cb +
                    s.snowpack + s.lag.snowpack + 
                    s.lag.yr.floralcount + 
                    s.year * s.summer.temp.cb + 
                    s.year * s.summer.precip.cb + s.year * s.lag.summer.precip.cb +
                    s.year * s.lag.snowpack +
                    s.year * s.lag.yr.floralcount,
                  weights = sqrt(yr.samples), 
                  data = data.year.out, na.action = na.omit)
vif(q2.total.no) 

#summary(q2.total.no) 
summary(q2.selected.total.no <- step(q2.total.no, direction = "both", 
                                     scope = list(lower = ~s.year)))    
# year, year:precip, year:lag precip
extractAIC(q2.selected.total.no)

# Total insect biomass (global model, stepwise)
# removed same factors as abundance model, with VIF > 3
q2.total.g <- lm((sum.resid.total.g/yr.samples) ~ s.year +
                   s.summer.temp.cb + s.lag.summer.temp.cb +
                   s.summer.precip.cb + s.lag.summer.precip.cb +
                   s.snowpack + s.lag.snowpack + 
                   s.yr.floralcount + s.lag.yr.floralcount + 
                   s.year * s.summer.precip.cb + s.year * s.lag.summer.precip.cb +
                   s.year * s.snowpack + s.year * s.lag.snowpack,
                 weights = sqrt(yr.samples), na.action = na.omit, data = data.year.out)
vif(q2.total.g)  

summary(q2.selected.total.g <- step(q2.total.g, direction = "both", 
                                    scope = list(lower = ~s.year)))        # year, year:precip, marginal year:lag snow              
vif(q2.selected.total.g)
extractAIC(q2.selected.total.g)

