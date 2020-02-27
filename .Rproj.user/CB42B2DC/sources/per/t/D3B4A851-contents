#Final Project - Christian Carson - Data Analysis in R
########Sea Lice Resistance - Mized Effects Models Template#########

# Date: 2/19/2020
# Last updated: 
# Name: Christian Carson
# Description:
  # Index:
#18-100: Setup and Workflow
#100-138: Analysis Setup
#138-219: Mixed Models for 2019 Resistance Report
#219-316: Dose response curves for 2019 Resistance Model and EC50 estimates
#
#
#
###########Setup for 2019 Resistance################
### INPUT ###
rm(list=ls())
graphics.off()

wd <- getwd()  # working directory

#create folders within wd, make sure to take off # if this is your first time loading in the data and starting the analysis
#dir.create("Analysis")
#dir.create("Plots")
#dir.create("Analysis/Data-generated")
#automated vers
folders <- c("Analysis","Plots", "Analysis/Data-generated", "Data-Raw")
# function to create folders below

for(i in 1:length(folders)){
  if(file.exists(folders[i]) == FALSE)
    dir.create(folders[i])}

# we also need to store the paths to all of these folders
path.analysis <- paste(wd, "/", folders[1], sep = "")
path.plots <- paste(wd, "/", folders[2], sep = "")
path.data.generated <- paste(wd, "/", folders[3], sep = "")
data.path <- paste(wd, "/", folders[4], sep = "")


### INITIAL SET-UP ###

#raw data folder will alread be in the zip with the data

#data import
SLICE.2019 <- read.csv(paste(data.path, "/", "SLICEResistanceCedarCoast2019.csv", sep = ""), stringsAsFactors = FALSE)


#remove unessisary variables
names(SLICE.2019)

remove.columns <- c("location","temp_0_hrs","temp_12_hrs","avg_temp_0_6_12_18_24hrs","comments","salmon_lifestage","temp_1_hr","temp_18_hrs","moribund","region","Trial_start","temp_6_hrs","temp_24_hrs","Ethanol_vial.")

# we write a for loop to remove these columns from the data frame
 for(i in remove.columns){
  SLICE.2019 <- SLICE.2019[, !(colnames(SLICE.2019) %in% i)]
}

View(SLICE.2019)
#structure data files
str(SLICE.2019)
head(SLICE.2019)
names(SLICE.2019)
summary(SLICE.2019)

#we want to ensure all the variable are in the right format using as. functions to create a factored list for each variables values. This will be helpful for things like adding dates as a random effects in the mixed effects models
#create collection date as.factor
SLICE.2019$date_of_collection <-  as.factor(SLICE.2019$date_of_collection)
#SLICE conc as.integer
SLICE.2019$slice_conc_PPB <- as.integer(SLICE.2019$slice_conc_PPB)
#dish as.factor
SLICE.2019$dish <-  as.factor(SLICE.2019$dish)
#sex as.factor
SLICE.2019$louse_sex <- as.factor(SLICE.2019$louse_sex)
#stage as.factor
SLICE.2019$louse_stage <- as.factor(SLICE.2019$louse_stage)
# alive as.factor
SLICE.2019$live <-  as.factor(SLICE.2019$live)
# dead as.factor
SLICE.2019$dead <-  as.factor(SLICE.2019$dead)


### LOAD PACKAGES ###
require(rlang)
require(ggplot2)
require(GGally)
require(reshape2)
require(lme4)
require(compiler)
require(parallel)
require(boot)
require(lattice)
require(dplyr)
#update.packages()


######Begin Analysis 2019 Resistance#####                                                                             
names(SLICE.2019)
#using dead_and_moribund and slice_conc_PPB to see the influence of various conc
#-entrations of slice on the probability of sea lice dieing

SLICE.M <- SLICE.2019[SLICE.2019$louse_sex=='M',]
SLICE.F <- SLICE.2019[SLICE.2019$louse_sex=='F',]

# vector to be used for changing the order of life stages from alphabetical to order of develpment
level_order <- c('P1', 'P2', 'A')

#visualize data
#female lice across stages
tmp <- melt(SLICE.F[, c("louse_stage", "slice_conc_PPB")], id.vars="louse_stage")
Plot.Females <- ggplot(tmp, aes(x = louse_stage, y = value)) +
  geom_jitter(alpha = .5) +
  ggtitle('Female Sea Lice')+
  ylab('EMB concentration')+
  xlab('Life stage')+
  facet_grid(variable ~ .) +
  theme_classic()
ggsave("FemaleSealiceNestedness",path = path.plots, plot = Plot.Females, device = "png")

#male lice across stages
tmp <- melt(SLICE.M[, c("louse_stage", "slice_conc_PPB")], id.vars="louse_stage")
Plot.Males <- ggplot(tmp, aes(x = louse_stage, y = value)) +
  geom_jitter(alpha = .5)  +
  ggtitle('Male Sea Lice')+
  ylab('EMB concentration')+
  xlab('Life stage')+
  facet_grid(variable ~ .) +
  theme_classic()
ggsave("MaleSealiceNestedness",path = path.plots, plot = Plot.Males, device = "png")

#so we know we are going to use a glm because our data is not normal and we are
#dealing with a binomial distribution. Lets start with our first model, with no
#fixed or random effects right now

######Models for 2019 Resistance######
#setup models -> our variables are slice_conc_PPB and dead_and_moribund

m.1 <- glm(dead_and_moribund ~ slice_conc_PPB, data = SLICE.2019, 
          family = binomial) 
summary(m.1)
logLik(m.1)

#there are several issues with this model, notably, it does not account for the
#fact there are multiple stages of lice which are known to be affected by 
#SLICE diffrenetly, but also, the sex lice also influences their resistance
#m2 should include both of these as explanitory variables

m.2 <- glm(dead_and_moribund ~ slice_conc_PPB + louse_sex + louse_stage, data = SLICE.2019, family = binomial) 
summary(m.2)
logLik(m.2)
#ok lets compare to our first model using ANOVA
anova(m.1,m.2)
#with AIC
AIC(m.1,m.2)

#Ok, so just to make sure m2 is the best model, lets do a step wise check to
#compare other models with variations in the fixed effects
#just accounting for sex difference
#m.3 (Concentration + Louse Sex)
m.3 <- glm(dead_and_moribund ~  slice_conc_PPB + louse_sex, data = SLICE.2019, family = binomial)
#just accounting for stage difference
# m.4 (Concentration + Louse Stage)
m.4 <- glm(dead_and_moribund ~  slice_conc_PPB + louse_stage, data = SLICE.2019, family = binomial)
#removing concentration all together and comparing to stage and sex
#m.5 (Louse Sex +Louse Stage )
m.5 <- glm(dead_and_moribund ~  louse_sex + louse_stage, data = SLICE.2019, family = binomial)
#vs just sex
#m.6 (Louse Sex)
m.6 <- glm(dead_and_moribund ~  louse_sex, data = SLICE.2019, family = binomial)
#vs just stage
#m.7 (Louse Stage)
m.7 <- glm(dead_and_moribund ~  louse_stage, data = SLICE.2019, family = binomial)

#lets compare all of these with an ANOVA and AIC
anova(m.1,m.2,m.3,m.4,m.5,m.6,m.7)
#AIC
AIC(m.1,m.2,m.3,m.4,m.5,m.6,m.7)
logLik(m.3)
logLik(m.4)
logLik(m.5)
logLik(m.6)
logLik(m.7)

#its clear that out of all of these fixed effects models, m.2 has the best fit to the data

#Now lets take a crack at implementing mixed effects, which include random and fixed effects. Our first one includes the fixed effects:SLICE_conc_PPB, louse_sex, louse_stage. Random effects (random intercepts): for now just date of collection to adjust for the vartiation sea lice per dish
?glmer
logLik(m.3)#note that we are now using glmer of lme4 because we are implemting mixed effects now. Look for date name using names()
#trying random effect of dish
m.8 <- glmer(dead_and_moribund ~  slice_conc_PPB + louse_sex + (1 | dish), data = SLICE.2019, family = binomial)
#summary and model comparison
summary(m.8)
anova(m.1,m.2,m.3,m.4,m.5,m.6,m.7, m.8)
#notice how anova did not include either of the random effects models 
AIC(m.1,m.2,m.3,m.4,m.5,m.6,m.7, m.8)
#I decided to got with model 2 because it best explained my data with the least amount of parameters. Though m.9 had a lower AIC, it was not significantly lower (2 values)

#lets plot m2 as a logistic regression, we know this is a logistic regression 
#because of the binomial nature of the response variable
m.2
summary(m.2)
#first, we need to get the intercepts of each of the parameters used so we can plot slopes later on
# assigining coefficients and random variance
intercept <- summary(m.2)$coefficients['(Intercept)', 1]
slopeconc <- summary(m.2)$coefficients['slice_conc_PPB', 1]
slopemale<- summary(m.2)$coefficients['louse_sexM', 1]
slopep1<- summary(m.2)$coefficients['louse_stageP1', 1]
slopep2<- summary(m.2)$coefficients['louse_stageP2', 1]
#print them all to view
intercept
slopeconc
slopemale
slopep1
slopep2

#######Graphing for 2019 Sensitivity Estimates - Dose Response Curves and EC50's#####
#setting up base plots
library(ggplot2)
male.base <- ggplot(data = SLICE.2019, aes(x=SLICE.2019$slice_conc_PPB, y=SLICE.2019$dead_and_moribund))  +
  #add response and explanatory base, fixed effects added later
  geom_point()+ 
  theme_bw()+
  #adding asthetics and titles
  labs(x= "EMB concentration (ppb)",y= expression("Probability of death"),title='Lethal Dose-Response Curve for Male Lice')+
  #adjust y lim so its all out of 1
  ylim(0,1)+
  #add .abline at .5 for EC50 visual = aka the concentration needed to get 50% of
  #the maximal response = death
  geom_hline(yintercept=0.5, size=1)+
  #no gridlines
  theme_classic()

male.base
##plotting males, we are just gonna do this by hand basically. We want to plot each fixed effect intercept as a slope using stat_function
?stat_function

stat_function

#start with assigning the male base and then begin adding each slope
#add pal 1's, the equation is provided within the function after stetting the functions dataframe.
#backround: https://stats.idre.ucla.edu/other/mult-pkg/introduction-to-generalized-linear-mixed-models/
#equation for slope based on:https://stats.idre.ucla.edu/stata/seminars/deciphering-interactions-in-logistic-regression/
#names for ref:
intercept
slopeconc
slopemale
slopep1
slopep2
summary(m.2)

males.prob <- male.base +
  stat_function(aes(color = "Adult"),size=1.3,
#equation:p = 1/1+e-(b0 +b1x+b2x+b3x)
#Adult
                fun = function(x){1/(1+exp(-intercept-slopeconc*(x)-slopemale*(1)-slopep1*(0)-slopep2*(0)))})+
#male p1
  stat_function(aes(color = "PAL I"),size=1.3,fun = function(x){1/(1+exp(-intercept-slopeconc*(x)-slopemale*(1)-slopep1*(1)-slopep2*(0)))})+
#male p2
stat_function(aes(color = "PAL II"),size=1.3,fun = function(x){1/(1+exp(-intercept-slopeconc*(x)-slopemale*(1)-slopep1*(0)-slopep2*(1)))})


males.prob
ggsave("DoseResponseCurveMaleLice",path = path.plots, plot = males.prob, device = "png")

#females
female.base <- ggplot(data = SLICE.2019, aes(x=SLICE.2019$slice_conc_PPB, y=SLICE.2019$dead_and_moribund))  +
  geom_point()+ 
  theme_bw()+
  labs(x= "EMB concentration (ppb)",y= expression("Probability of death"),title='Lethal Dose-Response Curve for Female Lice')+
  ylim(0,1)+
  geom_hline(yintercept=0.5, size=1)+
  theme_classic()
female.prob <- female.base + #equation:p = 1/1+e-(b0 +b1x+b2x+b3x)
          
  stat_function(aes(color = "PAL I"),size=1.3,fun = function(x){1/(1+exp(-intercept-slopeconc*(x)-slopemale*(0)-slopep1*(1)-slopep2*(0)))})+
  #male p2
  stat_function(aes(color = "PAL 2"),size=1.3,fun = function(x){1/(1+exp(-intercept-slopeconc*(x)-slopemale*(0)-slopep1*(0)-slopep2*(1)))})


female.prob
ggsave("DoseResponseCurveFemaleLice",path = path.plots, plot = female.prob, device = "png")

?ggsave
#EC50 estimates
#EC50 values
#assign 50% probaility 

p <- 0.5
#EC50 female pal 1
EC50.FP1 <- (log(p/(1-p)) - intercept-slopemale*(0)-slopep1*(1)-slopep2*(0)) / (slopeconc)
#EC50 female pal 2
EC50.FP2 <- (log(p/(1-p)) - intercept-slopemale*(0)-slopep1*(0)-slopep2*(1)) / (slopeconc)
#EC50 male pal 1
EC50.MP1 <- (log(p/(1-p)) - intercept-slopemale*(1)-slopep1*(1)-slopep2*(0)) / (slopeconc)
#EC50 male pal 2
EC50.MP2 <- (log(p/(1-p)) - intercept-slopemale*(1)-slopep1*(0)-slopep2*(1)) / (slopeconc)
#EC50 male adult
EC50.MADULT <- (log(p/(1-p)) - intercept-slopemale*(1)-slopep1*(0)-slopep2*(0)) / (slopeconc)

# EC50 VALUES
#Female Pal 1
EC50.FP1
#Female Pal 2
EC50.FP2
#Male Pal 1
EC50.MP1
#Male Pal 2
EC50.MP2
#Male Adult
EC50.MADULT

