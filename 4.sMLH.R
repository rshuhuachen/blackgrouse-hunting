### Effect of hunting on genetic diversity ###

#load packages
### Load packages ###
library(dplyr);library(tibble);library(ggplot2)
library(inbreedR); library(data.table)
library(lmerTest); library(DHARMa); library(performance);library(MuMIn)

#setwd("/Volumes/rchen2/Black Grouse PhD/Projects/Hunting_Microsats/Clean+final_analysis/ScriptsHuntingFinal/GithubDirectory") #set to directory of github repository

sitenames <- read.csv("data/details/Codes.pops.both.filtered_withcoord.csv") #information for each site
sitenames[1,1] <- "Koskenpää"
sitenames[5,1] <- "Nyrölä"

#density
dens <- read.csv("data/")


#Load in filtered data
males <- fread("data/cleandata/Unsplit.microsat.males.noLOCUS1+13.csv") #different format than structure files
males <- as.data.frame(males)
males <- left_join(males, sitenames[,c(1,2)], by = c("pop" = "pop_num")) #join to have full pop names
males <- males[,c(1,27,3:26)] #reorganise
names(males)[2] <- "pop"

females <- fread("data/cleandata/Unsplit.microsat.females.noLOCUS1+13.csv") 
females <- left_join(females, sitenames[,c(1,2)], by = c("pop" = "pop_num"))
females <- females[,c(1,27,3:26)] #reorganise
names(females)[2] <- "pop"

chicks <- fread("data/cleandata/Unsplit.microsat.chicks.noLOCUS1+13+14.csv")
chicks <- as.data.frame(chicks)
chicks <- left_join(chicks, sitenames[,c(1,2)], by = c("pop" = "pop_num"))
chicks <- chicks[,c(1,25,3:24)]
names(chicks)[2] <- "pop"

## Change formats to be loaded into inbreedR (binary where rownames are id's and no populations)
males.inb <- males
males.inb <- males.inb %>% remove_rownames %>% column_to_rownames(var="id")
males.inb <- males.inb[,-1]

females.inb <- females
females.inb <- females.inb %>% remove_rownames %>% column_to_rownames(var="id")
females.inb <- females.inb[,-1]

chicks.inb <- chicks
chicks.inb <- chicks.inb %>% remove_rownames %>% column_to_rownames(var="id")
chicks.inb <- chicks.inb[,-1]

# convert to inbreedR
males.inb <- convert_raw(males.inb)
females.inb <- convert_raw(females.inb)
chicks.inb <- convert_raw(chicks.inb)

#### Calculate sMLH ####

sMLH_females <- sMLH(females.inb) #sMLH
het_var_females <- var(sMLH_females, na.rm=TRUE) # variance in sMLH
summary(sMLH_females)
het_var_females

sMLH_males <- sMLH(males.inb) #sMLH
het_var_males <- var(sMLH_males, na.rm=TRUE) # variance in sMLH
summary(sMLH_males)
het_var_males

sMLH_chicks <- sMLH(chicks.inb)
het_var_chicks <- var(sMLH_chicks, na.rm=TRUE) # variance in sMLH
summary(sMLH_chicks)
het_var_chicks

# Add all to dataframe 

sMLH_males.df <- as.data.frame(sMLH_males)
colnames(sMLH_males.df)[1] <- "sMLH"
sMLH_males.df$age <- "adult"
sMLH_males.df <- tibble::rownames_to_column(sMLH_males.df, "id")
sMLH_males.df <- left_join(sMLH_males.df, males[,c(1,2)], by = "id")
sMLH_males.df$sex <- "male"

sMLH_females.df <- as.data.frame(sMLH_females)
colnames(sMLH_females.df)[1] <- "sMLH"
sMLH_females.df$age <- "adult"
sMLH_females.df <- tibble::rownames_to_column(sMLH_females.df, "id")
sMLH_females.df <- left_join(sMLH_females.df, females[,c(1,2)], by = "id")
sMLH_females.df$sex <- "female"

# add together with hunted status
load("data/rawdata/Fulldata_adults.RData")

sMLH_males.df <- left_join(sMLH_males.df, hunted.ad[,c(1,3,8)], by = c("id")) #select id, year, hunt
sMLH_males.df <- sMLH_males.df[,c(1,3,5,4,6,7,2)] #rearrange

sMLH_females.df <- left_join(sMLH_females.df, hunted.ad[,c(1,3,8)], by = c("id")) #select id, year, hunt
sMLH_females.df <- sMLH_females.df[,c(1,3,5,4,6,7,2)] #rearrange

# Chicks

sMLH_chicks.df <- as.data.frame(sMLH_chicks)
colnames(sMLH_chicks.df)[1] <- "sMLH"
sMLH_chicks.df$age <- "chick"
sMLH_chicks.df <- tibble::rownames_to_column(sMLH_chicks.df, "id")
sMLH_chicks.df <- left_join(sMLH_chicks.df, chicks[,c(1,2)], by = "id")

#add with hunted status
load("data/rawdata/Fulldata_chicks.RData")

sMLH_chicks.df <- left_join(sMLH_chicks.df, hunted.chick[,c(9,2,13,3)], by = c("id" = "ID")) #select sex, year, id,hunt
sMLH_chicks.df <- sMLH_chicks.df[,c(1,3,5,4,6,7,2)]
names(sMLH_males.df)
names(sMLH_females.df)
names(sMLH_chicks.df)[5] <- "year"

sMLH.all <- rbind(sMLH_males.df, sMLH_females.df, sMLH_chicks.df)

#ensure all the same:
levels(as.factor(sMLH.all$sex))
sMLH.all$sex[which(sMLH.all$sex == "F")] <- "f"
sMLH.all$sex[which(sMLH.all$sex == "M")] <- "m"
sMLH.all$sex[which(sMLH.all$sex == "NA")] <- NA
sMLH.all$sex[which(sMLH.all$sex == "female")] <- "f"
sMLH.all$sex[which(sMLH.all$sex == "male")] <- "m"

#add hunted/unhunted
msatsfactors <- c(1:6)
sMLH.all[msatsfactors] <- lapply(sMLH.all[factors], as.factor)

## Plotting
sMLH.all %>% filter(!is.na(sex)) %>%ggplot(aes(x = age, y = sMLH, col = sex)) + 
  geom_boxplot() + theme_classic ()

#### Or: load in sMLH ####
#save(sMLH.all, file = "data/cleandata/sMLH.all.RData")
load(file = "data/cleandata/sMLH.all.RData")

#### Modelling sMLH ####

# first divide up again in female adults, male adults and chicks
sMLH.males <- subset(sMLH.all, sex == "m" & age == "adult")
sMLH.females <- subset(sMLH.all, sex == "f" & age == "adult")
sMLH.chicks <- subset(sMLH.all, age == "chick")
sMLH.adults <- subset(sMLH.all, age != "chick")

## Look at distribution sMLH ##

ggplot(sMLH.males, aes(x=sMLH)) + geom_histogram() +theme_classic() +
  labs(title = "Histogram of sMLH for adult males")

ggplot(sMLH.females, aes(x=sMLH)) + geom_histogram() +theme_classic() +
  labs(title = "Histogram of sMLH for adult females")

ggplot(sMLH.chicks, aes(x=sMLH)) + geom_histogram() +theme_classic() +
  labs(title = "Histogram of sMLH for chicks")

# not perfectly normal, but good enough. also tried the below models with 
# both a poisson and gamma distribution, which had a bad fit

#male model
sMLH.model.males.lmer <- lmerTest::lmer(sMLH ~ hunt + pop*year + hunt*year + (1|pop)+ (1|year) , data = sMLH.males)

#failed to converge
#Rescale and center continuous parameters
numcols <- grep("^c\\.",names(sMLH.males))
sMLH.males_rescale <- sMLH.males
sMLH.males_rescale[,numcols] <- scale(sMLH.males_rescale[,numcols])
sMLH.model.males.lmer_rescale <- update(sMLH.model.males.lmer,data=sMLH.males_rescale)

# restart and bump up max interations
ss <- getME(sMLH.model.males.lmer_rescale,c("theta","fixef"))
sMLH.model.males.lmer_noerror <- update(sMLH.model.males.lmer_rescale,start=ss,control=lmerControl(optCtrl=list(maxfun=2e4)))

#fixed
coef(summary(sMLH.model.males.lmer_noerror))
VarCorr(sMLH.model.males.lmer_noerror)
simulateResiduals(fittedModel = sMLH.model.males.lmer_noerror, plot = T)
plot(sMLH.model.males.lmer_noerror)
r.squaredGLMM(sMLH.model.males.lmer_noerror)
icc(model = sMLH.model.males.lmer_noerror, by_group = TRUE)

#females
sMLH.model.females.lmer <- lmerTest::lmer(sMLH ~ hunt + pop*year + hunt*year+(1|pop)+ (1|year) , data = sMLH.females)
coef(summary(sMLH.model.females.lmer))
VarCorr(sMLH.model.females.lmer)
simulateResiduals(fittedModel = sMLH.model.females.lmer, plot = T)
plot(sMLH.model.females.lmer)
r.squaredGLMM(sMLH.model.females.lmer)
icc(model = sMLH.model.females.lmer, by_group = TRUE)

## chicks
sMLH.model.chicks.lmer <- lmerTest::lmer(sMLH ~ hunt + pop*year + hunt*year+
                                           (1|pop)+ (1|year) , data = sMLH.chicks)

# model failed to converge, still fails when taking out pop*year or hunt*year

#Rescale and center continuous parameters
numcols <- grep("^c\\.",names(sMLH.chicks))
sMLH.chicks_rescale <- sMLH.chicks
sMLH.chicks_rescale[,numcols] <- scale(sMLH.chicks_rescale[,numcols])
sMLH.model.chicks.lmer_rescale <- update(sMLH.model.chicks.lmer,data=sMLH.chicks_rescale)

# restart and bump up max interations
ss <- getME(sMLH.model.chicks.lmer_rescale,c("theta","fixef"))
sMLH.model.chicks.lmer_noerror <- update(sMLH.model.chicks.lmer_rescale,start=ss,control=lmerControl(optCtrl=list(maxfun=2e4)))

# fixed
coef(summary(sMLH.model.chicks.lmer_noerror))
VarCorr(sMLH.model.chicks.lmer_noerror)
simulateResiduals(fittedModel = sMLH.model.chicks.lmer_noerror, plot = T)
plot(sMLH.model.chicks.lmer_noerror)
r.squaredGLMM(sMLH.model.chicks.lmer_noerror)
icc(model = sMLH.model.chicks.lmer_noerror, by_group = TRUE)
