### Effect of hunting on genetic diversity ###

### Load packages ###
pacman::p_load(dplyr, tibble, ggplot2, inbreedR, data.table, lmerTest,
               DHARMa, performance, MuMIn)

sitenames <- read.csv("data/details/Codes.pops.both.filtered_withcoord.csv") #information for each site
sitenames[1,1] <- "Koskenpää"
sitenames[5,1] <- "Nyrölä"

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

chicks <- fread("data/cleandata/Unsplit.microsat.unrelated.chicks.noLOCUS1+13+14.csv")
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
sMLH.all[msatsfactors] <- lapply(sMLH.all[msatsfactors], as.factor)

## Plotting
sMLH.all %>% filter(!is.na(sex)) %>%ggplot(aes(x = age, y = sMLH, col = sex)) + 
  geom_boxplot() + theme_classic ()

#### Modelling sMLH ####
#load in density data
dens <- read.csv("data/details/DensityLekYear.csv")
dens$site[which(dens$site == "Koskenp\xe4\xe4")] <- "Koskenpää"
dens$site[which(dens$site == "Nyr\xf6l\xe4")] <- "Nyrölä"
dens$year <- as.factor(dens$year)

sMLH.all <- left_join(sMLH.all, dens[,c(1,2,4)], by = c("pop" = "site", "year" = "year"))
sMLH.all <- unique(sMLH.all)

sMLH.all$hunt <- relevel(as.factor(sMLH.all$hunt), ref = "unhunted")

#### Or: load in sMLH ####
save(sMLH.all, file = "data/cleandata/sMLH.all.RData")
load(file = "data/cleandata/sMLH.all.RData")

#### Modelling sMLH with density ####

## building models
sMLH.all.no.na <- subset(sMLH.all, !is.na(density))
sMLH.model.all.lmer <- lmerTest::lmer(sMLH ~ hunt*density + sex + age+ (1|pop), data = sMLH.all.no.na)
sMLH.model.all.lmer.noint <- lmerTest::lmer(sMLH ~ hunt + sex + age+ density + (1|pop), data = sMLH.all.no.na)
sMLH.model.all.lmer.null <- lmerTest::lmer(sMLH ~ density+ sex + age+ (1|pop), data = sMLH.all.no.na)

anova(sMLH.model.all.lmer.null, sMLH.model.all.lmer)
anova(sMLH.model.all.lmer.null, sMLH.model.all.lmer.noint)

coef(summary(sMLH.model.all.lmer))
VarCorr(sMLH.model.all.lmer)
simulateResiduals(fittedModel = sMLH.model.all.lmer, plot = T)
plot(sMLH.model.all.lmer)
r.squaredGLMM(sMLH.model.all.lmer)
icc(model = sMLH.model.all.lmer, by_group = TRUE)

compare_performance(sMLH.model.all.lmer.null, sMLH.model.all.lmer, rank = T)
compare_performance(sMLH.model.all.lmer.null, sMLH.model.all.lmer.noint, rank = T)
