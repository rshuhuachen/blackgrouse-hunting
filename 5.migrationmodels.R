## In this script, we will first convert our dataframe to fit the format for BA3, then run BA3,
## analyze the 5 different runs for the migration models
## as constructed with BayeSass v3, and then we will model the effect of hunting
## on both emigration and immigration rates

#load packages
library(data.table); library(tidyverse); library(tibble); library(MuMIn)
library(lme4); library(lmerTest); library(readxl); library(DHARMa);library(glmmTMB); library(performance)

# load in pop data
#setwd("P:\\Black Grouse PhD\\Projects\\Hunting_Microsats\\Clean+final_analysis\\ScriptsHuntingFinal\\blackgrouse-hunting\\")
pops <- read.csv("data/Codes.pops.both.filtered_withcoord.csv")
pops$pop_num <- as.factor(pops$pop_num)

## load in Matrix distances between sites ##
distance <- read_excel("data/CalculateDistanceSitesGenAlEx.xlsx", sheet = "MatrixForR")
names(distance)[1] <- "Site_A"

distance_long <- melt(distance)
names(distance_long) <- c("Site_A", "Site_B", "Distance")
distance_long <- subset(distance_long, distance_long$Site_A != distance_long$Site_B)

####### Reformatting for BA3 #######

male.stru <- fread("data/Microsat.males.noLOCUS1+13.forstructure.stru")
head(male.stru)
names(male.stru) <- c("indivID", "popID", "BG16", "BG18", "BG15", "BG19", "BG6", "TTT1", "TTD2",
                      "TTD3", "TUD6", "TUT3", "TUT4", "TTT2")

male.stru.a <- male.stru[seq(from = 1, by = 2, to = nrow(male.stru)-1),]
male.stru.b <- male.stru[seq(from = 2, by = 2, to = nrow(male.stru)),]

ba3.male.a <- melt(data = male.stru.a,
              id.vars = c("indivID", "popID"),
              variable.name = "locID",
              value.name = "allele1")

ba3.male.b <- melt(data = male.stru.b,
              id.vars = c("indivID", "popID"),
              variable.name = "locID",
              value.name = "allele2")

ba3.male <- left_join(ba3.male.a, ba3.male.b, by = c("indivID", "popID", "locID"))

head(ba3.male)

write.table(ba3, "data/data_males_ba3.txt",
            col.names = T, row.names = F, sep = " ", quote = F)


#females
female.stru <- fread("data/Microsat.females.noLOCUS1+13.forstructure.stru")
head(female.stru)
names(female.stru) <- c("indivID", "popID", "BG16", "BG18", "BG15", "BG19", "BG6", "TTT1", "TTD2",
                      "TTD3", "TUD6", "TUT3", "TUT4", "TTT2")

female.stru.a <- female.stru[seq(from = 1, by = 2, to = nrow(female.stru)-1),]
female.stru.b <- female.stru[seq(from = 2, by = 2, to = nrow(female.stru)),]

ba3.female.a <- melt(data = female.stru.a,
                   id.vars = c("indivID", "popID"),
                   variable.name = "locID",
                   value.name = "allele1")

ba3.female.b <- melt(data = female.stru.b,
                   id.vars = c("indivID", "popID"),
                   variable.name = "locID",
                   value.name = "allele2")

ba3.female <- left_join(ba3.female.a, ba3.female.b, by = c("indivID", "popID", "locID"))

head(ba3.female)

write.table(ba3, "data/data_females_ba3.txt",
            col.names = T, row.names = F, sep = " ", quote = F)

#### Running BA3 ####

## males
#install BA3, run through command line/terminal

system("mkdir /data/BA3runs/males_run1") #make directory per run
system("mkdir /data/BA3runs/males_run2") #make directory per run
system("mkdir /data/BA3runs/males_run3") #make directory per run
system("mkdir /data/BA3runs/males_run4") #make directory per run
system("mkdir /data/BA3runs/males_run5") #make directory per run

#5 runs with 5 different random seeds
system("~/bin/BA3/BA3MSAT -v -t -g -u -a 0.30 -f 0.40 -s 65323 -i 10000000 -b 1000000 -n 1000 -o males_run1.txt /data/data_males_ba3.txt") #change pathname to BA3, random seed generated
system("~/bin/BA3/BA3MSAT -v -t -g -u -a 0.30 -f 0.40 -s 76553 -i 10000000 -b 1000000 -n 1000 -o males_run2.txt /data/data_males_ba3.txt") #change pathname to BA3, random seed generated
system("~/bin/BA3/BA3MSAT -v -t -g -u -a 0.30 -f 0.40 -s 124643 -i 10000000 -b 1000000 -n 1000 -o males_run3.txt /data/data_males_ba3.txt") #change pathname to BA3, random seed generated
system("~/bin/BA3/BA3MSAT -v -t -g -u -a 0.30 -f 0.40 -s 885256 -i 10000000 -b 1000000 -n 1000 -o males_run4txt /data/data_males_ba3.txt") #change pathname to BA3, random seed generated
system("~/bin/BA3/BA3MSAT -v -t -g -u -a 0.30 -f 0.40 -s 235776 -i 10000000 -b 1000000 -n 1000 -o males_run5.txt /data/data_males_ba3.txt") #change pathname to BA3, random seed generated

## females
#install BA3, run through command line/terminal

system("mkdir /data/BA3runs/females_run1") #make directory per run
system("mkdir /data/BA3runs/females_run2") #make directory per run
system("mkdir /data/BA3runs/females_run3") #make directory per run
system("mkdir /data/BA3runs/females_run4") #make directory per run
system("mkdir /data/BA3runs/females_run5") #make directory per run

#5 runs with 5 different random seeds
system("~/bin/BA3/BA3MSAT -v -t -g -u -a 0.30 -f 0.40 -s 65323 -i 10000000 -b 1000000 -n 1000 -o females_run1.txt /data/data_females_ba3.txt") #change pathname to BA3, random seed generated
system("~/bin/BA3/BA3MSAT -v -t -g -u -a 0.30 -f 0.40 -s 76553 -i 10000000 -b 1000000 -n 1000 -o females_run2.txt /data/data_females_ba3.txt") #change pathname to BA3, random seed generated
system("~/bin/BA3/BA3MSAT -v -t -g -u -a 0.30 -f 0.40 -s 124643 -i 10000000 -b 1000000 -n 1000 -o females_run3.txt /data/data_females_ba3.txt") #change pathname to BA3, random seed generated
system("~/bin/BA3/BA3MSAT -v -t -g -u -a 0.30 -f 0.40 -s 885256 -i 10000000 -b 1000000 -n 1000 -o females_run4txt /data/data_females_ba3.txt") #change pathname to BA3, random seed generated
system("~/bin/BA3/BA3MSAT -v -t -g -u -a 0.30 -f 0.40 -s 235776 -i 10000000 -b 1000000 -n 1000 -o females_run5.txt /data/data_females_ba3.txt") #change pathname to BA3, random seed generated

#### Compare all 10 runs ####
setwd("data/BA3runs")
temp <- list.files(pattern = "*.txt")
myfiles = lapply(temp, fread, skip = 18, nrows = 12, header = F)

# formula for reshaping the dataframes

reshape_ba3 <- function(m) {
  m1 <- m[,c(1,2)]
  names(m1) <- c("pops", "migration")
  m2 <- m[,c(3,4)]
  names(m2) <- c("pops", "migration")
  m3 <- m[,c(5,6)]
  names(m3) <- c("pops", "migration")
  m4 <- m[,c(7,8)]
  names(m4) <- c("pops", "migration")
  m5 <- m[,c(9,10)]
  names(m5) <- c("pops", "migration")
  m6 <- m[,c(11,12)]
  names(m6) <- c("pops", "migration")
  m7 <- m[,c(13,14)]
  names(m7) <- c("pops", "migration")
  m8 <- m[,c(15,16)]
  names(m8) <- c("pops", "migration")
  m9 <- m[,c(17,18)]
  names(m9) <- c("pops", "migration")
  m10 <- m[,c(19,20)]
  names(m10) <- c("pops", "migration")
  m11 <- m[,c(21,22)]
  names(m11) <- c("pops", "migration")
  m12 <- m[,c(23,24)]
  names(m12) <- c("pops", "migration")
  mnew<-rbind(m1, m2, m3, m4, m5,m6,m7,m8,m9,m10,m11,m12)
  # Note that m[i][j] is the fraction of individuals in population i that are migrants derived from population j per generation
  # Pops:  0->1 9->10 10->11 11->12 1->2 2->3 3->4 4->5 5->6 6->7 7->8 8->9
  mnew$m_in <- c(rep(c(1:12), times = 12, each = 1)) # rename populations into actual names and split into migration_in
  mnew$m_out <- c(rep(c(1:12), times = 1, each = 12)) # rename populations into actual names and split into migration_out
  mnew <- separate(data = mnew, col = "migration", into = c("migration", "migration_SE"), 
                   sep = "[(]") #seperate migration and its SE
  mnew$migration_SE <- gsub(mnew$migration_SE, pattern = "[)]", replacement = "") #take out parenthesis
  mnew <- mnew[,c(4,5,2,3)]
  return(mnew)
}

# run it for all files
for (i in 1:length(myfiles)) {
  myfiles[[i]]<-reshape_ba3(myfiles[[i]])
}

#separate for males and females
maleruns <- myfiles[c(6:10)]
femaleruns <- myfiles[c(1:5)]

#separate per run to compare
male_run1 <- maleruns[[1]]
male_run2 <- maleruns[[2]]
male_run3 <- maleruns[[3]]
male_run4 <- maleruns[[4]]
male_run5 <- maleruns[[5]]

female_run1 <- femaleruns[[1]]
female_run2 <- femaleruns[[2]]
female_run3 <- femaleruns[[3]]
female_run4 <- femaleruns[[4]]
female_run5 <- femaleruns[[5]]

#### Compare runs ####
plot(male_run1$migration, male_run2$migration)
plot(male_run1$migration, male_run3$migration)
plot(male_run1$migration, male_run4$migration)
plot(male_run1$migration, male_run5$migration)
plot(male_run2$migration, male_run3$migration)
plot(male_run2$migration, male_run4$migration)
plot(male_run2$migration, male_run5$migration)
plot(male_run3$migration, male_run4$migration)
plot(male_run3$migration, male_run5$migration)
plot(male_run4$migration, male_run5$migration)
# all runs correspond

plot(female_run1$migration, female_run2$migration)
plot(female_run1$migration, female_run3$migration)
plot(female_run1$migration, female_run4$migration)
plot(female_run1$migration, female_run5$migration)
plot(female_run2$migration, female_run3$migration)
plot(female_run2$migration, female_run4$migration)
plot(female_run2$migration, female_run5$migration)
plot(female_run3$migration, female_run4$migration)
plot(female_run3$migration, female_run5$migration)
plot(female_run4$migration, female_run5$migration)

# all runs correspond

## Going to pick run 5 for both

#### Clean migration files ####
# to do: add ESS, corrected migration value, hunted status, distance between sites

#first add ESS
ESS_male_run5 <- read.delim("tracer_summary/tracer_summary_males5.txt", sep = "\t")
ESS_female_run5 <- read.csv("tracer_summary/tracer_summary_females5.txt", sep = "\t")

ESS_male_run5.df <- t(ESS_male_run5)
colnames(ESS_male_run5.df) <- ESS_male_run5.df[1,]
ESS_male_run5.df <- as.data.frame(ESS_male_run5.df[-c(1,2),])#take out col names and row for log prob
ESS_male_run5.df <- rownames_to_column(ESS_male_run5.df, "m")

ESS_female_run5.df <- t(ESS_female_run5)
colnames(ESS_female_run5.df) <- ESS_female_run5.df[1,]
ESS_female_run5.df <- as.data.frame(ESS_female_run5.df[-c(1,2),])#take out col names and row for log prob
ESS_female_run5.df <- rownames_to_column(ESS_female_run5.df, "m")

#males
male_run5 <- male_run5 %>% arrange(m_in) #now, the two files are sorted the same way, so can just cbind
male_run5_clean <- cbind(male_run5, ESS_male_run5.df[,c(11)]) #select ESS

head(male_run5_clean)
names(male_run5_clean)[5] <- "ESS"
#females
female_run5 <- female_run5 %>% arrange(m_in) #now, the two files are sorted the same way, so can just cbind
female_run5_clean <- cbind(female_run5, ESS_female_run5.df[,c(11)]) #select ESS

head(female_run5_clean)
names(female_run5_clean)[5] <- "ESS"

#add migration_ESSc column, which turns to NA if ESS < 200
#male
male_run5_clean$ESS <- as.numeric(male_run5_clean$ESS)
male_run5_clean$migration <- as.numeric(male_run5_clean$migration)
male_run5_clean$migration_SE <- as.numeric(male_run5_clean$migration_SE)
male_run5_clean$migration_ESSc <- case_when(male_run5_clean$ESS < 200 ~ as.numeric(NA),
                                            male_run5_clean$ESS >= 200 ~ as.numeric(male_run5_clean$migration))

#female
female_run5_clean$ESS <- as.numeric(female_run5_clean$ESS)
female_run5_clean$migration <- as.numeric(female_run5_clean$migration)
female_run5_clean$migration_SE <- as.numeric(female_run5_clean$migration_SE)
female_run5_clean$migration_ESSc <- case_when(female_run5_clean$ESS < 200 ~ as.numeric(NA),
                                            female_run5_clean$ESS >= 200 ~ as.numeric(female_run5_clean$migration))


# then, add hunted status and full site names
male_run5_clean$m_in<-as.factor(male_run5_clean$m_in)
male_run5_clean$m_out<-as.factor(male_run5_clean$m_out)

male_run5_clean <- left_join(male_run5_clean, pops[,c(1:3)], by = c("m_in" = "pop_num"))
names(male_run5_clean)[7] <- "pop_in"
names(male_run5_clean)[8] <- "hunt_in"

male_run5_clean <- left_join(male_run5_clean, pops[,c(1:3)], by = c("m_out" = "pop_num"))
names(male_run5_clean)[9] <- "pop_out"
names(male_run5_clean)[10] <- "hunt_out"

female_run5_clean$m_in<-as.factor(female_run5_clean$m_in)
female_run5_clean$m_out<-as.factor(female_run5_clean$m_out)

female_run5_clean <- left_join(female_run5_clean, pops[,c(1:3)], by = c("m_in" = "pop_num"))
names(female_run5_clean)[7] <- "pop_in"
names(female_run5_clean)[8] <- "hunt_in"

female_run5_clean <- left_join(female_run5_clean, pops[,c(1:3)], by = c("m_out" = "pop_num"))
names(female_run5_clean)[9] <- "pop_out"
names(female_run5_clean)[10] <- "hunt_out"

## add distance
male_run5_clean <- left_join(male_run5_clean, distance_long, by = c("pop_in" = "Site_A", "pop_out" = "Site_B"))
female_run5_clean <- left_join(female_run5_clean, distance_long, by = c("pop_in" = "Site_A", "pop_out" = "Site_B"))

write.table(male_run5_clean, file = "data/run5_males_clean.csv", sep = ",", row.names = F,quote = F)
write.table(female_run5_clean, file = "data/run5_females_clean.csv", sep = ",", row.names = F,quote = F)

#### Load in clean runs ####

#this file has been cleaned up, hunted status added, distances between sites and corrected migration value excluding those with ESS < 200

male_run5_clean <- read.csv("data/run5_males_clean.csv")
female_run5_clean <- read.csv("data/run5_females_clean.csv")

### Plotting migration rates from run 5 ####
#first, exclude the 'non-migration rates' which are those where pop in = pop out
male_run5_clean <- subset(male_run5_clean, m_in != m_out)
female_run5_clean <- subset(female_run5_clean, m_in != m_out)

theme_set(theme_classic())
# males: out hunted vs out unhunted but exclude nonmigration rates

ttest_males_in_rates <- t.test(male_run5_clean$migration_ESSc[which(male_run5_clean$hunt_in == "hunted" & male_run5_clean$m_in != male_run5_clean$m_out)], male_run5_clean$migration_ESSc[which(male_run5_clean$hunt_in == "unhunted" & male_run5_clean$m_in != male_run5_clean$m_out)])

ggplot(male_run5_clean, aes(x = hunt_in, y = migration_ESSc)) + geom_boxplot(aes(fill = "Type")) +
  labs(title = "Migration rates IN for males") + xlab("Hunted status") + ylab("migration rates in")+
  labs(subtitle = 
         paste("Excluding values with ESS < 200 
Student t-test comparing hunted/unhunted site filtered migration rate in: 
t = ", round(ttest_males_in_rates$statistic,2), "and p-value = ", round(ttest_males_in_rates$p.value, 2)))+ 
  theme(text = element_text(family = "Arial", size = 22),
        plot.subtitle = element_text(size = 14),
        legend.position = "none")+
  scale_fill_manual(values = c("cyan3"))

# males: in hunted vs in unhunted but exclude nonmigration_ESSc rates

ttest_males_out_rates <- t.test(male_run5_clean$migration_ESSc[which(male_run5_clean$hunt_out == "hunted")], male_run5_clean$migration_ESSc[which(male_run5_clean$hunt_out == "unhunted")])

ggplot(male_run5_clean, aes(x = hunt_out, y = migration_ESSc, fill = "blue")) + geom_boxplot(aes(fill = "Type")) +
  labs(title = "Migration rates OUT for males") + xlab("Hunted status") + ylab("migration rates out")+
  labs(subtitle = 
         paste("Excluding values with ESS < 200 
Student t-test comparing hunted/unhunted site corrected migration rate out: 
t = ", round(ttest_males_out_rates$statistic,2), "and p-value = ", round(ttest_males_out_rates$p.value, 3)))+ 
  theme(text = element_text(family = "Arial", size = 22),
        plot.subtitle = element_text(size = 14),
        legend.position = "none")+
  scale_fill_manual(values = c("cyan3"))

# females: out hunted vs out unhunted but exclude nonmigration_ESSc rates
ttest_females_in_rates <- t.test(female_run5_clean$migration_ESSc[which(female_run5_clean$hunt_in == "hunted" & female_run5_clean$m_in != female_run5_clean$m_out)], female_run5_clean$migration_ESSc[which(female_run5_clean$hunt_in == "unhunted" & female_run5_clean$m_in != female_run5_clean$m_out)])

ggplot(female_run5_clean, aes(x = hunt_in, y = migration_ESSc, fill = "red")) + geom_boxplot() +
  labs(title = "Migration rates IN for females") + xlab("Hunted status") + ylab("migration rates in")+
  labs(subtitle = 
         paste("Excluding values with ESS < 200 
Student t-test comparing hunted/unhunted site corrected migration rate in: 
t = ", round(ttest_females_in_rates$statistic,2), "and p-value = ", round(ttest_females_in_rates$p.value, 2)))+ 
  theme(text = element_text(family = "Arial", size = 22),
        plot.subtitle = element_text(size = 14),
        legend.position = "none")

# females: in hunted vs in unhunted but exclude nonmigration_ESSc rates

ttest_females_out_rates <- t.test(female_run5_clean$migration_ESSc[which(female_run5_clean$hunt_out == "hunted" & female_run5_clean$m_in != female_run5_clean$m_out)], female_run5_clean$migration_ESSc[which(female_run5_clean$hunt_out == "unhunted" & female_run5_clean$m_in != female_run5_clean$m_out)])

ggplot(female_run5_clean, aes(x = hunt_out, y = migration_ESSc, fill = "red")) + geom_boxplot() +
  labs(title = "Migration rates OUT for females") + xlab("Hunted status") + ylab("migration rates out")+
  labs(subtitle = 
         paste("Excluding values with ESS < 200
Student t-test comparoutg hunted/unhunted site corrected migration rate out: 
t = ", round(ttest_females_out_rates$statistic,2), "and p-value = ", round(ttest_females_out_rates$p.value, 2)))+ 
  theme(text = element_text(family = "Arial", size = 22),
        plot.subtitle = element_text(size = 14),
        legend.position = "none")

#### Modelling migration ####

names(female_run5_clean) == names(male_run5_clean)
female_run5_clean$sex <- "Female"
male_run5_clean$sex <- "Male"

#change levels hunted/unhunted
male_run5_clean$hunt_in <- relevel(as.factor(male_run5_clean$hunt_in), ref = "unhunted")
male_run5_clean$hunt_out <- relevel(as.factor(male_run5_clean$hunt_out), ref = "unhunted")

female_run5_clean$hunt_in <- relevel(as.factor(female_run5_clean$hunt_in), ref = "unhunted")
female_run5_clean$hunt_out <- relevel(as.factor(female_run5_clean$hunt_out), ref = "unhunted")

#combine in one df
migration_both <- rbind(male_run5_clean, female_run5_clean)

#immigration

model.both.in <- glmmTMB(migration_ESSc~ hunt_in + Distance + sex + (1|pop_out) + (1|pop_in), 
                         data = migration_both, 
                         family = Gamma(link = "log"))
summary(model.both.in) 
plot(model.both.in)
simulateResiduals(fittedModel = model.both.in, plot = T)

model.both.in.interaction <- glmmTMB(migration_ESSc~ hunt_in + Distance*sex + (1|pop_out) + (1|pop_in), 
                         data = migration_both, 
                         family = Gamma(link = "log"))
summary(model.both.in.interaction) 
plot(model.both.in.interaction)
simulateResiduals(fittedModel = model.both.in.interaction, plot = T)

#out
model.both.out <- glmmTMB(migration_ESSc~ hunt_out + Distance + sex + (1|pop_out) + (1|pop_in), 
                          data = migration_both, 
                          family = Gamma(link = "log"))
summary(model.both.out) 
plot(model.both.out)
simulateResiduals(fittedModel = model.both.out, plot = T)

model.both.out.interaction <- glmmTMB(migration_ESSc~ hunt_out + Distance*sex + (1|pop_out) + (1|pop_in), 
                          data = migration_both, 
                          family = Gamma(link = "log"))
summary(model.both.out.interaction) 
plot(model.both.out.interaction)
simulateResiduals(fittedModel = model.both.out.interaction, plot = T)

## model performance
icc(model = model.both.in, by_group = TRUE)
icc(model = model.both.in.interaction, by_group = TRUE)
icc(model = model.both.out, by_group = TRUE)
icc(model = model.both.out.interaction, by_group = TRUE)


r.squaredGLMM(model.both.in)
r.squaredGLMM(model.both.in.interaction)
r.squaredGLMM(model.both.out)
r.squaredGLMM(model.both.out.interaction)

compare_performance(model.both.in, model.both.in.interaction, rank=T)
compare_performance(model.both.out, model.both.out.interaction, rank=T)
