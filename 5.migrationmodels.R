## In this script, we will first convert our dataframe to fit the format for BA3, then run BA3,
## analyze the 5 different runs for the migration models
## as constructed with BayeSass v3, and then we will model the effect of hunting
## on both emigration and immigration rates

#load packages
library(data.table); library(stringr); library(tidyr); library(dplyr); library(ggplot2)
library(lme4); library(lmerTest); library(readxl); library(DHARMa);library(glmmTMB); library(performance)

####### Reformatting for BA3 #######

male.stru <- fread("/Volumes/rchen2/Black Grouse PhD/Projects/Hunting_Microsats/Clean+final_analysis/Clean_data/Microsat.males.noLOCUS1+13.forstructure.stru")
head(male.stru)
names(male.stru) <- c("indivID", "popID", "BG16", "BG18", "BG15", "BG19", "BG6", "TTT1", "TTD2",
                      "TTD3", "TUD6", "TUT3", "TUT4", "TTT2")

male.stru.a <- structure[seq(from = 1, by = 2, to = nrow(male.stru)-1),]
male.stru.b <- structure[seq(from = 2, by = 2, to = nrow(male.stru)),]

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
female.stru <- fread("/Volumes/rchen2/Black Grouse PhD/Projects/Hunting_Microsats/Clean+final_analysis/Clean_data/Microsat.females.noLOCUS1+13.forstructure.stru")
head(female.stru)
names(female.stru) <- c("indivID", "popID", "BG16", "BG18", "BG15", "BG19", "BG6", "TTT1", "TTD2",
                      "TTD3", "TUD6", "TUT3", "TUT4", "TTT2")

female.stru.a <- structure[seq(from = 1, by = 2, to = nrow(female.stru)-1),]
female.stru.b <- structure[seq(from = 2, by = 2, to = nrow(female.stru)),]

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
male_run6 <- maleruns[[1]]
male_run7 <- maleruns[[2]]
male_run8 <- maleruns[[3]]
male_run9 <- maleruns[[4]]
male_run10 <- maleruns[[5]]

female_run6 <- femaleruns[[1]]
female_run7 <- femaleruns[[2]]
female_run8 <- femaleruns[[3]]
female_run9 <- femaleruns[[4]]
female_run10 <- femaleruns[[5]]

#### Compare runs ####
plot(male_run6$migration, male_run7$migration)
plot(male_run6$migration, male_run8$migration)
plot(male_run6$migration, male_run9$migration)
plot(male_run6$migration, male_run10$migration)
plot(male_run7$migration, male_run8$migration)
plot(male_run7$migration, male_run9$migration)
plot(male_run7$migration, male_run10$migration)
plot(male_run8$migration, male_run9$migration)
plot(male_run8$migration, male_run10$migration)
plot(male_run9$migration, male_run10$migration)
# seems like run 6, 8 and 10 correspond
# see if it's due to convergence: 
# with tracer count ESS < 100 for run 6 = 37
# with tracer count ESS < 100 for run 7 = 43

plot(female_run6$migration, female_run7$migration)
plot(female_run6$migration, female_run8$migration)
plot(female_run6$migration, female_run9$migration)
plot(female_run6$migration, female_run10$migration)
plot(female_run7$migration, female_run8$migration)
plot(female_run7$migration, female_run9$migration)
plot(female_run7$migration, female_run10$migration)
plot(female_run8$migration, female_run9$migration)
plot(female_run8$migration, female_run10$migration)
plot(female_run9$migration, female_run10$migration)

#seems like run 6, 7, 8, 10 (just not 9)

## Going to pick run 10 for both

#### Load in run 10 ####

#this file has been cleaned up, hunted status added, distances between sites and corrected migration value excluding those with ESS < 200

male_run10_select <- read.csv("data/run10_males_select_edited.csv")
female_run10_select <- read.csv("data/run10_females_select_edited.csv")

#correct site names that don't load in correctly
male_run10_select$pop_in[which(male_run10_select$pop_in == "Koskenp\xe4\xe4")] <- "Koskenpää"
male_run10_select$pop_in[which(male_run10_select$pop_in == "Nyr\xf6l\xe4")] <- "Nyrölä"
male_run10_select$pop_out[which(male_run10_select$pop_out == "Koskenp\xe4\xe4")] <- "Koskenpää"
male_run10_select$pop_out[which(male_run10_select$pop_out == "Nyr\xf6l\xe4")] <- "Nyrölä"

female_run10_select$pop_in[which(female_run10_select$pop_in == "Koskenp\xe4\xe4")] <- "Koskenpää"
female_run10_select$pop_in[which(female_run10_select$pop_in == "Nyr\xf6l\xe4")] <- "Nyrölä"
female_run10_select$pop_out[which(female_run10_select$pop_out == "Koskenp\xe4\xe4")] <- "Koskenpää"
female_run10_select$pop_out[which(female_run10_select$pop_out == "Nyr\xf6l\xe4")] <- "Nyrölä"

### Plotting migration rates from run 10 ####

theme_set(theme_classic())
# males: out hunted vs out unhunted but exclude nonmigration rates

ttest_males_in_rates <- t.test(male_run10_select$migration_ESSc[which(male_run10_select$hunt_in == "hunted" & male_run10_select$m_in != male_run10_select$m_out)], male_run10_select$migration_ESSc[which(male_run10_select$hunt_in == "unhunted" & male_run10_select$m_in != male_run10_select$m_out)])

ggplot(male_run10_select, aes(x = hunt_in, y = migration_ESSc)) + geom_boxplot(aes(fill = "Type")) +
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

ttest_males_out_rates <- t.test(male_run10_select$migration_ESSc[which(male_run10_select$hunt_out == "hunted")], male_run10_select$migration_ESSc[which(male_run10_select$hunt_out == "unhunted")])

ggplot(male_run10_select, aes(x = hunt_out, y = migration_ESSc, fill = "blue")) + geom_boxplot(aes(fill = "Type")) +
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
ttest_females_in_rates <- t.test(female_run10_select$migration_ESSc[which(female_run10_select$hunt_in == "hunted" & female_run10_select$m_in != female_run10_select$m_out)], female_run10_select$migration_ESSc[which(female_run10_select$hunt_in == "unhunted" & female_run10_select$m_in != female_run10_select$m_out)])

ggplot(female_run10_select, aes(x = hunt_in, y = migration_ESSc, fill = "red")) + geom_boxplot() +
  labs(title = "Migration rates IN for females") + xlab("Hunted status") + ylab("migration rates in")+
  labs(subtitle = 
         paste("Excluding values with ESS < 200 
Student t-test comparing hunted/unhunted site corrected migration rate in: 
t = ", round(ttest_females_in_rates$statistic,2), "and p-value = ", round(ttest_females_in_rates$p.value, 2)))+ 
  theme(text = element_text(family = "Arial", size = 22),
        plot.subtitle = element_text(size = 14),
        legend.position = "none")

# females: in hunted vs in unhunted but exclude nonmigration_ESSc rates

ttest_females_out_rates <- t.test(female_run10_select$migration_ESSc[which(female_run10_select$hunt_out == "hunted" & female_run10_select$m_in != female_run10_select$m_out)], female_run10_select$migration_ESSc[which(female_run10_select$hunt_out == "unhunted" & female_run10_select$m_in != female_run10_select$m_out)])

ggplot(female_run10_select, aes(x = hunt_out, y = migration_ESSc, fill = "red")) + geom_boxplot() +
  labs(title = "Migration rates OUT for females") + xlab("Hunted status") + ylab("migration rates out")+
  labs(subtitle = 
         paste("Excluding values with ESS < 200
Student t-test comparoutg hunted/unhunted site corrected migration rate out: 
t = ", round(ttest_females_out_rates$statistic,2), "and p-value = ", round(ttest_females_out_rates$p.value, 2)))+ 
  theme(text = element_text(family = "Arial", size = 22),
        plot.subtitle = element_text(size = 14),
        legend.position = "none")

#### Modelling migration ####

names(female_run10_select) == names(male_run10_select)
female_run10_select$sex <- "Female"
male_run10_select$sex <- "Male"

migration_both <- rbind(male_run10_select, female_run10_select)

#immigration

model.both.in <- glmmTMB(migration_ESSc~ hunt_in + Distance + sex + (1|pop_out) + (1|pop_in), 
                         data = migration_both, 
                         family = Gamma(link = "log"))
summary(model.both.in) 
plot(model.both.in)
simulateResiduals(fittedModel = model.both.in, plot = T)

#out
model.both.out <- glmmTMB(migration_ESSc~ hunt_out + Distance + sex + (1|pop_out) + (1|pop_in), 
                          data = migration_both, 
                          family = Gamma(link = "log"))
summary(model.both.out) 
plot(model.both.out)
simulateResiduals(fittedModel = model.both.out, plot = T)
