## In this script, we will first convert our dataframe to fit the format for BA3, then run BA3,
## analyze the 5 different runs for the migration models
## as constructed with BayeSass v3, and then we will model the effect of hunting
## on both emigration and immigration rates

#load packages
library(data.table); library(tidyverse); library(tibble); library(MuMIn)
library(lme4); library(lmerTest); library(readxl); library(DHARMa);library(glmmTMB); library(performance)

# load in pop data

pops <- read.csv("data/details/Codes.pops.both.filtered_withcoord.csv")
pops$pop_num <- as.factor(pops$pop_num)
pops[1,1] <- "Koskenpää"
pops[5,1] <- "Nyrölä"

## load in Matrix distances between sites ##
distance <- read_excel("data/details/CalculateDistanceSitesGenAlEx.xlsx", sheet = "MatrixForR")
names(distance)[1] <- "Site_A"

distance_long <- melt(distance)
names(distance_long) <- c("Site_A", "Site_B", "Distance")
distance_long <- subset(distance_long, distance_long$Site_A != distance_long$Site_B)

###### Before BA3 - pop assignment test #####

assign <- read_xlsx("analyses/genalex/GenAlEx.adults.assignment.noLOCUS1+13.xlsx", sheet = "ASS",
                    skip = 11)
assign <- assign[-c(13,14),] #take out total
assign$Total <- assign$`Self Pop`+assign$`Other Pop`
assign$prop_migrated <- assign$`Other Pop`/assign$Total

#add hunted/unhunted
assign <- left_join(assign, pops[,c(1,3)], by = c("Pop" = "pop"))
assign$hunt <- as.factor(assign$hunt)
#t.test
t.test(assign$prop_migrated[which(assign$hunt == "hunted")], 
                            assign$prop_migrated[which(assign$hunt == "unhunted")])

#model
model_assign <- lm(prop_migrated ~ hunt + Total, data = assign)
summary(model_assign)

######## BA3 ##########
#######################

####### Reformatting for BA3 #######

all.stru <- fread("data/cleandata/Microsat.adults.noLOCUS1+13.forstructure.stru")
head(all.stru)
names(all.stru) <- c("indivID", "popID", "BG16", "BG18", "BG15", "BG19", "BG6", "TTT1", "TTD2",
                      "TTD3", "TUD6", "TUT3", "TUT4", "TTT2")

all.stru[all.stru==-9] <- 0
all.stru.a <- all.stru[seq(from = 1, by = 2, to = nrow(all.stru)-1),]
all.stru.b <- all.stru[seq(from = 2, by = 2, to = nrow(all.stru)),]

ba3.all.a <- melt(data = all.stru.a,
              id.vars = c("indivID", "popID"),
              variable.name = "locID",
              value.name = "allele1")

ba3.all.b <- melt(data = all.stru.b,
              id.vars = c("indivID", "popID"),
              variable.name = "locID",
              value.name = "allele2")

ba3.all <- left_join(ba3.all.a, ba3.all.b, by = c("indivID", "popID", "locID"))

head(ba3.all)

write.table(ba3.all, "analyses/migrationanalysis/data_all_ba3.txt",
            col.names = T, row.names = F, sep = " ", quote = F)

#### Running BA3 ####
#make directory per run
system(paste0("mkdir ", getwd(), "/analyses/migrationanalysis/BA3runs/run1")) 
system(paste0("mkdir ", getwd(), "/analyses/migrationanalysis/BA3runs/run2"))
system(paste0("mkdir ", getwd(), "/analyses/migrationanalysis/BA3runs/run3"))
system(paste0("mkdir ", getwd(), "/analyses/migrationanalysis/BA3runs/run4"))
system(paste0("mkdir ", getwd(), "/analyses/migrationanalysis/BA3runs/run5"))

#5 runs with 5 different random seeds
pathba3 <- "/prj/blackgrouse/bin/"
pathba3 <- "/Users/vistor/Documents/Work/Bielefeld/PhD/Software/BA3-migration/"
pathba3 <- "~/"
system(paste0(pathba3, "BA3/BA3MSAT -v -t -g -u -a 0.30 -f 0.40 -s 65323 -i 10000000 -b 1000000 -n 1000 -o run1.txt ", getwd(), "/analyses/migrationanalysis/data_females_ba3.txt")) 
system(paste0(pathba3, "BA3/BA3MSAT -v -t -g -u -a 0.30 -f 0.40 -s 76553 -i 10000000 -b 1000000 -n 1000 -o run2.txt ", getwd(), "/analyses/migrationanalysis/data_all_ba3.txt")) 
system(paste0(pathba3, "BA3/BA3MSAT -v -t -g -u -a 0.30 -f 0.40 -s 124643 -i 10000000 -b 1000000 -n 1000 -o run3.txt ", getwd(), "/analyses/migrationanalysis/data_all_ba3.txt")) 
system(paste0(pathba3, "BA3/BA3MSAT -v -t -g -u -a 0.30 -f 0.40 -s 885256 -i 10000000 -b 1000000 -n 1000 -o run4txt ", getwd(), "/analyses/migrationanalysis/data_all_ba3.txt")) 
system(paste0(pathba3, "BA3/BA3MSAT -v -t -g -u -a 0.30 -f 0.40 -s 235776 -i 10000000 -b 1000000 -n 1000 -o run5.txt ", getwd(), "/analyses/migrationanalysis/data_all_ba3.txt")) 

#### Compare all 10 runs ####
temp <- list.files(path = "analyses/migrationanalysis/BA3runs/", pattern = ".txt", full.names=T)
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


#separate per run to compare
male_run1 <- maleruns[[1]]
male_run2 <- maleruns[[2]]
male_run3 <- maleruns[[3]]
male_run4 <- maleruns[[4]]
male_run5 <- maleruns[[5]]

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

## Going to pick run 5 for both

#### Clean migration files ####
setwd("data")
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

write.table(male_run5_clean, file = "data/migrationanalysis/run5_males_clean.csv", sep = ",", row.names = F,quote = F)
write.table(female_run5_clean, file = "data/migrationanalysis/run5_females_clean.csv", sep = ",", row.names = F,quote = F)

#### Load in clean runs ####

#this file has been cleaned up, hunted status added, distances between sites and corrected migration value excluding those with ESS < 200

male_run5_clean <- read.csv("data/migrationanalysis/run5_males_clean.csv")
female_run5_clean <- read.csv("data/migrationanalysis/run5_females_clean.csv")

## write out as matrix for supplements ##
males_run5_m <- male_run5_clean[,c(7,9,6)]
males_run5_m <- pivot_wider(males_run5_m, names_from = "pop_out", values_from = "migration_ESSc")
males_run5_se<- male_run5_clean[,c(7,9,4)]
males_run5_se <- pivot_wider(males_run5_se, names_from = "pop_out", values_from = "migration_SE")
males_run5_m_raw<- male_run5_clean[,c(7,9,3)]
males_run5_m_raw <- pivot_wider(males_run5_m_raw, names_from = "pop_out", values_from = "migration")

females_run5_m <- female_run5_clean[,c(7,9,6)]
females_run5_m <- pivot_wider(females_run5_m, names_from = "pop_out", values_from = "migration_ESSc")
females_run5_se<- female_run5_clean[,c(7,9,4)]
females_run5_se <- pivot_wider(females_run5_se, names_from = "pop_out", values_from = "migration_SE")
females_run5_m_raw<- female_run5_clean[,c(7,9,3)]
females_run5_m_raw <- pivot_wider(females_run5_m_raw, names_from = "pop_out", values_from = "migration")

write.csv(males_run5_m, "data/migrationanalysis/males.migration.csv", row.names = F)
write.csv(males_run5_se, "data/migrationanalysis/males.migration.se.csv", row.names = F)
write.csv(males_run5_m_raw, "data/migrationanalysis/males.migration.raw.csv", row.names = F)

write.csv(females_run5_m, "data/migrationanalysis/females.migration.csv", row.names = F)
write.csv(females_run5_se, "data/migrationanalysis/females.migration.se.csv", row.names = F)
write.csv(females_run5_m_raw, "data/migrationanalysis/females.migration.raw.csv", row.names = F)

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

model.both.in.null <- glmmTMB(migration_ESSc~ Distance + sex + (1|pop_out) + (1|pop_in), 
                         data = migration_both, 
                         family = Gamma(link = "log"))

anova(model.both.in.null, model.both.in)

summary(model.both.in) 
simulateResiduals(fittedModel = model.both.in, plot = T)

#with interaction
model.both.in.interaction <- glmmTMB(migration_ESSc~ hunt_in + Distance*sex + (1|pop_out) + (1|pop_in), 
                         data = migration_both, 
                         family = Gamma(link = "log"))
model.both.in.interaction.null <- glmmTMB(migration_ESSc~ Distance*sex + (1|pop_out) + (1|pop_in), 
                                     data = migration_both, 
                                     family = Gamma(link = "log"))
anova(model.both.in.interaction, model.both.in.interaction.null)

summary(model.both.in.interaction) 
simulateResiduals(fittedModel = model.both.in.interaction, plot = T)

#out
model.both.out <- glmmTMB(migration_ESSc~ hunt_out + Distance + sex + (1|pop_out) + (1|pop_in), 
                          data = migration_both, 
                          family = Gamma(link = "log"))
model.both.out.null <- glmmTMB(migration_ESSc~ Distance + sex + (1|pop_out) + (1|pop_in), 
                          data = migration_both, 
                          family = Gamma(link = "log"))

anova(model.both.out, model.both.out.null)

summary(model.both.out) 
simulateResiduals(fittedModel = model.both.out, plot = T)

#with interaction
model.both.out.interaction <- glmmTMB(migration_ESSc~ hunt_out + Distance*sex + (1|pop_out) + (1|pop_in), 
                          data = migration_both, 
                          family = Gamma(link = "log"))
model.both.out.interaction.null <- glmmTMB(migration_ESSc~ Distance*sex + (1|pop_out) + (1|pop_in), 
                                      data = migration_both, 
                                      family = Gamma(link = "log"))
anova(model.both.out.interaction, model.both.out.interaction.null)

summary(model.both.out.interaction) 
simulateResiduals(fittedModel = model.both.out.interaction, plot = T)

## assess model performance
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
