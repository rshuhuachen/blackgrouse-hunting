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

######## BA3 ##########
#######################

####### Reformatting for BA3 with BG20 #######
all.stru.nohwe <- fread("data/cleandata/Microsat.adults.noLOCUS1.forstructure.stru")
head(all.stru.nohwe)
names(all.stru.nohwe) <- c("indivID", "popID", "BG16", "BG18", "BG15", "BG19", "BG6", "TTT1", "TTD2",
                     "TTD3", "TUD6", "TUT3", "TUT4", "BG20", "TTT2")

all.stru.nohwe[all.stru.nohwe==-9] <- 0
all.stru.nohwe.a <- all.stru.nohwe[seq(from = 1, by = 2, to = nrow(all.stru.nohwe)-1),]
all.stru.nohwe.b <- all.stru.nohwe[seq(from = 2, by = 2, to = nrow(all.stru.nohwe)),]

ba3.all.nohwe.a <- melt(data = all.stru.nohwe.a,
                  id.vars = c("indivID", "popID"),
                  variable.name = "locID",
                  value.name = "allele1")

ba3.all.nohwe.b <- melt(data = all.stru.nohwe.b,
                  id.vars = c("indivID", "popID"),
                  variable.name = "locID",
                  value.name = "allele2")

ba3.nohwe.all <- left_join(ba3.all.nohwe.a, ba3.all.nohwe.b, by = c("indivID", "popID", "locID"))

head(ba3.nohwe.all)

write.table(ba3.nohwe.all, "analyses/migrationanalysis/data_all_ba3_allloci.txt",
            col.names = F, row.names = F, sep = " ", quote = F)

#side note: mac can give segmentation error when running BA3 as this file is not in a location it can access (github dir?)
#if this happens, just move the .txt file with data to another folder and adjust the location accordingly in the BA3 commands below

#### Running BA3 ####
# do 5 runs with 5 different random seeds

#(1) migration rates; (2) individual migrant ancestries; (3) allele frequencies; (4) inbreeding coefficients; (5) missing genotypes
pathba3 <- "/Users/vistor/Documents/Work/Bielefeld/PhD/Software/BA3-migration/" #path to BA3
system(paste0(pathba3, "BA3/BA3MSAT -v -t -g -u -a 0.30 -f 0.40 -s 65323 -i 10000000 -b 1000000 -n 500 -o run1_nohwe.txt ", getwd(), "/analyses/migrationanalysis/data_all_ba3_allloci.txt")) 
system(paste0(pathba3, "BA3/BA3MSAT -v -t -g -u -a 0.30 -f 0.40 -s 235776 -i 10000000 -b 1000000 -n 500 -o run2_nohwe.txt ", getwd(), "/analyses/migrationanalysis/data_all_ba3_allloci.txt")) 
system(paste0(pathba3, "BA3/BA3MSAT -v -t -g -u -a 0.30 -f 0.40 -s 124643 -i 10000000 -b 1000000 -n 500 -o run3_nohwe.txt ", getwd(), "/analyses/migrationanalysis/data_all_ba3_allloci.txt")) 
system(paste0(pathba3, "BA3/BA3MSAT -v -t -g -u -a 0.30 -f 0.40 -s 885256 -i 10000000 -b 1000000 -n 500 -o run4_nohwe.txt ", getwd(), "/analyses/migrationanalysis/data_all_ba3_allloci.txt")) 
system(paste0(pathba3, "BA3/BA3MSAT -v -t -g -u -a 0.30 -f 0.40 -s 76553 -i 10000000 -b 1000000 -n 500 -o run5_nohwe.txt ", getwd(), "/analyses/migrationanalysis/data_all_ba3_allloci.txt")) 

#### Compare all 5 runs ####
temp <- list.files(path = "analyses/migrationanalysis/BA3runs/", pattern = "*nohwe.txt", full.names=T)
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
run1 <- myfiles[[1]]
run2 <- myfiles[[2]]
run3 <- myfiles[[3]]
run4 <- myfiles[[4]]
run5 <- myfiles[[5]]

#### Compare runs ####
plot(run1$migration, run2$migration)
plot(run1$migration, run3$migration)
plot(run1$migration, run4$migration)
plot(run1$migration, run5$migration)
plot(run2$migration, run3$migration)
plot(run2$migration, run4$migration)
plot(run2$migration, run5$migration)
plot(run3$migration, run4$migration)
plot(run3$migration, run5$migration)
plot(run4$migration, run5$migration)

# all runs correspond
## Going to pick run 5

#### Clean migration files ####
# to do: add ESS, corrected migration value, hunted status, distance between sites

#first add ESS
ESS_run5 <- read.delim("analyses/migrationanalysis/BA3runs/run5_nohwe_BA3_tracersummary.txt", sep = "\t")

ESS_run5.df <- t(ESS_run5)
colnames(ESS_run5.df) <- ESS_run5.df[1,]
ESS_run5.df <- as.data.frame(ESS_run5.df[-c(1,2),])#take out col names and row for log prob
ESS_run5.df <- rownames_to_column(ESS_run5.df, "m")

run5 <- run5 %>% arrange(m_in) #now, the two files are sorted the same way, so can just cbind
run5_clean <- cbind(run5, ESS_run5.df[,c(11)]) #select ESS

head(run5_clean)
names(run5_clean)[5] <- "ESS"

#add migration_ESSc column, which turns to NA if ESS < 200

run5_clean$ESS <- as.numeric(run5_clean$ESS)
run5_clean$migration <- as.numeric(run5_clean$migration)
run5_clean$migration_SE <- as.numeric(run5_clean$migration_SE)
run5_clean$migration_ESSc <- case_when(run5_clean$ESS < 200 ~ as.numeric(NA),
                                            run5_clean$ESS >= 200 ~ as.numeric(run5_clean$migration))

# then, add hunted status and full site names
run5_clean$m_in<-as.factor(run5_clean$m_in)
run5_clean$m_out<-as.factor(run5_clean$m_out)

run5_clean <- left_join(run5_clean, pops[,c(1:3)], by = c("m_in" = "pop_num"))
names(run5_clean)[7] <- "pop_in"
names(run5_clean)[8] <- "hunt_in"

run5_clean <- left_join(run5_clean, pops[,c(1:3)], by = c("m_out" = "pop_num"))
names(run5_clean)[9] <- "pop_out"
names(run5_clean)[10] <- "hunt_out"

## add distance
run5_clean <- left_join(run5_clean, distance_long, by = c("pop_in" = "Site_A", "pop_out" = "Site_B"))

write.table(run5_clean, file = "analyses/migrationanalysis/run5_clean.csv", sep = ",", row.names = F,quote = F)

#### Load in clean runs ####

#this file has been cleaned up, hunted status added, distances between sites and corrected migration value excluding those with ESS < 200

run5_clean <- read.csv("analyses/migrationanalysis/run5_clean.csv")

## write out as matrix for supplements ##
run5_m <- run5_clean[,c(7,9,6)]
run5_m <- pivot_wider(run5_m, names_from = "pop_out", values_from = "migration_ESSc")
run5_se<- run5_clean[,c(7,9,4)]
run5_se <- pivot_wider(run5_se, names_from = "pop_out", values_from = "migration_SE")
run5_m_raw<- run5_clean[,c(7,9,3)]
run5_m_raw <- pivot_wider(run5_m_raw, names_from = "pop_out", values_from = "migration")

write.csv(run5_m, "analyses/migrationanalysis/migration.csv", row.names = F)
write.csv(run5_se, "analyses/migrationanalysis/migration.se.csv", row.names = F)
write.csv(run5_m_raw, "analyses/migrationanalysis/migration.raw.csv", row.names = F)

### Plotting migration rates from run 5 ####
#first, exclude the 'non-migration rates' which are those where pop in = pop out
run5_clean <- subset(run5_clean, m_in != m_out)

theme_set(theme_classic())
# out hunted vs out unhunted but exclude nonmigration rates

ttest_in_rates <- t.test(run5_clean$migration_ESSc[which(run5_clean$hunt_in == "hunted" & run5_clean$m_in != run5_clean$m_out)], run5_clean$migration_ESSc[which(run5_clean$hunt_in == "unhunted" & run5_clean$m_in != run5_clean$m_out)])

ggplot(run5_clean, aes(x = hunt_in, y = migration_ESSc)) + geom_boxplot(aes(fill = "Type")) +
  labs(title = "Migration rates IN") + xlab("Hunted status") + ylab("migration rates in")+
  labs(subtitle = 
         paste("Excluding values with ESS < 200 
Student t-test comparing hunted/unhunted site filtered migration rate in: 
t = ", round(ttest_in_rates$statistic,2), "and p-value = ", round(ttest_in_rates$p.value, 2)))+ 
  theme(text = element_text(family = "Arial", size = 22),
        plot.subtitle = element_text(size = 14),
        legend.position = "none")+
  scale_fill_manual(values = c("cyan3"))

# in hunted vs in unhunted but exclude nonmigration_ESSc rates

ttest_out_rates <- t.test(run5_clean$migration_ESSc[which(run5_clean$hunt_out == "hunted")], run5_clean$migration_ESSc[which(run5_clean$hunt_out == "unhunted")])

ggplot(run5_clean, aes(x = hunt_out, y = migration_ESSc, fill = "blue")) + geom_boxplot(aes(fill = "Type")) +
  labs(title = "Migration rates OUT") + xlab("Hunted status") + ylab("migration rates out")+
  labs(subtitle = 
         paste("Excluding values with ESS < 200 
Student t-test comparing hunted/unhunted site corrected migration rate out: 
t = ", round(ttest_out_rates$statistic,2), "and p-value = ", round(ttest_out_rates$p.value, 3)))+ 
  theme(text = element_text(family = "Arial", size = 22),
        plot.subtitle = element_text(size = 14),
        legend.position = "none")+
  scale_fill_manual(values = c("cyan3"))

#### Modelling migration ####

#change levels hunted/unhunted
run5_clean$hunt_in <- relevel(as.factor(run5_clean$hunt_in), ref = "unhunted")
run5_clean$hunt_out <- relevel(as.factor(run5_clean$hunt_out), ref = "unhunted")

#immigration

model.in <- glmmTMB(migration_ESSc~ hunt_in + Distance +  (1|pop_out) + (1|pop_in), 
                         data = run5_clean, 
                         family = Gamma(link = "log"))

model.in.null <- glmmTMB(migration_ESSc~ Distance +  (1|pop_out) + (1|pop_in), 
                         data = run5_clean, 
                         family = Gamma(link = "log"))

anova(model.in, model.in.null)

summary(model.in) 
simulateResiduals(fittedModel = model.in, plot = T)
icc(model.in, by_group=T)
r.squaredGLMM(model.in)

testUniformity(model.in)
testOutliers(model.in)
testOverdispersion(model.in)

#out

model.out <- glmmTMB(migration_ESSc~ hunt_out + Distance + (1|pop_out) + (1|pop_in), 
                          data = run5_clean, 
                          family = Gamma(link = "log"))
model.out.null <- glmmTMB(migration_ESSc~ Distance + (1|pop_out) + (1|pop_in), 
                          data = run5_clean, 
                          family = Gamma(link = "log"))

anova(model.out, model.out.null)

summary(model.out) 
simulateResiduals(fittedModel = model.out, plot = T)
icc(model.out, by_group = T)
r.squaredGLMM(model.out)

testUniformity(model.out)
testOutliers(model.out)
testOverdispersion(model.out)

compare_performance(model.in, model.in.null, rank=T)
compare_performance(model.out, model.out.null, rank=T)

#### Additional: assignment test #####

### Done with geneclass2, rannala&mountain(97) criterion, 1000 simulated individuals with Paetkau et al 2004 algorithm. 0.01 p value threshold
library(data.table);library(dplyr)
assign<-fread("analyses/migrationanalysis/geneclass_results.csv", skip=13)
migrants <-subset(assign, probability <= 0.01)
migrants$id_n <- row.names(migrants)
migrants <- migrants[,c(20,2:17)]
names(migrants) <- c("id_n", "site", "home_max", "pval", "1", "2", "3", "4", 
                     "5", "6", "7", "8", "9", "10", "11", "12", "n_loci")

#add location of lowest -log(L)
#transform from wide to long
migrants_long <- melt(setDT(migrants), id.vars = c("id_n", "site", "home_max", "pval", "n_loci"), variable.name = "potential_source")
migrants_long_min <- migrants_long %>% group_by(id_n) %>% filter(value == min(value))
migrants_long_min$site <- as.factor(migrants_long_min$site)
migrants_long_min$potential_source <- as.factor(migrants_long_min$potential_source)
## add in hunted or not
pops <- read.csv("data/details/Codes.pops.both.filtered_withcoord.csv")
pops$pop_num <- as.factor(pops$pop_num)
pops[1,1] <- "Koskenpää"
pops[5,1] <- "Nyrölä"

#merge
migrants_long_min <- left_join(migrants_long_min, pops[,c(1:3)], by = c("site" = "pop_num"))
migrants_long_min <- left_join(migrants_long_min, pops[,c(1:3)], by = c("potential_source" = "pop_num"))
names(migrants_long_min)[8] <- "site_name"
names(migrants_long_min)[9] <- "site_hunted"
names(migrants_long_min)[10] <- "source_name"
names(migrants_long_min)[11] <- "source_hunted"

migrants_long_min %>% group_by(site_hunted, source_hunted) %>% count()
write.csv(migrants_long_min, "analyses/migrationanalysis/geneclass_migrants.csv",row.names = F, quote=F)

# left out assignment test as likely our data are not suitable for an individual based migration analysis with high gene flow, with only 12 msats
