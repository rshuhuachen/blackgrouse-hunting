#### Population Structure #####

# In this script, we will 1) calculate genetic summary statistics, 2) make a PCA, 
# 3) calculate pairwise Fst

## load libraries
pacman::p_load(data.table, tidyverse, hierfstat, plot.matrix, lme4, adegenet, forcats, ape, readxl, tibble)

all.raw <- read.structure ("data/cleandata/Microsat.adults.plus.unrelated.chicks.noLOCUS1.forstructure.stru", n.ind = 2078, n.loc = 13, onerowperind = F,
                       col.lab = 1, col.pop = 2, col.others = NULL,
                       row.marknames = 0, NA.char = "-9", pop = NULL, sep = NULL,
                       ask = F, quiet = FALSE)

all <- read.structure ("data/cleandata/Microsat.adults.plus.unrelated.chicks.noLOCUS1+13.forstructure.stru", n.ind = 2078, n.loc = 12, onerowperind = F,
                                col.lab = 1, col.pop = 2, col.others = NULL,
                                row.marknames = 0, NA.char = "-9", pop = NULL, sep = NULL,
                                ask = F, quiet = FALSE)

pops <- read.csv("data/details/Codes.pops.both.filtered_withcoord.csv")
pops[1,1] <- "Koskenpää"
pops[5,1] <- "Nyrölä"

# males.data <- fread("data/cleandata/Unsplit.microsat.males.noLOCUS1+13.csv")
# females.data <- fread("data/cleandata/Unsplit.microsat.females.noLOCUS1+13.csv")
# chicks.data <- fread("data/cleandata/Unsplit.microsat.unrelated.chicks.noLOCUS1+13+14.csv")

#### Summary statistics ####

# This is based on all individuals and all loci
basicstat.all <- basic.stats(all.raw, diploid = TRUE, digits = 2) 
allelic.richness.all <- allelic.richness(all.raw, diploid = TRUE)

#Ar
allelic.richness.all.df <- as.data.frame(allelic.richness.all$Ar)

#get mean Ho and Hx per pop
x.pop = seppop(all.raw) 
summary.by.pop = lapply(x.pop, summary) 
Hobs.ls = rep(NA, length(summary.by.pop)) 
for (i in 1:length(summary.by.pop)){ 
  Hobs.ls[i] = mean(summary.by.pop[[i]]$Hobs, na.rm=TRUE) 
} 
Hobs.ls

Hexp.ls = rep(NA, length(summary.by.pop)) 
for (i in 1:length(summary.by.pop)){ 
  Hexp.ls[i] = mean(summary.by.pop[[i]]$Hexp, na.rm=TRUE) 
} 
Hexp.ls

#### Compare genetic diversity between hunted vs unhunted ####
pops$pop_num <- as.character(pops$pop_num)
measures <- data.frame(pop = names(summary.by.pop), Hobs = Hobs.ls, Hexp = Hexp.ls, Ar = colMeans(allelic.richness.all.df, na.rm = T))
measures <- left_join(measures, pops[,c(2,3)], by = c("pop" = "pop_num"))

t.test(measures$Hobs[which(measures$hunt == "hunted")], measures$Hobs[which(measures$hunt == "unhunted")])
t.test(measures$Hexp[which(measures$hunt == "hunted")], measures$Hexp[which(measures$hunt == "unhunted")])
t.test(measures$Ar[which(measures$hunt == "hunted")], measures$Ar[which(measures$hunt == "unhunted")])

# let's also do it including all data per locus, and model with locus as random factor
Hobs.all <- NULL

for (i in 1:length(summary.by.pop)){
  Hobs.i <- data.frame(locus= c(names(summary.by.pop[[i]]$Hobs)),
                       pop = names(summary.by.pop[i]),
                       Hobs = c(summary.by.pop[[i]]$Hobs))
  Hobs.all <- rbind(Hobs.all, Hobs.i)
  rownames(Hobs.all) <- NULL
}

Hexp.all <- NULL

for (i in 1:length(summary.by.pop)){
  Hexp.i <- data.frame(locus= c(names(summary.by.pop[[i]]$Hexp)),
                       pop = names(summary.by.pop[i]),
                       Hexp = c(summary.by.pop[[i]]$Hexp))
  Hexp.all <- rbind(Hexp.all, Hexp.i)
  rownames(Hexp.all) <- NULL
}

Ar.all <- NULL

for (i in 1:ncol(allelic.richness.all.df)){
  Ar.i <- data.frame(locus= c(rownames(allelic.richness.all.df[i])),
                       pop = colnames(allelic.richness.all.df[i]),
                       Ar = c(allelic.richness.all.df[,i]))
  Ar.all <- rbind(Ar.all, Ar.i)
  rownames(Ar.all) <- NULL
}

measures.all <- left_join(Hobs.all, Hexp.all) %>% left_join(Ar.all) %>% left_join(pops[,c(2,3)], by = c("pop" = "pop_num"))

measures.all
measures.all$locus <- as.factor(measures.all$locus)
measures.all$pop <- as.factor(measures.all$pop)
measures.all$hunt <- as.factor(measures.all$hunt)

Ho_model <- lmerTest::lmer(Hobs ~ hunt + (1|locus) + (1|pop), data = subset(measures.all, locus != "L12")) # excluding BG20
Ho_model_null <- lmerTest::lmer(Hobs ~ (1|locus) + (1|pop), data = subset(measures.all, locus != "L12")) # excluding BG20

He_model <- lmerTest::lmer(Hexp ~ hunt + (1|locus) + (1|pop), data = subset(measures.all, locus != "L12"))
He_model_null <- lmerTest::lmer(Hexp ~ (1|locus) + (1|pop), data = subset(measures.all, locus != "L12"))

Ar_model <- lmerTest::lmer(Ar ~ hunt + (1|locus) + (1|pop), data = subset(measures.all, locus != "L12"))
Ar_model_null <- lmerTest::lmer(Ar ~ (1|locus) + (1|pop), data = subset(measures.all, locus != "L12"))

#with BG20
# Ho_model <- lmerTest::lmer(Hobs ~ hunt + (1|locus) + (1|pop), data = measures.all) # excluding BG20
# Ho_model_null <- lmerTest::lmer(Hobs ~ (1|locus) + (1|pop), data = measures.all) # excluding BG20
# 
# He_model <- lmerTest::lmer(Hexp ~ hunt + (1|locus) + (1|pop), data =  measures.all)
# He_model_null <- lmerTest::lmer(Hexp ~ (1|locus) + (1|pop), data =  measures.all)
# 
# Ar_model <- lmerTest::lmer(Ar ~ hunt + (1|locus) + (1|pop), data =  measures.all)
# Ar_model_null <- lmerTest::lmer(Ar ~ (1|locus) + (1|pop), data =  measures.all)

#LRT
anova(Ho_model, Ho_model_null)
anova(He_model, He_model_null)
anova(Ar_model, Ar_model_null)

#### Calculate Fst ####

## Calculated based on all adults + unrelated chicks HWE filtered

#Total
all.hfstat <- genind2hierfstat(all)
#calculate stats
basicstat.all <- basic.stats(all, diploid = TRUE, digits = 2) # Hobs, mean gene diversities, Fis and Fst

# per locus
fst.all.perlocus <- basicstat.all$perloc$Fst
fst.all.perlocus <- data.frame(Locus = seq(from = 1, to = 12), Fst = fst.all.perlocus)

# Pairwise Fst
fst.all <- pairwise.neifst(all.hfstat)
head(fst.all)

# Fst per population
boxplot(fst.all, col=funky(nPop(all)), las=3,
        xlab="Population", ylab="Fst", main = "Pairwise Fst values per population") 

#Bootstrap
boot.fst.all <- boot.ppfst(all.hfstat, nboot = 1000) #confidence interval for Fst, upper triangle is upper limit, lower triangle lower limit

#create a long dataframe for pairwise Fst
boot.fst.all.UL <- boot.fst.all$ul
boot.fst.all.LL <- boot.fst.all$ll

flat.matrix <- function(d){
  data.frame(i=rep(row.names(d),ncol(d)),
             j=rep(colnames(d),each=nrow(d)),
             score=as.vector(d))
}

fst.all.flat <- flat.matrix(fst.all)
names(fst.all.flat) <- c("site.x", "site.y", "Fst")

boot.fst.all.LL.flat <- flat.matrix(boot.fst.all.LL)
names(boot.fst.all.LL.flat) <- c("site.x", "site.y", "LL")

boot.fst.all.UL.flat <- flat.matrix(boot.fst.all.UL)
names(boot.fst.all.UL.flat) <- c("site.x", "site.y", "UL")

pairwise.fst.all <- left_join(fst.all.flat, boot.fst.all.LL.flat, by = c("site.x", "site.y"))
pairwise.fst.all <- left_join(pairwise.fst.all, boot.fst.all.UL.flat, by = c("site.x", "site.y"))

pairwise.fst.all <- subset(pairwise.fst.all, site.x != "-9" & site.y != "-9")
pairwise.fst.all <- subset(pairwise.fst.all, !is.na(Fst) & !is.na(UL))
pairwise.fst.all <- pairwise.fst.all %>% mutate(Significance = case_when(
  UL > 0 & LL > 0 ~ "significant",
  UL > 0 & LL < 0 ~ "insignificant",
  UL < 0 & LL < 0 ~ "significant" ))

pairwise.fst.all <- pairwise.fst.all %>% mutate(Sig = case_when(
  UL > 0 & LL > 0 ~ "",
  UL > 0 & LL < 0 ~ "NS",
  UL < 0 & LL < 0 ~ "" ))

write.csv(pairwise.fst.all, "tables/Pairwise_Fst_all.csv", row.names = F, quote=F)

#### Spatial autocorrelation #####

# add site names rather than numbers
males.data <- left_join(males.data, pops[,-3], by = c("pop" = "pop_num"))
females.data <- left_join(females.data, pops[,-3], by = c("pop" = "pop_num"))
chicks.data <- left_join(chicks.data, pops[,-3], by = c("pop" = "pop_num"))

males.data <- males.data[,c(1,27,3:26,28,29)] #select only pop in letters not numbers
names(males.data)[2] <- "SITE"
names(males.data)[1]<- "CODE"

females.data <- females.data[,c(1,27,3:26,28,29)] #select only pop in letters not numbers
names(females.data)[2] <- "SITE"
names(females.data)[1]<- "CODE"

chicks.data <- chicks.data[,c(1,25,3:24,26,27)] #select only pop in letters not numbers
names(chicks.data)[2] <- "SITE"
names(chicks.data)[1]<- "CODE"

#no coordinates for Lauttasuo (adults)
# both in adults and chicks, some have missing locations and thus missing coordinates

### Change missing data to 0 for GenAIEx 
males.data[is.na(males.data)] <- 0
males.data[(males.data == -9)] <- 0

females.data[is.na(females.data)] <- 0
females.data[(females.data == -9)] <- 0

chicks.data[is.na(chicks.data)] <- 0
chicks.data[(chicks.data == -9)] <- 0

summary(as.factor(males.data$SITE))
summary(as.factor(females.data$SITE))
summary(as.factor(chicks.data$SITE))

#combine adult data for sex-biased dispersal tests
males.df <- males.data %>% add_column(sex = "M")
females.df <- females.data %>% add_column(sex = "F")
adults.data <- rbind(males.df, females.df)

#have to take out all rows with missing data
adults.sexbias <- adults.data[,c(1,29,3:26)]
adults.sexbias[(adults.sexbias == 0)] <- NA
adults.sexbias <- adults.sexbias[,-c(15,16,25,26)] #take out two loci with lots of missing data 
adults.sexbias.complete <- adults.sexbias[complete.cases(adults.sexbias),] #1861 N with complete cases

## Then format exactly as GenAIEx requires in excel manually
# A1: NO OF LOCI, B2: NO OF SAMPLES, C1: NO OF POPULATIONS, D1-N1: SIZE OF EACH POPULATION
# A2: optional title, D2: F2: pop labels
# headers A3:N3: CODE, SITE, LOCUS1, LOCUS 1, LOCUS 2, LOCUS 2, [EMPTY COL], X, Y
# Data starts on C4 with microsats

write.table(males.data, file = "analyses/genalex/msat_males_withcoord.csv",quote=F, row.names=F)
write.table(females.data, file = "analyses/genalex/msat_females_withcoord.csv",quote=F, row.names=F)
write.table(chicks.data, file = "analyses/genalex/msat_chicks_withcoord.csv",quote=F, row.names=F)
write.table(adults.data, file = "analyses/genalex/msat_adults_withcoord.csv",quote=F, row.names=F)
write.table(adults.sexbias.complete, file = "analyses/genalex/msat_adults_sexbias_withcoord.csv",quote=F, row.names=F)






