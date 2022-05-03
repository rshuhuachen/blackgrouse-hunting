#### Population Structure #####

# In this script, we will 1) calculate genetic summary statistics, 2) make a PCA, 
# 3) calculate pairwise Fst

## load libraries
library(data.table); library(tidyverse); library(hierfstat); 
library(plot.matrix); library(lme4); library(adegenet); library(forcats)
library(ape)

males.stru <- read.structure("data/cleandata/Microsat.males.noLOCUS1+13.forstructure.stru", n.ind = 1065, n.loc = 12, onerowperind = F,
                             col.lab = 1, col.pop = 2, col.others = NULL,
                             row.marknames = 0, NA.char = "-9", pop = NULL, sep = NULL,
                             ask = F, quiet = FALSE)

females.stru <- read.structure("data/cleandata/Microsat.females.noLOCUS1+13.forstructure.stru", n.ind = 813, n.loc = 12, onerowperind = F,
                               col.lab = 1, col.pop = 2, col.others = NULL,
                               row.marknames = 0, NA.char = "-9", pop = NULL, sep = NULL,
                               ask = F, quiet = FALSE)

chicks.stru <- read.structure("data/cleandata/Microsat.chicks.noLOCUS1+13+14.forstructure.stru", n.ind = 1370, n.loc = 11, onerowperind = F,
                              col.lab = 1, col.pop = 2, col.others = NULL,
                              row.marknames = 0, NA.char = "-9", pop = NULL, sep = NULL,
                              ask = F, quiet = FALSE)

all <- read.structure ("data/cleandata/Microsat.all.stru", n.ind = 3248, n.loc = 14, onerowperind = F,
                       col.lab = 1, col.pop = 2, col.others = NULL,
                       row.marknames = 0, NA.char = "-9", pop = NULL, sep = NULL,
                       ask = F, quiet = FALSE)

#### Summary statistics ####

# This is based on all individuals and all loci
basicstat.all <- basic.stats(all, diploid = TRUE, digits = 2) 
allelic.richness.all <- allelic.richness(all, diploid = TRUE)

#Ar
allelic.richness.all
allelic.richness.all.df <- as.data.frame(allelic.richness.all$Ar)

#get mean Ho and Hx per pop
x.pop = seppop(all) 
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

#### PCA ####

## males
x.males <- tab(males.stru, freq=TRUE, NA.method="mean")
pca.males <- dudi.pca(x.males, center=TRUE, scale=FALSE) #selected 3 axes

s.class(pca.males$li, fac=pop(males.stru), col=funky(15), sub = "PCA males")

# percentages of variation explained
eig.perc.males <- 100*pca.males$eig/sum(pca.males$eig)
head(eig.perc.males) 

## females
x.females <- tab(females.stru, freq=TRUE, NA.method="mean")
pca.females <- dudi.pca(x.females, center=TRUE, scale=FALSE) #selected 3 axes
s.class(pca.females$li, fac=pop(females.stru), col=funky(15), sub = "PCA females") 

# percentages of variation explained
eig.perc.females <- 100*pca.females$eig/sum(pca.females$eig)
head(eig.perc.females) 

## chicks
x.chicks <- tab(chicks.stru, freq=TRUE, NA.method="mean")
pca.chicks <- dudi.pca(x.chicks, center=TRUE, scale=FALSE) #3 axes

s.class(pca.chicks$li, fac=pop(chicks.stru), col=funky(16), sub = "PCA chicks")

# percentages of variation explained
eig.perc.chicks <- 100*pca.chicks$eig/sum(pca.chicks$eig)
head(eig.perc.chicks) 

#### Calculate Fst ####

## Males 
#convert to hfstat object
males.hfstat <- genind2hierfstat(males.stru)
#calculate stats
basicstat.males <- basic.stats(males.stru, diploid = TRUE, digits = 2) # Hobs, mean gene diversities, Fis and Fst

# per locus
fst.males.perlocus <- basicstat.males$perloc$Fst
fst.males.perlocus <- data.frame(Locus = seq(from = 1, to = 12), Fst = fst.males.perlocus)

# Pairwise Fst
fst.males <- pairwise.neifst(males.hfstat)
head(fst.males)

# Fst per population
boxplot(fst.males, col=funky(nPop(males.stru)), las=3,
        xlab="Population", ylab="Fst", main = "Pairwise Fst values per population only males") #pop 6 only has 1 sample

#Bootstrap
boot.fst.males <- boot.ppfst(males.hfstat, nboot = 1000) #confidence interval for Fst, upper triangle is upper limit, lower triangle lower limit

#create a long dataframe for pairwise Fst
boot.fst.males.UL <- boot.fst.males$ul
boot.fst.males.LL <- boot.fst.males$ll

flat.matrix <- function(d){
  data.frame(i=rep(row.names(d),ncol(d)),
             j=rep(colnames(d),each=nrow(d)),
             score=as.vector(d))
}

fst.males.flat <- flat.matrix(fst.males)
names(fst.males.flat) <- c("site.x", "site.y", "Fst")

boot.fst.males.LL.flat <- flat.matrix(boot.fst.males.LL)
names(boot.fst.males.LL.flat) <- c("site.x", "site.y", "LL")

boot.fst.males.UL.flat <- flat.matrix(boot.fst.males.UL)
names(boot.fst.males.UL.flat) <- c("site.x", "site.y", "UL")

pairwise.fst.males <- left_join(fst.males.flat, boot.fst.males.LL.flat, by = c("site.x", "site.y"))
pairwise.fst.males <- left_join(pairwise.fst.males, boot.fst.males.UL.flat, by = c("site.x", "site.y"))

pairwise.fst.males <- subset(pairwise.fst.males, site.x != "-9" & site.y != "-9")
pairwise.fst.males <- subset(pairwise.fst.males, !is.na(Fst) & !is.na(UL))
pairwise.fst.males <- pairwise.fst.males %>% mutate(Significance = case_when(
  UL > 0 & LL > 0 ~ "significant",
  UL > 0 & LL < 0 ~ "insignificant",
  UL < 0 & LL < 0 ~ "significant" ))


## Females 
#convert to hfstat object
females.hfstat <- genind2hierfstat(females.stru)
#calculate stats
basicstat.females <- basic.stats(females.stru, diploid = TRUE, digits = 2) # Hobs, mean gene diversities, Fis and Fst

# per locus
fst.females.perlocus <- basicstat.females$perloc$Fst
fst.females.perlocus <- data.frame(Locus = seq(from = 1, to = 12), Fst = fst.females.perlocus)

# Pairwise Fst
fst.females <- pairwise.neifst(females.hfstat)
head(fst.females)

# Fst per population
boxplot(fst.females, col=funky(nPop(females.stru)), las=3,
        xlab="Population", ylab="Fst", main = "Pairwise Fst values per population females") #pop 6 only has 1 sample

#Bootstrap
boot.fst.females <- boot.ppfst(females.hfstat, nboot = 1000) #confidence interval for Fst, upper triangle is upper limit, lower triangle lower limit

#create long dataframe
boot.fst.females.UL <- boot.fst.females$ul
boot.fst.females.LL <- boot.fst.females$ll

fst.females.flat <- flat.matrix(fst.females)
names(fst.females.flat) <- c("site.x", "site.y", "Fst")

boot.fst.females.LL.flat <- flat.matrix(boot.fst.females.LL)
names(boot.fst.females.LL.flat) <- c("site.x", "site.y", "LL")

boot.fst.females.UL.flat <- flat.matrix(boot.fst.females.UL)
names(boot.fst.females.UL.flat) <- c("site.x", "site.y", "UL")

pairwise.fst.females <- left_join(fst.females.flat, boot.fst.females.LL.flat, by = c("site.x", "site.y"))
pairwise.fst.females <- left_join(pairwise.fst.females, boot.fst.females.UL.flat, by = c("site.x", "site.y"))

pairwise.fst.females <- subset(pairwise.fst.females, site.x != "-9" & site.y != "-9")
pairwise.fst.females <- subset(pairwise.fst.females, !is.na(Fst) & !is.na(UL))
pairwise.fst.females <- pairwise.fst.females %>% mutate(Significance = case_when(
  UL > 0 & LL > 0 ~ "significant",
  UL > 0 & LL < 0 ~ "insignificant",
  UL < 0 & LL < 0 ~ "significant" ))

## Chicks
chicks.hfstat <- genind2hierfstat(chicks.stru)
#calculate stats
basicstat.chicks <- basic.stats(chicks.stru, diploid = TRUE, digits = 2) # Hobs, mean gene diversities, Fis and Fst
# per locus
fst.chicks.perlocus <- basicstat.chicks$perloc$Fst
fst.chicks.perlocus <- data.frame(Locus = seq(from = 1, to = 11), Fst = fst.chicks.perlocus)

## Pairwise Fst
fst.chicks <- pairwise.neifst(chicks.hfstat)

#Fst per population for chicks
boxplot(fst.chicks, col=funky(nPop(chicks.stru)), las=3,
        xlab="Population", ylab="Fst")

# Bootstrap
boot.fst.chicks <- boot.ppfst(chicks.hfstat, nboot = 1000) #confidence interval for Fst, upper triangle is upper limit, lower triangle lower limit
boot.fst.chicks.UL <- boot.fst.chicks$ul
boot.fst.chicks.LL <- boot.fst.chicks$ll

#create long df
fst.chicks.flat <- flat.matrix(fst.chicks)
names(fst.chicks.flat) <- c("site.x", "site.y", "Fst")

boot.fst.chicks.LL.flat <- flat.matrix(boot.fst.chicks.LL)
names(boot.fst.chicks.LL.flat) <- c("site.x", "site.y", "LL")

boot.fst.chicks.UL.flat <- flat.matrix(boot.fst.chicks.UL)
names(boot.fst.chicks.UL.flat) <- c("site.x", "site.y", "UL")

pairwise.fst.chicks <- left_join(fst.chicks.flat, boot.fst.chicks.LL.flat, by = c("site.x", "site.y"))
pairwise.fst.chicks <- left_join(pairwise.fst.chicks, boot.fst.chicks.UL.flat, by = c("site.x", "site.y"))

pairwise.fst.chicks <- subset(pairwise.fst.chicks, site.x != "-9" & site.y != "-9")
pairwise.fst.chicks <- subset(pairwise.fst.chicks, !is.na(Fst) & !is.na(UL))
pairwise.fst.chicks <- pairwise.fst.chicks %>% mutate(Significance = case_when(
  UL > 0 & LL > 0 ~ "significant",
  UL > 0 & LL < 0 ~ "insignificant",
  UL < 0 & LL < 0 ~ "significant" ))
