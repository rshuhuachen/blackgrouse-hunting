#### Population Structure #####

# In this script, we will 1) calculate genetic summary statistics, 2) make a PCA, 
# 3) calculate pairwise Fst

## load libraries
library(data.table); library(tidyverse); library(hierfstat); 
library(plot.matrix); library(lme4); library(adegenet); library(forcats)
library(ape)

all.raw <- read.structure ("data/cleandata/Microsat.adults.plus.unrelated.chicks.noLOCUS1.forstructure.stru", n.ind = 2078, n.loc = 13, onerowperind = F,
                       col.lab = 1, col.pop = 2, col.others = NULL,
                       row.marknames = 0, NA.char = "-9", pop = NULL, sep = NULL,
                       ask = F, quiet = FALSE)

all <- read.structure ("data/cleandata/Microsat.adults.plus.unrelated.chicks.noLOCUS1+13.forstructure.stru", n.ind = 2078, n.loc = 12, onerowperind = F,
                                col.lab = 1, col.pop = 2, col.others = NULL,
                                row.marknames = 0, NA.char = "-9", pop = NULL, sep = NULL,
                                ask = F, quiet = FALSE)

pops <- read.csv("data/details/Codes.pops.both.filtered_withcoord.csv")

#### Summary statistics ####

# This is based on all individuals and all loci
basicstat.all <- basic.stats(all.raw, diploid = TRUE, digits = 2) 
allelic.richness.all <- allelic.richness(all.raw, diploid = TRUE)

#Ar
allelic.richness.all
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

Ho_model <- lmerTest::lmer(Hobs ~ hunt + (1|locus) + (1|pop), data = subset(measures.all, locus != "L13")) # excluding BG20
Ho_model_null <- lmerTest::lmer(Hobs ~ (1|locus) + (1|pop), data = subset(measures.all, locus != "L13")) # excluding BG20

He_model <- lmerTest::lmer(Hexp ~ hunt + (1|locus) + (1|pop), data = subset(measures.all, locus != "L13"))
He_model_null <- lmerTest::lmer(Hexp ~ (1|locus) + (1|pop), data = subset(measures.all, locus != "L13"))

Ar_model <- lmerTest::lmer(Ar ~ hunt + (1|locus) + (1|pop), data = subset(measures.all, locus != "L13"))
Ar_model_null <- lmerTest::lmer(Ar ~ (1|locus) + (1|pop), data = subset(measures.all, locus != "L13"))

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

write.csv(pairwise.fst.all, "data/tables/Pairwise_Fst_all.csv", row.names = F, quote=F)

