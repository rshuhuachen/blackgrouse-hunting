#________________________________________________________________________
####### POPULATION STRUCTURE IN HUNTED VS UNHUNTED LEKKING SITES #######
#________________________________________________________________________

##### In this script, we will import the structure data files for adults and chicks 
##### separately and plot some basic distribution and statistics on the microsatellite data,        
##### plot PCA's, test for HWE 

### Load packages ###
install.packages(pacman)
pacman::p_load(dplyr, tibble, adegenet, pegas, data.table, hierfstat)

all <- read.structure("data/rawdata/Microsat.adults.plus.unrelated.chicks.forstructure.stru", n.ind = 2078, n.loc = 14, onerowperind = F,
                         col.lab = 1, col.pop = 2, col.others = NULL,
                         row.marknames = 0, NA.char = "-9", pop = NULL, sep = NULL,
                         ask = F, quiet = FALSE) #just need this for HWE


#### Summary ####
summary(all)

#### Testing for Hardy-Weinberg equilibrium ####

# First all together
allHWE.all <- pegas::hw.test(all, B = 1000) #B = 1000 for 1000 Monte Carlo permutations 
allHWE.all 

# Then per population using a loop
allpop <- seppop(all)
# Run loop
allHWE = NULL
for(i in 1:length(allpop)) {
  hwt <- pegas::hw.test(allpop[[i]], B=1000)
  smry <- summary(allpop[[i]])
  
  Hobs <- smry[[6]]
  Hexp <- smry[[7]]
  pexact <- hwt[,4] #hw.test does chi2 test and exact test. We use p-values of exact test which are given in 4th col
  qval.FDR <- p.adjust(pexact, method = "fdr")
  qval.bon <- p.adjust(pexact, method = "bonferroni")
  allHWE <- as.data.frame(cbind(allHWE, Hobs, Hexp, pexact, qval.FDR, qval.bon))
  
}

sites<-rep(names(allpop[1:length(allpop)]),each=5)
allHWE <- rbind(allHWE, sites)
allHWE <- allHWE[c(nrow(allHWE),1:(nrow(allHWE)-1)),]
rownames(allHWE)[1] <- "Site"

allHWE.t <- as.data.frame(t(allHWE))
nums <- c(2:15)
allHWE.t[nums] <- lapply(allHWE.t[nums], as.numeric)

allHWE.t[,c(2:15)]<- round(allHWE.t[,c(2:15)], 2)
head(allHWE.t)
View(allHWE.t[c(seq(from=4, by = 5, to = 59)),])

write.table(allHWE.t, "tables/HWE_all.csv", quote=F)

# locus names: BG10, BG16, BG18, BG15, BG19, BG6, TTT1, TTD2, TTD3, TUD6, TUT3a, TUT4, BG20, TTT2
#exclude locus 1, and 13 as FDR-corrected q-values indicate violation
#of assumption of HWE in over 70% of the sites

raw <- read.table("data/rawdata/Microsat.adults.plus.unrelated.chicks.forstructure.stru") #total of 14 loci, of which the last one (L14, TTT2 only sampled in adults)
data_noL1 <- raw[,-c(3)] #this locus resulted in a very odd PCA plot in female adults only so we excluded it from the entire paper
data_noL13 <- raw[,-c(15)]
data_noL14 <- raw[,-c(16)]
data_noL1L13 <- raw[,-c(3,15)]
data_noL1L14 <- raw[,-c(3,16)]
data_noL13L14 <- raw[,-c(15,16)]
data_noL1L13L14 <-   raw[,-c(3,15,16)]

write.table(data_noL1, "data/cleandata/Microsat.adults.plus.unrelated.chicks.noLOCUS1.forstructure.stru", col.names = F, row.names = F, quote=F, sep = " ")
write.table(data_noL13, "data/cleandata/Microsat.adults.plus.unrelated.chicks.noLOCUS13.forstructure.stru", col.names = F, row.names = F, quote=F, sep = " ")
write.table(data_noL14, "data/cleandata/Microsat.adults.plus.unrelated.chicks.noLOCUS14.forstructure.stru", col.names = F, row.names = F, quote=F, sep = " ")
write.table(data_noL1L13, "data/cleandata/Microsat.adults.plus.unrelated.chicks.noLOCUS1+13.forstructure.stru", col.names = F, row.names = F, quote=F, sep = " ")
write.table(data_noL1L14, "data/cleandata/Microsat.adults.plus.unrelated.chicks.noLOCUS1+14.forstructure.stru", col.names = F, row.names = F, quote=F, sep = " ")
write.table(data_noL13L14, "data/cleandata/Microsat.adults.plus.unrelated.chicks.noLOCUS13+14.forstructure.stru", col.names = F, row.names = F, quote=F, sep = " ")
write.table(data_noL1L13L14, "data/cleandata/Microsat.adults.plus.unrelated.chicks.noLOCUS1+13+14.forstructure.stru", col.names = F, row.names = F, quote=F, sep = " ")

### Explore PCA with different HWE filters and calculate Fis
data_stru <- read.structure("data/rawdata/Microsat.adults.plus.unrelated.chicks.forstructure.stru", n.ind = 2078, n.loc = 14, onerowperind = F,
                                 col.lab = 1, col.pop = 2, col.others = NULL,
                                 row.marknames = 0, NA.char = "-9", pop = NULL, sep = NULL,
                                 ask = F, quiet = FALSE)

data_noL1_stru <- read.structure("data/cleandata/Microsat.adults.plus.unrelated.chicks.noLOCUS1.forstructure.stru", n.ind = 2078, n.loc = 13, onerowperind = F,
                             col.lab = 1, col.pop = 2, col.others = NULL,
                             row.marknames = 0, NA.char = "-9", pop = NULL, sep = NULL,
                             ask = F, quiet = FALSE)

data_noL13_stru <- read.structure("data/cleandata/Microsat.adults.plus.unrelated.chicks.noLOCUS13.forstructure.stru", n.ind = 2078, n.loc = 13, onerowperind = F,
                                 col.lab = 1, col.pop = 2, col.others = NULL,
                                 row.marknames = 0, NA.char = "-9", pop = NULL, sep = NULL,
                                 ask = F, quiet = FALSE)
data_noL14_stru <- read.structure("data/cleandata/Microsat.adults.plus.unrelated.chicks.noLOCUS14.forstructure.stru", n.ind = 2078, n.loc = 13, onerowperind = F,
                                 col.lab = 1, col.pop = 2, col.others = NULL,
                                 row.marknames = 0, NA.char = "-9", pop = NULL, sep = NULL,
                                 ask = F, quiet = FALSE)
data_noL1L13_stru <- read.structure("data/cleandata/Microsat.adults.plus.unrelated.chicks.noLOCUS1+13.forstructure.stru", n.ind = 2078, n.loc = 12, onerowperind = F,
                                 col.lab = 1, col.pop = 2, col.others = NULL,
                                 row.marknames = 0, NA.char = "-9", pop = NULL, sep = NULL,
                                 ask = F, quiet = FALSE)
data_noL1L14_stru <- read.structure("data/cleandata/Microsat.adults.plus.unrelated.chicks.noLOCUS1+14.forstructure.stru", n.ind = 2078, n.loc = 12, onerowperind = F,
                                 col.lab = 1, col.pop = 2, col.others = NULL,
                                 row.marknames = 0, NA.char = "-9", pop = NULL, sep = NULL,
                                 ask = F, quiet = FALSE)
data_noL13L14_stru <- read.structure("data/cleandata/Microsat.adults.plus.unrelated.chicks.noLOCUS13+14.forstructure.stru", n.ind = 2078, n.loc = 12, onerowperind = F,
                                 col.lab = 1, col.pop = 2, col.others = NULL,
                                 row.marknames = 0, NA.char = "-9", pop = NULL, sep = NULL,
                                 ask = F, quiet = FALSE)
data_noL1L13L14_stru <- read.structure("data/cleandata/Microsat.adults.plus.unrelated.chicks.noLOCUS1+13+14.forstructure.stru", n.ind = 2078, n.loc = 11, onerowperind = F,
                                 col.lab = 1, col.pop = 2, col.others = NULL,
                                 row.marknames = 0, NA.char = "-9", pop = NULL, sep = NULL,
                                 ask = F, quiet = FALSE)

#PCA plots
basic.stats(data_stru)
x <- tab(data_stru, freq=TRUE, NA.method="mean")
pca <- dudi.pca(x, center=TRUE, scale=FALSE) #selected 3 axes
s.class(pca$li, fac=pop(data_stru), col=funky(15), sub = "PCA all")

x_noL1 <- tab(data_noL1_stru, freq=TRUE, NA.method="mean")
pca.noL1 <- dudi.pca(x_noL1, center=TRUE, scale=FALSE) #selected 3 axes
s.class(pca.noL1$li, fac=pop(data_noL1_stru), col=funky(15), sub = "PCA all no L1")

x_noL13 <- tab(data_noL13_stru, freq=TRUE, NA.method="mean")
pca.noL13 <- dudi.pca(x_noL13, center=TRUE, scale=FALSE) #selected 3 axes
s.class(pca.noL13$li, fac=pop(data_noL13_stru), col=funky(15), sub = "PCA all no L13")

x_noL14 <- tab(data_noL14_stru, freq=TRUE, NA.method="mean")
pca.noL14 <- dudi.pca(x_noL14, center=TRUE, scale=FALSE) #selected 3 axes
s.class(pca.noL14$li, fac=pop(data_noL14_stru), col=funky(15), sub = "PCA all no L14")

x_noL1L13 <- tab(data_noL1L13_stru, freq=TRUE, NA.method="mean")
pca.noL1L13 <- dudi.pca(x_noL1L13, center=TRUE, scale=FALSE) #selected 3 axes
s.class(pca.noL1L13$li, fac=pop(data_noL1L13_stru), col=funky(15), sub = "PCA all no L1 or L13")

x_noL1L14 <- tab(data_noL1L14_stru, freq=TRUE, NA.method="mean")
pca.noL1L14 <- dudi.pca(x_noL1L14, center=TRUE, scale=FALSE) #selected 3 axes
s.class(pca.noL1L14$li, fac=pop(data_noL1L14_stru), col=funky(15), sub = "PCA all no L1 or L14")

x_noL13L14 <- tab(data_noL13L14_stru, freq=TRUE, NA.method="mean")
pca.noL13L14 <- dudi.pca(x_noL13L14, center=TRUE, scale=FALSE) #selected 3 axes
s.class(pca.noL13L14$li, fac=pop(data_noL13L14_stru), col=funky(15), sub = "PCA all no L13 or L14")

basic.stats(data_noL1L13L14_stru)
x_noL1L13L14 <- tab(data_noL1L13L14_stru, freq=TRUE, NA.method="mean")
pca.noL1L13L14 <- dudi.pca(x_noL1L13L14, center=TRUE, scale=FALSE) #selected 3 axes
s.class(pca.noL1L13L14$li, fac=pop(data_noL1L13L14_stru), col=funky(15), sub = "PCA all no L1 or L13 or L14")

